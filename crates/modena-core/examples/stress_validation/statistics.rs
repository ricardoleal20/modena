//! Inference on the second half only. Eight seed labels are paired across timesteps.
use super::protocol::*;
use chacha20::ChaCha12Rng;
use rand::{RngExt, SeedableRng};
use serde::Serialize;
pub const T7: f64 = 2.3646242510102993;
pub fn mean(v: &[f64]) -> R<f64> {
    finite(v)?;
    if v.is_empty() {
        return Err("empty series".into());
    }
    let z = v.iter().sum::<f64>() / v.len() as f64;
    finite(&[z])?;
    Ok(z)
}
pub fn replica_se(x: &[f64]) -> R<f64> {
    let mu = mean(x)?;
    if x.len() != 8 {
        return Err("eight seeds required".into());
    }
    let se = (x.iter().map(|v| (v - mu).powi(2)).sum::<f64>() / 56.).sqrt();
    finite(&[se])?;
    Ok(se)
}
/// Geyer's initial positive pairs: tau_int = -1/2 + sum (rho[2k]+rho[2k+1]).
pub fn tau(v: &[f64]) -> R<f64> {
    let m = mean(v)?;
    let x: Vec<_> = v.iter().map(|a| a - m).collect();
    finite(&x)?;
    if x.len() < 4 {
        return Err("too few autocorrelation samples".into());
    }
    let c0 = x.iter().map(|a| a * a).sum::<f64>();
    finite(&[c0])?;
    if c0 == 0. {
        return Ok(0.5);
    }
    let cov = |d: usize| x.iter().zip(&x[d..]).map(|(a, b)| a * b).sum::<f64>() / c0;
    let mut t = -0.5;
    for lag in (0..x.len() - 1).step_by(2) {
        let pair = cov(lag) + cov(lag + 1);
        finite(&[pair])?;
        if pair <= 0. {
            break;
        }
        t += pair
    }
    finite(&[t])?;
    Ok(t.max(0.5))
}
pub fn first_half_tau(x: &[Vec<f64>]) -> R<f64> {
    if x.len() != 8 || x[0].len() < 8 || x.iter().any(|r| r.len() != x[0].len()) {
        return Err("expected eight equally long replicas".into());
    }
    let mut t = 0.5_f64;
    for r in x {
        finite(r)?;
        t = t.max(tau(&r[..r.len() / 2])?);
    }
    Ok(t)
}
#[derive(Debug, Serialize)]
pub struct Inference {
    pub estimate: Estimate,
    pub seed_means: Vec<f64>,
    pub bootstrap_ci95: Option<[f64; 2]>,
    pub tau_samples: f64,
    pub enough: bool,
}
pub fn infer(rows: &[Vec<f64>], block: usize, t: f64) -> R<Inference> {
    if rows.len() != 8 || block == 0 || !t.is_finite() || t <= 0. {
        return Err("invalid blocking layout".into());
    }
    let mut blocks = vec![];
    let mut means = vec![];
    for row in rows {
        finite(row)?;
        if row.len() != rows[0].len() {
            return Err("unequal replicas".into());
        }
        let b = row[row.len() / 2..]
            .chunks_exact(block)
            .map(mean)
            .collect::<R<Vec<_>>>()?;
        // A retained half shorter than a single block still supplies the observed seed mean.
        means.push(if b.is_empty() {
            mean(&row[row.len() / 2..])?
        } else {
            mean(&b)?
        });
        blocks.push(b);
    }
    let nb = blocks[0].len();
    let mu = mean(&means)?;
    let se = replica_se(&means)?;
    let ci = [mu - T7 * se, mu + T7 * se];
    finite(&ci)?;
    let est = Estimate {
        mean: mu,
        se,
        ci95: ci,
        n_eff: 8. * (nb * block) as f64 / (2. * t),
        block_samples: block,
        retained_blocks: nb * 8,
    };
    let bootstrap = if nb > 0 {
        let mut rng = ChaCha12Rng::seed_from_u64(7001);
        let mut draws = Vec::with_capacity(2000);
        // Independent resampling of contiguous blocks within each independent seed.
        for _ in 0..2000 {
            let mut total = 0.;
            for b in &blocks {
                for _ in 0..nb {
                    total += b[rng.random_range(0..nb)]
                }
            }
            draws.push(total / (8 * nb) as f64);
        }
        finite(&draws)?;
        draws.sort_by(f64::total_cmp);
        Some([draws[49], draws[1949]])
    } else {
        None
    };
    let enough = nb >= 20 && est.n_eff >= 1000.;
    Ok(Inference {
        estimate: est,
        seed_means: means,
        bootstrap_ci95: bootstrap,
        tau_samples: t,
        enough,
    })
}
pub fn estimate_replicates(samples: &[Vec<f64>], sample_interval: f64) -> R<Estimate> {
    if sample_interval != 0.1 {
        return Err("expected .1 nominal sample clock".into());
    }
    let t = first_half_tau(samples)?;
    Ok(infer(samples, (10. * t).ceil() as usize, t)?.estimate)
}
#[derive(Serialize)]
pub struct Fit {
    pub intercept: f64,
    pub se: f64,
    pub ci95: [f64; 2],
    pub slope: f64,
    pub weighted_lack_of_fit: f64,
    pub residual_df: usize,
    pub paired_intercepts: [f64; 8],
    pub drop_shift: f64,
    pub drop_paired_se: f64,
    pub drop_coarsest_intercept_se: f64,
    pub bias_ratio_orders: Vec<[f64; 3]>,
    pub weak_order_resolved: bool,
}
fn regression(rows: &[Inference], dt: &[f64]) -> R<(f64, f64, [f64; 8], f64)> {
    finite(dt)?;
    if rows.len() != dt.len()
        || rows.len() < 3
        || rows
            .iter()
            .any(|r| r.seed_means.len() != 8 || !r.estimate.mean.is_finite())
    {
        return Err("h² regression roster".into());
    }
    let mut sw = 0.;
    let mut sx = 0.;
    let mut sxx = 0.;
    for (r, h) in rows.iter().zip(dt) {
        if !r.estimate.se.is_finite() || r.estimate.se <= 0. {
            return Err("unresolved between-seed uncertainty".into());
        }
        let w = 1. / r.estimate.se.powi(2);
        sw += w;
        sx += w * h * h;
        sxx += w * h.powi(4);
    }
    let det = sw * sxx - sx * sx;
    finite(&[det])?;
    if det <= 0. {
        return Err("degenerate fit".into());
    }
    let mut paired = [0.; 8];
    let mut slope = 0.;
    for (r, h) in rows.iter().zip(dt) {
        let x = h * h;
        let w = 1. / r.estimate.se.powi(2);
        for (y, mu) in paired.iter_mut().zip(&r.seed_means) {
            *y += w * (sxx - sx * x) / det * mu;
        }
        slope += w * (sw * x - sx) / det * r.estimate.mean;
    }
    let a = mean(&paired)?;
    let lack = rows
        .iter()
        .zip(dt)
        .map(|(r, h)| ((r.estimate.mean - a - slope * h * h) / r.estimate.se).powi(2))
        .sum();
    finite(&[slope, lack])?;
    Ok((a, slope, paired, lack))
}
pub fn fit(rows: &[Inference], dt: &[f64], exact: f64) -> R<Fit> {
    if dt != [0.04, 0.02, 0.01, 0.005] || !exact.is_finite() {
        return Err("unregistered fit".into());
    }
    let (a, b, p, chi) = regression(rows, dt)?;
    let (a1, _, p1, _) = regression(&rows[1..], &dt[1..])?;
    let se = replica_se(&p)?;
    let se1 = replica_se(&p1)?;
    let d: Vec<_> = p.iter().zip(p1).map(|(a, b)| a - b).collect();
    let mut orders = vec![];
    for r in rows.windows(2) {
        let u = r[0].estimate.mean - exact;
        let v = r[1].estimate.mean - exact;
        if u * v <= 0. || u.abs() <= 2. * r[0].estimate.se || v.abs() <= 2. * r[1].estimate.se {
            continue;
        }
        let power = (u / v).log2();
        let influence: Vec<_> = r[0]
            .seed_means
            .iter()
            .zip(&r[1].seed_means)
            .map(|(a, b)| {
                ((a - r[0].estimate.mean) / u - (b - r[1].estimate.mean) / v)
                    / std::f64::consts::LN_2
            })
            .collect();
        let half = T7 * replica_se(&influence)?;
        let o = [power, power - half, power + half];
        finite(&o)?;
        orders.push(o);
    }
    let resolved = rows.iter().all(|r| r.enough)
        && !orders.is_empty()
        && orders.iter().all(|r| r[1] >= 1. && r[2] <= 3.);
    let ci = [a - T7 * se, a + T7 * se];
    finite(&ci)?;
    finite(&[a - a1, se.hypot(se1)])?;
    Ok(Fit {
        intercept: a,
        se,
        ci95: ci,
        slope: b,
        weighted_lack_of_fit: chi,
        residual_df: rows.len() - 2,
        paired_intercepts: p,
        drop_shift: a - a1,
        drop_paired_se: replica_se(&d)?,
        drop_coarsest_intercept_se: se1,
        bias_ratio_orders: orders,
        weak_order_resolved: resolved,
    })
}

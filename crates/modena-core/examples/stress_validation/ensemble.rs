//! Disjoint production BAOAB trajectories, with continued physical/RNG states across duration rungs.
use super::{
    fixtures::{common, dyad, inertia, invert, relative, rotate, Zero, CENTER},
    protocol::*,
    statistics::*,
};
use chacha20::ChaCha12Rng;
use modena_core::constraints::RigidWater;
use modena_core::stress::{ConstraintSolverOptions, ModelEvaluation};
use modena_core::{
    Constraint, Constraints, Integrator, LangevinBAOAB, NeighborList, SimBox, System, Vec3,
};
use rand::SeedableRng;
use rand_distr::{Distribution, StandardNormal};
use std::{
    fs::{File, OpenOptions},
    io::{BufRead, BufReader, BufWriter, Write},
    path::Path,
};
const DT: [f64; 4] = [0.04, 0.02, 0.01, 0.005];
const IDEAL: f64 = 1. / 54.;
fn normal(r: &mut ChaCha12Rng) -> Vec3 {
    Vec3(std::array::from_fn(|_| StandardNormal.sample(r)))
}
pub fn gas_fixture(water: bool, charged: bool, stream: u64) -> R<ValidationFixture> {
    let (base, _) = if water {
        common::water([15.9994, 1.008, 1.008])
    } else {
        common::rotor(1., 3.)
    };
    let n = base.n();
    let m = &base.masses;
    let total = m.iter().sum::<f64>();
    let mut s = System::new(SimBox::cubic(12.));
    s.resize(32 * n);
    let mut constraints = Constraints::new();
    let mut orientation = ChaCha12Rng::seed_from_u64(901);
    let mut maxwell = ChaCha12Rng::seed_from_u64(902);
    maxwell.set_stream(stream);
    let mut scales = vec![];
    for j in 0..32 {
        let p = Vec3::new(
            1.5 + 3. * (j % 4) as f64,
            1.5 + 3. * ((j / 4) % 4) as f64,
            3. + 6. * (j / 16) as f64,
        );
        let mut q: [f64; 4] = std::array::from_fn(|_| StandardNormal.sample(&mut orientation));
        let norm = q.iter().map(|x| x * x).sum::<f64>().sqrt();
        q = q.map(|a| a / norm);
        let r: Vec<_> = (0..n).map(|i| rotate(q, base.pos(i) - CENTER)).collect();
        let vcm = normal(&mut maxwell) * (1. / total.sqrt());
        // Independent Maxwell atom momenta yield body angular momentum with covariance I.
        let l = (0..n).fold(Vec3::zero(), |a, i| {
            a + r[i].cross(normal(&mut maxwell)) * m[i].sqrt()
        });
        let w = if water {
            let axes = [
                Vec3::new(1., 0., 0.),
                Vec3::new(0., 1., 0.),
                Vec3::new(0., 0., 1.),
            ];
            invert(axes.map(|a| inertia(&r, m, a)), l)?
        } else {
            l * (1. / (m[0] * m[1] / total * 1.3_f64.powi(2)))
        };
        let start = j * n;
        for i in 0..n {
            s.set_pos(start + i, p + r[i]);
            s.set_vel(start + i, vcm + w.cross(r[i]));
            s.masses[start + i] = m[i];
            if charged {
                s.charges[start + i] = [-0.8, 0.4, 0.4][i];
            }
        }
        if water {
            constraints.waters.push(RigidWater {
                first: start,
                d_oh: 1.,
                d_hh: 2. * (109.47_f64.to_radians() / 2.).sin(),
            });
            if charged {
                for (i, j) in [(0, 1), (0, 2), (1, 2)] {
                    scales.push((start + i, start + j, 0.));
                }
            }
        } else {
            constraints.bonds.push(Constraint {
                i: start,
                j: start + 1,
                d: 1.3,
            });
        }
    }
    let id = if charged {
        if stream == 0 {
            "interacting32".into()
        } else {
            format!("interacting32-{stream}")
        }
    } else if water {
        "gas-water".into()
    } else {
        "gas-rotor".into()
    };
    let mut f = ValidationFixture::from_system(id, &s, constraints);
    if charged {
        f.kappa = Some(0.8);
        f.real_cutoff = Some(3.5);
        f.reciprocal_cutoff = Some(8.);
        f.pair_scales = scales;
    }
    Ok(f)
}
pub fn sample(
    s: &System,
    n: usize,
    seed: u64,
    step: u64,
    dt: f64,
    endpoint: [f64; 6],
    interval: [f64; 6],
) -> R<Sample> {
    let mut com_energy = 0.;
    let mut rotation_energy = 0.;
    let mut orient = [0.; 6];
    for start in (0..s.n()).step_by(n) {
        let (r, v, vcm) = relative(s, start, n);
        let m = &s.masses[start..start + n];
        com_energy += 0.5 * m.iter().sum::<f64>() * vcm.norm2();
        rotation_energy += 0.5 * v.iter().zip(m).map(|(v, m)| v.norm2() * m).sum::<f64>();
        let axis = if n == 2 {
            r[1] - r[0]
        } else {
            (r[1] - r[0]).cross(r[2] - r[0])
        };
        let e = axis * (1. / axis.norm());
        for (x, a) in orient.iter_mut().zip(dyad(e, e)) {
            *x += a / (s.n() / n) as f64
        }
    }
    let v = Sample {
        seed,
        step,
        time: step as f64 * dt,
        dt,
        endpoint,
        interval,
        com_energy,
        rotation_energy,
        orientation_second: orient,
    };
    v.validate()?;
    Ok(v)
}
pub fn sample_step(i: usize, k: u64) -> u64 {
    match i {
        0 => (5 * k).div_ceil(2),
        1 => 5 * k,
        2 => 10 * k,
        3 => 20 * k,
        _ => unreachable!(),
    }
}
struct Replica {
    s: System,
    integrator: LangevinBAOAB,
    nl: NeighborList,
    step: u64,
    seed: u64,
    n: usize,
}
fn extend(
    out: &Path,
    kind: &str,
    i: usize,
    r: &mut Replica,
    duration: u64,
    prior: u64,
    work: &mut u64,
) -> R<()> {
    let dt = DT[i];
    let burn = (50. / dt).round() as u64;
    let target = burn + sample_step(i, 10 * duration);
    let path = out.join(format!("{kind}-{i}-{}.jsonl", r.seed));
    let mut file = BufWriter::new(
        if prior == 0 {
            File::create_new(path)
        } else {
            OpenOptions::new().append(true).open(path)
        }
        .map_err(err)?,
    );
    let mut k = 10 * prior + 1;
    let mut next = burn + sample_step(i, k);
    while r.step < target {
        let t = r
            .integrator
            .step_with_stress_options(
                &mut r.s,
                &Zero,
                &mut r.nl,
                dt,
                &ModelEvaluation::zero(),
                ConstraintSolverOptions::TIGHT,
            )
            .map_err(err)?
            .stress
            .ok_or("missing BAOAB stress")?;
        r.step += 1;
        *work += 32;
        if r.step == next {
            let z = sample(
                &r.s,
                r.n,
                r.seed,
                r.step,
                dt,
                tensor(&t.endpoint.total)?,
                tensor(&t.interval.total)?,
            )?;
            serde_json::to_writer(&mut file, &z).map_err(err)?;
            writeln!(file).map_err(err)?;
            k += 1;
            next = burn + sample_step(i, k);
        }
        if r.step.is_multiple_of(256) || r.step == target {
            file.flush().map_err(err)?;
            let p = out.join("gas-progress.tmp");
            std::fs::write(&p, format!("{work}\n")).map_err(err)?;
            std::fs::rename(p, out.join("gas-progress")).map_err(err)?;
        }
    }
    let snap = serde_json::json!({"step":r.step,"seed":r.seed,"dt":dt,"q":r.s.positions,"v":r.s.velocities,"forces":r.s.forces,"rng_controller":r.integrator.checkpoint_state().map_err(err)?,"molecule_steps_cumulative":work});
    write_new(
        &out.join(format!(
            "continuation-{kind}-{i}-{}-{duration}.json",
            r.seed
        )),
        &snap,
    )
}
// Even integer degrees of freedom: exact finite-sum survival function of chi-square.
fn chi_quantile(df: usize, p: f64) -> f64 {
    let mut lo = 0.;
    let mut hi = 4. * df as f64;
    for _ in 0..80 {
        let x = (hi + lo) * 0.25;
        let mut term = 1.;
        let mut sum = 1.;
        for n in 1..df / 2 {
            term *= x / n as f64;
            sum += term
        }
        if 1. - (-x).exp() * sum < p {
            lo = 2. * x
        } else {
            hi = 2. * x
        }
    }
    (hi + lo) / 2.
}
#[derive(Clone, Copy)]
enum Obs {
    Pressure(bool, usize),
    Energy(bool, usize),
    Orientation(usize),
}
impl Obs {
    fn parameters(self, water: bool) -> (String, f64, f64, bool) {
        match self {
            Self::Pressure(interval, k) => (
                format!(
                    "{}_{}",
                    if interval { "interval" } else { "endpoint" },
                    ["xx", "yy", "zz", "xy", "xz", "yz", "scalar"][k]
                ),
                if k < 3 || k == 6 { IDEAL } else { 0. },
                0.01 * IDEAL,
                true,
            ),
            Self::Orientation(k) => (
                format!("orientation_{k}"),
                if k < 3 { 1. / 3. } else { 0. },
                0.02,
                false,
            ),
            Self::Energy(rot, i) => {
                let n = if rot && !water { 64. } else { 96. };
                let exact = [n, 2. * n, 0.1, 0.5, 0.9][i];
                (
                    format!(
                        "{}_{}",
                        if rot { "rotation" } else { "com" },
                        ["mean", "variance", "cdf10", "cdf50", "cdf90"][i]
                    ),
                    exact,
                    if i == 0 {
                        0.01 * n
                    } else if i == 1 {
                        0.03 * 2. * n
                    } else {
                        0.02
                    },
                    false,
                )
            }
        }
    }
    fn values(self, rows: &[Sample], water: bool) -> Vec<f64> {
        match self {
            Self::Pressure(interval, k) => rows
                .iter()
                .map(|r| {
                    let t = if interval { r.interval } else { r.endpoint };
                    if k == 6 {
                        t[..3].iter().sum::<f64>() / 3.
                    } else {
                        t[k]
                    }
                })
                .collect(),
            Self::Orientation(k) => rows.iter().map(|r| r.orientation_second[k]).collect(),
            Self::Energy(rot, i) => {
                let q = if i >= 2 {
                    chi_quantile(if rot && !water { 64 } else { 96 }, [0.1, 0.5, 0.9][i - 2])
                } else {
                    0.
                };
                rows.iter()
                    .map(|r| {
                        let x = 2. * if rot { r.rotation_energy } else { r.com_energy };
                        if i >= 2 {
                            f64::from(x <= q)
                        } else {
                            x
                        }
                    })
                    .collect()
            }
        }
    }
}
fn squared_deviations(row: &mut [f64]) -> R<()> {
    let n = row.len();
    if n < 2 {
        return Err("variance requires two samples".into());
    }
    let m = mean(row)?;
    for v in &mut *row {
        *v = (*v - m).powi(2) * n as f64 / (n - 1) as f64
    }
    finite(row)
}
fn analysis(out: &Path, id: &str, duration: u64) -> R<(Vec<GateResult>, bool)> {
    let water = id == "gas-water";
    let mut obs = vec![];
    for interval in [false, true] {
        for k in 0..7 {
            obs.push(Obs::Pressure(interval, k))
        }
    }
    for rot in [false, true] {
        for i in 0..5 {
            obs.push(Obs::Energy(rot, i))
        }
    }
    for k in 0..6 {
        obs.push(Obs::Orientation(k))
    }
    let mut by_obs: Vec<Vec<Inference>> = obs.iter().map(|_| vec![]).collect();
    let mut precise = true;
    let mut gates = vec![];
    for (i, dt) in DT.into_iter().enumerate() {
        let mut data = vec![];
        for seed in 1001..=1008 {
            let path = out.join(format!("{id}-{i}-{seed}.jsonl"));
            let mut series = vec![];
            for (k, line) in BufReader::new(File::open(path).map_err(err)?)
                .lines()
                .enumerate()
            {
                let row: Sample = serde_json::from_str(&line.map_err(err)?).map_err(err)?;
                row.validate()?;
                if row.step != (50. / dt).round() as u64 + sample_step(i, k as u64 + 1)
                    || row.time != row.step as f64 * dt
                    || row.dt != dt
                    || row.seed != seed
                {
                    return Err("gas sample identity/time".into());
                }
                series.push(row);
            }
            if series.len() != 10 * duration as usize {
                return Err("gas sample count".into());
            }
            data.push(series);
        }
        let mut metrics = vec![];
        let mut tt = vec![];
        for &o in &obs {
            let mut cols: Vec<_> = data.iter().map(|s| o.values(s, water)).collect();
            if matches!(o, Obs::Energy(_, 1)) {
                for r in &mut cols {
                    let mid = r.len() / 2;
                    let (a, b) = r.split_at_mut(mid);
                    squared_deviations(a)?;
                    squared_deviations(b)?;
                }
            }
            tt.push(first_half_tau(&cols)?);
            metrics.push(cols);
        }
        let block = (10. * tt.iter().copied().fold(0.5_f64, f64::max)).ceil() as usize;
        for (k, o) in obs.iter().enumerate() {
            let (name, _, _, pressure) = o.parameters(water);
            if matches!(o, Obs::Energy(_, 1)) {
                for (s, r) in data.iter().zip(&mut metrics[k]) {
                    let n = r.len() / 2;
                    let nb = (r.len() - n) / block;
                    let len = if nb == 0 { r.len() - n } else { nb * block };
                    let raw = o.values(&s[n..n + len], water);
                    r[n..n + len].copy_from_slice(&raw);
                    squared_deviations(&mut r[n..n + len])?;
                }
            }
            let est = infer(&metrics[k], block, tt[k])?;
            let ok = est.enough && (!pressure || T7 * est.estimate.se <= 0.005 * IDEAL);
            precise &= ok;
            let mut pre = gate(
                &format!("{id}_{duration}_{dt}_{name}_precision"),
                &[
                    est.estimate.n_eff.min(1000.),
                    (est.estimate.retained_blocks / 8).min(20) as f64,
                ],
                &[1000., 20.],
                &[0., 0.],
            )?;
            if !ok {
                pre.status = CaseStatus::Inconclusive
            }
            gates.push(pre);
            by_obs[k].push(est);
        }
    }
    for (o, rows) in obs.iter().zip(by_obs) {
        let (name, exact, practical, pressure) = o.parameters(water);
        if rows.iter().any(|r| r.estimate.se == 0.) {
            write_new(
                &out.join(format!("stats-{id}-{duration}-{name}.json")),
                &serde_json::json!({"dt":DT,"rows":rows,"fit":null,"reason":"zero between-seed weight variance"}),
            )?;
            let mut g = gate(
                &format!("{id}_{duration}_{name}_unresolved_fit"),
                &[rows[3].estimate.mean],
                &[exact],
                &[practical],
            )?;
            g.status = CaseStatus::Inconclusive;
            gates.push(g);
            precise = false;
            continue;
        }
        let fit = fit(&rows, &DT, exact)?;
        write_new(
            &out.join(format!("stats-{id}-{duration}-{name}.json")),
            &serde_json::json!({"dt":DT,"rows":rows,"fit":fit,"target":exact,"practical":practical}),
        )?;
        let mut g = gate(
            &format!("{id}_{duration}_{name}_intercept"),
            &[fit.intercept],
            &[exact],
            &[practical.min(3. * fit.se)],
        )?;
        let enough = rows.iter().all(|r| r.enough);
        let ci_ok = !pressure || T7 * fit.se <= 0.005 * IDEAL;
        precise &= ci_ok;
        if !enough || !ci_ok {
            g.status = CaseStatus::Inconclusive
        }
        gates.push(g);
        if pressure {
            let last = &rows[3].estimate;
            let mut fine = gate(
                &format!("{id}_{duration}_{name}_finest"),
                &[last.mean],
                &[exact],
                &[(3. * last.se).min(practical)],
            )?;
            let mut drop = gate(
                &format!("{id}_{duration}_{name}_drop_coarse"),
                &[fit.drop_shift],
                &[0.],
                &[(2. * fit.drop_paired_se).min(0.005 * IDEAL)],
            )?;
            if !enough || T7 * last.se > 0.005 * IDEAL {
                fine.status = CaseStatus::Inconclusive
            }
            if !enough {
                drop.status = CaseStatus::Inconclusive
            }
            gates.extend([fine, drop]);
        }
    }
    Ok((gates, precise))
}
pub fn run_ensemble(out: &Path) -> R<Vec<GateResult>> {
    let p: ValidationProtocol = read(&out.join("protocol.json"))?;
    p.validate()?;
    let mut work = 0;
    let mut gates = vec![];
    let mut groups = vec![];
    for (id, n) in [("gas-rotor", 2), ("gas-water", 3)] {
        let f = p
            .fixtures
            .iter()
            .find(|f| f.id == id)
            .ok_or("missing gas fixture")?;
        let mut reps = vec![];
        for _ in DT {
            let mut row = vec![];
            for seed in 1001..=1008 {
                row.push(Replica {
                    s: f.system()?,
                    integrator: LangevinBAOAB::with_constraints(
                        1.,
                        1.,
                        seed,
                        f.constraints.clone(),
                    ),
                    nl: NeighborList::new(1.),
                    step: 0,
                    seed,
                    n,
                });
            }
            reps.push(row);
        }
        groups.push((id, reps, false));
    }
    let mut prior = 0;
    for duration in [400, 800, 1600] {
        if groups.iter().all(|(_, _, done)| *done) {
            break;
        }
        if prior != 0 {
            let python =
                std::env::var_os("B07_HOST_PYTHON").ok_or("missing pinned host interpreter")?;
            let launcher = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
                .join("../../bench/b07_validation/campaign.py");
            let status = std::process::Command::new(python)
                .arg(launcher)
                .args(["verify", "--out"])
                .arg(out)
                .status()
                .map_err(err)?;
            if !status.success() {
                return Err("epoch changed before duration extension".into());
            }
        }
        for (id, reps, done) in &mut groups {
            if *done {
                continue;
            }
            for (i, row) in reps.iter_mut().enumerate() {
                for r in row {
                    extend(out, id, i, r, duration, prior, &mut work)?;
                }
            }
            let (rows, precise) = analysis(out, id, duration)?;
            write_new(&out.join(format!("{id}-{duration}-gates.json")), &rows)?;
            if precise || duration == 1600 {
                gates.extend(rows);
                *done = true;
            }
        }
        prior = duration;
    }
    Ok(gates)
}

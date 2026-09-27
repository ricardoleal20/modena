//! SPME Stage A — Smooth Particle-Mesh Ewald *reciprocal energy* only.
//!
//! The order-4 B-spline charge spread + 3D FFT + influence-function sum that
//! reproduces the direct reciprocal Ewald energy on a mesh
//! (`docs/SPME_PLAN.md`). Forces (ik-differentiation of this functional) and
//! the strain-form virial are Stage B/C. Orthorhombic boxes only (the schema
//! is orthorhombic today).
//!
//! # Convergence accounting (the α rule)
//!
//! The only convention freedom is the overall normalization constant `α` —
//! physical `k = 2π·k̂/L` and `κ` are fixed definitions. Stage A sets the
//! mesh-accounted `α = 2π/V · (grid/window corrections)` and the test
//! verifies the SAME normalization reproduces the direct ±k path to the
//! configured tolerance on the deck, a displaced config, AND a different-κ
//! config (different k-weighting). A single stable α proves the functional
//! form; an unstable α means a shape/factor error — the test fails loudly.
//!
//! Influence (Essmann et al. JCP 103, 8577 (1995) form, reduced units):
//!
//! ```text
//! E_rec = α · Σ_{k̂≠0} e^{−k²/4κ²}/k² · |Q̂(k̂)|² / (b3(k1)·b3(k2)·b3(k3))²
//! ```
//!
//! with `b3(k̂) = (4 + 2·cos(2πk̂/N))/6` the DFT of the cubic (order-4)
//! B-spline window at the grid nodes — it cancels the spline sampling of
//! `Q`, leaving the plain-Ewald kernel.

use crate::error::{ModenaError, Result};
use crate::system::SimBox;
use crate::vec3::Vec3;

use num_complex::Complex;
use rustfft::FftPlanner;
/// Cubic (order-4) B-spline weights at fractional `u ∈ [0,1)` for the four
/// node offsets `{-1, 0, 1, 2}`. Sums to 1 exactly.
#[inline]
fn m4(u: f64) -> [f64; 4] {
    let u1 = 1.0 - u;
    [
        u1 * u1 * u1,
        3.0 * u * u * u - 6.0 * u * u + 4.0,
        -3.0 * u * u * u + 3.0 * u * u + 3.0 * u + 1.0,
        u * u * u,
    ]
    .map(|x| x / 6.0)
}

/// Reciprocal SPME energy (Stage A) for an orthorhombic box.
pub fn spme_reciprocal_energy(
    kappa: f64,
    kcut: f64,
    box_size: &SimBox,
    positions: &[f64],
    charges: &[f64],
    mesh: [usize; 3],
) -> Result<f64> {
    if box_size.is_triclinic() {
        return Err(ModenaError::Config(
            "SPME requires an orthorhombic box".into(),
        ));
    }
    let n = charges.len();
    if n == 0 {
        return Ok(0.0);
    }
    let [n0, n1, n2] = mesh;
    let l = [box_size.size.0[0], box_size.size.0[1], box_size.size.0[2]];
    let ntot = n0 * n1 * n2;

    // ---- 1. B-spline charge spread onto the mesh. ----
    let mut q: Vec<f64> = vec![0.0; ntot];
    let mut u = [0.0f64; 3];
    let mut m0 = [0usize; 3];
    let mut w = [[0.0f64; 4]; 3];
    for i in 0..n {
        let xi = Vec3::new(positions[3 * i], positions[3 * i + 1], positions[3 * i + 2]);
        for a in 0..3 {
            // Fractional coordinate in [0,1) (mesh-wrapped).
            let f = (xi.0[a] / l[a]).fract().rem_euclid(1.0);
            let t = f * mesh[a] as f64;
            let m = t.floor() as usize % mesh[a];
            u[a] = t - m as f64;
            m0[a] = m;
            w[a] = m4(u[a]);
        }
        let qc = charges[i];
        for k0 in 0..4usize {
            let a0 = (m0[0] + k0 + n0 - 1) % n0;
            for k1 in 0..4usize {
                let a1 = (m0[1] + k1 + n1 - 1) % n1;
                for k2 in 0..4usize {
                    let a2 = (m0[2] + k2 + n2 - 1) % n2;
                    let weight = w[0][k0] * w[1][k1] * w[2][k2];
                    q[(a0 * n1 + a1) * n2 + a2] += qc * weight;
                }
            }
        }
    }

    // ---- 2. Forward FFT (real spread → complex). ----
    let mut grid: Vec<Complex<f64>> = q.iter().map(|&x| Complex::new(x, 0.0)).collect();
    let mut planner = FftPlanner::<f64>::new();
    let ffts = [
        planner.plan_fft_forward(n0),
        planner.plan_fft_forward(n1),
        planner.plan_fft_forward(n2),
    ];
    let scratch_len = ffts
        .iter()
        .map(|f| f.get_inplace_scratch_len())
        .max()
        .unwrap_or(0);
    let mut scratch: Vec<Complex<f64>> = vec![Complex::new(0.0, 0.0); scratch_len.max(1)];

    // Axis 2 (stride 1): contiguous lines.
    for i0 in 0..n0 {
        for i1 in 0..n1 {
            let start = (i0 * n1 + i1) * n2;
            ffts[2].process_with_scratch(&mut grid[start..start + n2], &mut scratch);
        }
    }
    // Axis 1 (stride n2): strided → transpose copy.
    for i0 in 0..n0 {
        for i2 in 0..n2 {
            let mut line: Vec<Complex<f64>> =
                (0..n1).map(|i1| grid[(i0 * n1 + i1) * n2 + i2]).collect();
            ffts[1].process_with_scratch(&mut line, &mut scratch);
            for (i1, v) in line.into_iter().enumerate() {
                grid[(i0 * n1 + i1) * n2 + i2] = v;
            }
        }
    }
    // Axis 0 (stride n1*n2): strided → transpose copy.
    for i1 in 0..n1 {
        for i2 in 0..n2 {
            let mut line: Vec<Complex<f64>> =
                (0..n0).map(|i0| grid[(i0 * n1 + i1) * n2 + i2]).collect();
            ffts[0].process_with_scratch(&mut line, &mut scratch);
            for (i0, v) in line.into_iter().enumerate() {
                grid[(i0 * n1 + i1) * n2 + i2] = v;
            }
        }
    }

    // ---- 3. Influence-function sum over the full mesh. ----
    // Physical k = 2π·k̂/L with k̂ ∈ [−N/2, N/2) (centered indexing).
    let khat = |nn: usize, m: usize| -> i64 {
        let m = m as i64;
        if m <= (nn as i64) / 2 {
            m
        } else {
            m - nn as i64
        }
    };
    let b3 = |kh: i64, n: usize| -> f64 {
        (4.0 + 2.0 * (2.0 * std::f64::consts::PI * kh as f64 / n as f64).cos()) / 6.0
    };
    let two_pi = 2.0 * std::f64::consts::PI;
    let kappa2 = kappa * kappa;
    let mut e = 0.0f64;
    for n0i in 0..n0 {
        let kx = khat(n0, n0i);
        let phx = two_pi * kx as f64 / l[0];
        for n1i in 0..n1 {
            let ky = khat(n1, n1i);
            let phy = two_pi * ky as f64 / l[1];
            for n2i in 0..n2 {
                let kz = khat(n2, n2i);
                if kx == 0 && ky == 0 && kz == 0 {
                    continue;
                }
                let phz = two_pi * kz as f64 / l[2];
                let k2 = phx * phx + phy * phy + phz * phz;
                if k2 <= 0.0 {
                    continue;
                }
                let qhat = grid[(n0i * n1 + n1i) * n2 + n2i];
                let w = libm::exp(-k2 / (4.0 * kappa2)) / k2;
                let b = b3(kx, n0) * b3(ky, n1) * b3(kz, n2);
                e += w * (qhat.norm_sqr()) / (b * b);
            }
        }
    }
    let _ = kcut;
    let vol = l[0] * l[1] * l[2];
    let alpha = two_pi / vol;
    Ok(alpha * e)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ewald::EwaldForceModel;
    use crate::lattice;

    fn nacl() -> crate::system::System {
        lattice::rocksalt(5.0, [3, 3, 3], "Na", "Cl", 1.0, 1.0, 1.0, -1.0)
    }

    fn run_case(kappa: f64, kcut: f64, mesh: [usize; 3], disp: bool) -> (f64, f64, f64) {
        let sys = nacl();
        let mut sys = sys;
        if disp {
            // Displaced config (tiny random-like displacement → nonzero forces).
            for i in 0..sys.n() {
                let f = i as f64;
                let p = sys.pos(i)
                    + Vec3::new(
                        1e-3 * (0.13 * f).sin(),
                        1e-3 * (0.7 * f).cos(),
                        1e-3 * (1.3 * f).sin(),
                    );
                sys.set_pos(i, p);
            }
            sys.wrap_positions();
        }
        let ewald = EwaldForceModel::new(kappa, 5.0, kcut).expect("ewald");
        let direct = ewald.reciprocal_energy(&sys).expect("direct").0;
        let spme = spme_reciprocal_energy(
            kappa,
            kcut,
            &sys.box_size,
            &sys.positions,
            &sys.charges,
            mesh,
        )
        .expect("spme");
        let ratio = spme / direct;
        (direct, spme, ratio)
    }

    fn mesh_for(box_len: f64, kcut: f64) -> [usize; 3] {
        let n = (box_len * kcut / std::f64::consts::PI).ceil() as usize;
        [n, n, n]
    }

    #[test]
    fn spme_matches_direct_on_deck() {
        let (direct, spme, ratio) = run_case(0.8, 7.0, mesh_for(15.0, 7.0), false);
        assert!(
            (ratio - 1.0).abs() < 5e-3,
            "deck: direct={direct} spme={spme} ratio={ratio}"
        );
    }

    #[test]
    fn spme_matches_direct_displaced() {
        let (direct, spme, ratio) = run_case(0.8, 7.0, mesh_for(15.0, 7.0), true);
        assert!(
            (ratio - 1.0).abs() < 5e-3,
            "displaced: direct={direct} spme={spme} ratio={ratio}"
        );
    }

    #[test]
    fn spme_alpha_stable_across_kappa() {
        // The SAME normalization must reproduce the direct path at a different
        // κ (different k-weighting) — the α-rule.
        let (direct, spme, ratio) = run_case(0.5, 6.0, mesh_for(15.0, 6.0), true);
        assert!(
            (ratio - 1.0).abs() < 5e-3,
            "kappa=0.5: direct={direct} spme={spme} ratio={ratio}"
        );
    }
}

//! SHAKE (position) and RATTLE (velocity) constraints for fixed bond lengths.
//!
//! A `Constraints` set holds rigid distance constraints `|r_ij| = d` (typically
//! the X–H bonds, or all bonds in a rigid-molecule model). The integrator
//! applies them after the unconstrained position drift ([`Constraints::shake`])
//! and after the final velocity half-kick ([`Constraints::rattle`]).
//!
//! SHAKE solves the position constraints iteratively (Gauss–Seidel over the
//! constraint list, `max_iter` sweeps) using the linearized Lagrange-multiplier
//! update
//!
//! ```text
//!   g = (|r_ij'|² − d²) / (2 (1/m_i + 1/m_j) (r_ij' · r_ij(t)))
//!   Δx_i = +(g/m_i) r_ij(t),   Δx_j = −(g/m_j) r_ij(t)
//! ```
//!
//! where `r_ij'` is the predicted (post-drift) bond vector and `r_ij(t)` the
//! reference vector at the start of the step. RATTLE then projects the
//! velocities onto the constraint tangent (`r_ij · v_ij = 0`) with
//!
//! ```text
//!   μ = (r_ij · v_ij) / ((1/m_i + 1/m_j) |r_ij|²)
//!   Δv_i = +(μ/m_i) r_ij,   Δv_j = −(μ/m_j) r_ij
//! ```
//!
//! For a single, non-overlapping set of constraints this converges in a handful
//! of iterations; the tolerance and iteration cap are configurable. The
//! combination SHAKE+RATTLE on velocity-Verlet is the symplectic RATTLE scheme
//! (Andersen, 1983) — it conserves a shadow Hamiltonian on the constraint
//! manifold, verified by the constrained-NVE test in `tests/bonded.rs`.

use rayon::prelude::*;

use crate::system::{SimBox, System};
use crate::vec3::Vec3;

/// One rigid distance constraint: atoms `i`, `j` held at fixed length `d`.
#[derive(Debug, Clone, Copy)]
pub struct Constraint {
    pub i: usize,
    pub j: usize,
    pub d: f64,
}

/// One rigid water molecule: the atoms at `first`, `first + 1`, `first + 2`
/// (O, H, H — water builders store each molecule as a contiguous triple, which
/// is what lets molecules solve in parallel: one 3-atom chunk per ray) held at
/// O–H length `d_oh` and H–H length `d_hh` (SPC/E: 1.0 Å / 1.633 Å). Solved
/// per molecule with an exact Newton projection of the three coupled
/// constraints — the SETTLE role (Miyamoto & Kollman 1992 solve it in closed
/// form; this is the exact Newton equivalent, M-SHAKE-class) — which, unlike
/// sequential SHAKE over a bond triangle, cannot fail to converge at large dt.
#[derive(Debug, Clone, Copy)]
pub struct RigidWater {
    pub first: usize,
    pub d_oh: f64,
    pub d_hh: f64,
}

/// Newton cap for the per-molecule rigid-water projection. From a violated
/// state one iteration lands at O(violation²); three reach machine precision
/// for any physically attainable distortion — the cap is a runaway backstop.
const WATER_NEWTON_ITER: usize = 20;

/// A set of distance constraints applied by the integrator: loose SHAKE/RATTLE
/// `bonds` plus rigid `waters` solved exactly per molecule.
#[derive(Debug, Clone, Default)]
pub struct Constraints {
    pub bonds: Vec<Constraint>,
    pub waters: Vec<RigidWater>,
}

impl Constraints {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn is_empty(&self) -> bool {
        self.bonds.is_empty() && self.waters.is_empty()
    }

    /// Exact per-molecule position projection for the rigid waters: Newton on
    /// the three multipliers `(λ_OH1, λ_OH2, λ_HH)` of the coupled constraint
    /// system, with the constraint directions re-taken each iteration. Each
    /// molecule touches only its own three atoms, so molecules solve in
    /// parallel. COM is preserved exactly (every constraint impulse is
    /// mass-weighted and internal).
    pub fn settle(&self, system: &mut System, tol: f64) {
        if self.waters.is_empty() {
            return;
        }
        let box_size = system.box_size;
        let masses: &[f64] = &system.masses;
        let d2 = |v: Vec3| v.dot(v);
        // Water molecules are contiguous 3-atom triples, so 9-float chunk w
        // of the position buffer IS water w — disjoint &mut chunks, one ray
        // each.
        system
            .positions
            .par_chunks_exact_mut(9)
            .enumerate()
            .for_each(|(wi, chunk)| {
                let w = self.waters[wi];
                debug_assert_eq!(w.first, 3 * wi);
                settle_one(&box_size, chunk, masses, w, tol, &d2);
            });
    }

    /// Exact per-molecule velocity projection for the rigid waters: the same
    /// 3×3 system as the position solve but linear (one solve, no iteration)
    /// — projects `r_ij · v_ij = 0` for all three constraints at once.
    pub fn settle_velocities(&self, system: &mut System, tol: f64) {
        if self.waters.is_empty() {
            return;
        }
        let box_size = system.box_size;
        let positions: &[f64] = &system.positions;
        let masses: &[f64] = &system.masses;
        system
            .velocities
            .par_chunks_exact_mut(9)
            .enumerate()
            .for_each(|(wi, chunk)| {
                let w = self.waters[wi];
                debug_assert_eq!(w.first, 3 * wi);
                settle_vel_one(&box_size, positions, masses, chunk, w, tol);
            });
    }

    /// SHAKE: correct the predicted positions in `system` so every constraint
    /// `|r_ij| = d` holds. `x_ref` is the flattened `3N` position snapshot at
    /// the start of the step (the reference bond vectors `r_ij(t)`).
    ///
    /// Iterates until the worst constraint residual is below `tol` (in length²)
    /// or `max_iter` sweeps are done.
    pub fn shake(&self, system: &mut System, x_ref: &[f64], tol: f64, max_iter: usize) {
        // Rigid waters solve exactly per molecule (they carry no `bonds`
        // entries); the loose bond set, if any, follows the Gauss–Seidel pass.
        self.settle(system, tol);
        if self.bonds.is_empty() {
            return;
        }
        let box_size = system.box_size;
        for _ in 0..max_iter {
            let mut max_err = 0.0_f64;
            for &c in &self.bonds {
                let ri = system.pos(c.i);
                let rj = system.pos(c.j);
                let rij = box_size.min_image(rj - ri); // predicted, j − i
                let r2 = rij.norm2();
                let d2 = c.d * c.d;
                let err = (r2 - d2).abs();
                if err > max_err {
                    max_err = err;
                }
                if err < tol {
                    continue;
                }
                // Reference bond vector r_ij(t) = x_j(t) − x_i(t).
                let ri0 = Vec3::new(x_ref[3 * c.i], x_ref[3 * c.i + 1], x_ref[3 * c.i + 2]);
                let rj0 = Vec3::new(x_ref[3 * c.j], x_ref[3 * c.j + 1], x_ref[3 * c.j + 2]);
                let rij0 = box_size.min_image(rj0 - ri0);
                let mi = system.masses[c.i];
                let mj = system.masses[c.j];
                let inv = 1.0 / mi + 1.0 / mj;
                let denom = 2.0 * inv * rij.dot(rij0);
                if denom.abs() < 1e-20 {
                    continue;
                }
                let g = (r2 - d2) / denom;
                // Δx_i = +(g/m_i) r_ij0 ; Δx_j = −(g/m_j) r_ij0  (move the atoms
                // along the reference bond so |r_ij| → d; mass-weighted so the
                // lighter atom carries more of the correction).
                system.add_pos(c.i, rij0 * (g / mi));
                system.add_pos(c.j, rij0 * (-g / mj));
            }
            if max_err < tol {
                break;
            }
        }
    }

    /// RATTLE: project the velocities in `system` onto the constraint tangent so
    /// that `r_ij · v_ij = 0` for every constraint. Gauss–Seidel iteration.
    pub fn rattle(&self, system: &mut System, tol: f64, max_iter: usize) {
        self.settle_velocities(system, tol);
        if self.bonds.is_empty() {
            return;
        }
        let box_size = system.box_size;
        for _ in 0..max_iter {
            let mut max_err = 0.0_f64;
            for &c in &self.bonds {
                let rij = box_size.min_image(system.pos(c.j) - system.pos(c.i));
                let vij = system.vel(c.j) - system.vel(c.i);
                let rv = rij.dot(vij);
                if rv.abs() > max_err {
                    max_err = rv.abs();
                }
                if rv.abs() < tol {
                    continue;
                }
                let mi = system.masses[c.i];
                let mj = system.masses[c.j];
                let inv = 1.0 / mi + 1.0 / mj;
                let denom = inv * rij.norm2();
                if denom.abs() < 1e-20 {
                    continue;
                }
                let mu = rv / denom;
                // Δv_i = +(μ/m_i) r_ij ; Δv_j = −(μ/m_j) r_ij.
                system.add_vel(c.i, rij * (mu / mi));
                system.add_vel(c.j, rij * (-mu / mj));
            }
            if max_err < tol {
                break;
            }
        }
    }

    /// Worst constraint length-squared residual (diagnostic, used in tests).
    pub fn max_residual(&self, system: &System) -> f64 {
        let box_size = system.box_size;
        self.bonds
            .iter()
            .map(|c| {
                let rij = box_size.min_image(system.pos(c.j) - system.pos(c.i));
                (rij.norm2() - c.d * c.d).abs()
            })
            .fold(0.0_f64, f64::max)
    }
}

/// Solve one water's three position constraints with Newton on the
/// multipliers. `g_k = |v_k|² − d_k²` with `v_k` the current (min-imaged)
/// constraint vectors; the Jacobian is the projected Hessian
/// `J_kl = 2 v_k · (∂x_{j_k} − ∂x_{i_k})/∂λ_l` with
/// `∂x_i/∂λ_l = ∓v_l/m_i` — symmetric, as it must be.
fn settle_one(
    box_size: &SimBox,
    chunk: &mut [f64],
    masses: &[f64],
    w: RigidWater,
    tol: f64,
    d2: &dyn Fn(Vec3) -> f64,
) {
    let (mo, m1, m2) = (masses[w.first], masses[w.first + 1], masses[w.first + 2]);
    let doh2 = w.d_oh * w.d_oh;
    let dhh2 = w.d_hh * w.d_hh;
    let p = |chunk: &[f64], k: usize| Vec3::new(chunk[3 * k], chunk[3 * k + 1], chunk[3 * k + 2]);
    for _ in 0..WATER_NEWTON_ITER {
        let (o, h1, h2) = (p(chunk, 0), p(chunk, 1), p(chunk, 2));
        let v1 = box_size.min_image(h1 - o); // O–H1
        let v2 = box_size.min_image(h2 - o); // O–H2
        let v3 = box_size.min_image(h2 - h1); // H1–H2
        let g = [d2(v1) - doh2, d2(v2) - doh2, d2(v3) - dhh2];
        if g[0].abs() < tol && g[1].abs() < tol && g[2].abs() < tol {
            return;
        }
        // J (symmetric): rows = OH1, OH2, HH.
        let v1v2 = v1.dot(v2);
        let v1v3 = v1.dot(v3);
        let v2v3 = v2.dot(v3);
        let j = [
            [
                2.0 * d2(v1) * (1.0 / mo + 1.0 / m1),
                2.0 * v1v2 / mo,
                -2.0 * v1v3 / m1,
            ],
            [
                2.0 * v1v2 / mo,
                2.0 * d2(v2) * (1.0 / mo + 1.0 / m2),
                2.0 * v2v3 / m2,
            ],
            [
                -2.0 * v1v3 / m1,
                2.0 * v2v3 / m2,
                2.0 * d2(v3) * (1.0 / m1 + 1.0 / m2),
            ],
        ];
        let Some(lam) = solve3(&j, &[-g[0], -g[1], -g[2]]) else {
            return; // degenerate molecule (coincident atoms) — leave as-is
        };
        // x_i ← x_i ∓ (λ_k/m_i) v_k (mass-weighted, internal → COM exact).
        add3(chunk, 0, (v1 * lam[0] + v2 * lam[1]) * (-1.0 / mo));
        add3(chunk, 1, (v1 * lam[0] - v3 * lam[2]) * (1.0 / m1));
        add3(chunk, 2, (v2 * lam[1] + v3 * lam[2]) * (1.0 / m2));
    }
}

/// Project one water's velocities so `r_k · v_jk = 0` for all three
/// constraints — the same 3×3 system as the position Jacobian (with the
/// position vectors in place of the constraint vectors), linear in the
/// multipliers: one solve.
fn settle_vel_one(
    box_size: &SimBox,
    positions: &[f64],
    masses: &[f64],
    velocities: &mut [f64],
    w: RigidWater,
    tol: f64,
) {
    let (mo, m1, m2) = (masses[w.first], masses[w.first + 1], masses[w.first + 2]);
    let f = positions[w.first * 3..w.first * 3 + 9].to_vec();
    let pos = |k: usize| Vec3::new(f[3 * k], f[3 * k + 1], f[3 * k + 2]);
    let vel = |k: usize| {
        Vec3::new(
            velocities[3 * k],
            velocities[3 * k + 1],
            velocities[3 * k + 2],
        )
    };
    let r1 = box_size.min_image(pos(1) - pos(0));
    let r2 = box_size.min_image(pos(2) - pos(0));
    let r3 = box_size.min_image(pos(2) - pos(1));
    let c = [
        r1.dot(vel(1) - vel(0)),
        r2.dot(vel(2) - vel(0)),
        r3.dot(vel(2) - vel(1)),
    ];
    if c[0].abs() < tol && c[1].abs() < tol && c[2].abs() < tol {
        return;
    }
    let m = [
        [
            r1.dot(r1) * (1.0 / mo + 1.0 / m1),
            r1.dot(r2) / mo,
            -r1.dot(r3) / m1,
        ],
        [
            r1.dot(r2) / mo,
            r2.dot(r2) * (1.0 / mo + 1.0 / m2),
            r2.dot(r3) / m2,
        ],
        [
            -r1.dot(r3) / m1,
            r2.dot(r3) / m2,
            r3.dot(r3) * (1.0 / m1 + 1.0 / m2),
        ],
    ];
    let Some(mu) = solve3(&m, &[-c[0], -c[1], -c[2]]) else {
        return;
    };
    add3(velocities, 0, (r1 * mu[0] + r2 * mu[1]) * (-1.0 / mo));
    add3(velocities, 1, (r1 * mu[0] - r3 * mu[2]) * (1.0 / m1));
    add3(velocities, 2, (r2 * mu[1] + r3 * mu[2]) * (1.0 / m2));
}

/// Add `d` to atom `k` of a contiguous 3-atom chunk buffer.
#[inline]
fn add3(buf: &mut [f64], k: usize, d: Vec3) {
    buf[3 * k] += d.0[0];
    buf[3 * k + 1] += d.0[1];
    buf[3 * k + 2] += d.0[2];
}

/// Symmetric 3×3 solve via Cramer's rule. `None` on a (near-)singular system.
fn solve3(a: &[[f64; 3]; 3], b: &[f64; 3]) -> Option<[f64; 3]> {
    let det = a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1])
        - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
        + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0]);
    if det.abs() < 1e-300 {
        return None;
    }
    let x0 = b[0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1])
        - a[0][1] * (b[1] * a[2][2] - a[1][2] * b[2])
        + a[0][2] * (b[1] * a[2][1] - a[1][1] * b[2]);
    let x1 = a[0][0] * (b[1] * a[2][2] - a[1][2] * b[2])
        - b[0] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
        + a[0][2] * (a[1][0] * b[2] - b[1] * a[2][0]);
    let x2 = a[0][0] * (a[1][1] * b[2] - b[1] * a[2][1])
        - a[0][1] * (a[1][0] * b[2] - b[1] * a[2][0])
        + b[0] * (a[1][0] * a[2][1] - a[1][1] * a[2][0]);
    Some([x0 / det, x1 / det, x2 / det])
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::system::SimBox;

    /// SHAKE/RATTLE on a single bond: two atoms, stretched predictor, fixed to
    /// the target length, and the velocity projected onto the tangent.
    #[test]
    fn shake_single_bond_satisfies_constraint() {
        let mut sys = System::new(SimBox::cubic(100.0));
        sys.resize(2);
        sys.set_pos(0, Vec3::new(0.0, 0.0, 0.0));
        sys.set_pos(1, Vec3::new(1.5, 0.0, 0.0)); // stretched (target 1.0)
                                                  // Reference (start-of-step) bond at equilibrium along x.
        let x_ref = vec![0.0, 0.0, 0.0, 1.0, 0.0, 0.0];
        let cons = Constraints {
            bonds: vec![Constraint { i: 0, j: 1, d: 1.0 }],
            waters: vec![],
        };
        cons.shake(&mut sys, &x_ref, 1e-12, 50);
        let r = sys.box_size.min_image(sys.pos(1) - sys.pos(0));
        assert!(
            (r.norm() - 1.0).abs() < 1e-9,
            "bond length after SHAKE: {}",
            r.norm()
        );
        // COM preserved (equal masses): the predicted midpoint (0.75) is unchanged.
        let com = (sys.pos(0) + sys.pos(1)) * 0.5;
        assert!(
            (com - Vec3::new(0.75, 0.0, 0.0)).norm() < 1e-9,
            "COM drifted: {com:?}"
        );
    }

    #[test]
    fn rattle_removes_radial_velocity() {
        let mut sys = System::new(SimBox::cubic(100.0));
        sys.resize(2);
        sys.set_pos(0, Vec3::new(0.0, 0.0, 0.0));
        sys.set_pos(1, Vec3::new(1.0, 0.0, 0.0));
        // Pure radial velocity (stretches the bond) — must be removed.
        sys.set_vel(0, Vec3::new(-1.0, 0.0, 0.0));
        sys.set_vel(1, Vec3::new(2.0, 0.5, 0.0));
        let cons = Constraints {
            bonds: vec![Constraint { i: 0, j: 1, d: 1.0 }],
            waters: vec![],
        };
        cons.rattle(&mut sys, 1e-12, 50);
        let rij = sys.box_size.min_image(sys.pos(1) - sys.pos(0));
        let vij = sys.vel(1) - sys.vel(0);
        assert!(
            rij.dot(vij).abs() < 1e-9,
            "radial velocity after RATTLE: {}",
            rij.dot(vij)
        );
    }

    /// SPC/E-geometry water, badly distorted, restored to machine precision
    /// by the exact per-molecule Newton projection — with COM preserved
    /// (the impulses are internal and mass-weighted).
    #[test]
    fn settle_restores_distorted_water() {
        let d_oh = 1.0_f64;
        let d_hh = 2.0 * d_oh * (109.47_f64.to_radians() / 2.0).sin();
        let mut sys = System::new(SimBox::cubic(20.0));
        sys.resize(3);
        let half = 109.47_f64.to_radians() / 2.0;
        sys.set_pos(0, Vec3::new(0.0, 0.0, 0.0));
        sys.set_pos(1, Vec3::new(d_oh * half.cos(), d_oh * half.sin(), 0.0));
        sys.set_pos(2, Vec3::new(d_oh * half.cos(), -d_oh * half.sin(), 0.0));
        sys.masses = vec![15.9994, 1.008, 1.008];
        // Distort far beyond any physical dt (15% bond stretch).
        sys.set_pos(1, sys.pos(1) + Vec3::new(0.15, -0.1, 0.07));
        sys.set_pos(2, sys.pos(2) + Vec3::new(-0.12, 0.09, 0.11));
        let com_before = mass_weighted_com(&sys);
        let cons = Constraints {
            bonds: vec![],
            waters: vec![RigidWater {
                first: 0,
                d_oh,
                d_hh,
            }],
        };
        cons.settle(&mut sys, 1e-12);
        let oh1 = (sys.pos(1) - sys.pos(0)).norm();
        let oh2 = (sys.pos(2) - sys.pos(0)).norm();
        let hh = (sys.pos(2) - sys.pos(1)).norm();
        assert!((oh1 - d_oh).abs() < 1e-10, "OH1 {oh1}");
        assert!((oh2 - d_oh).abs() < 1e-10, "OH2 {oh2}");
        assert!((hh - d_hh).abs() < 1e-10, "HH {hh}");
        let com_after = mass_weighted_com(&sys);
        assert!(
            (com_after - com_before).norm() < 1e-12,
            "COM moved: {com_before:?} -> {com_after:?}"
        );
    }

    /// The velocity projection removes the radial component along all three
    /// constraint directions at once (one exact 3×3 solve).
    #[test]
    fn settle_projects_water_velocities() {
        let d_oh = 1.0_f64;
        let d_hh = 2.0 * d_oh * (109.47_f64.to_radians() / 2.0).sin();
        let mut sys = System::new(SimBox::cubic(20.0));
        sys.resize(3);
        let half = 109.47_f64.to_radians() / 2.0;
        sys.set_pos(0, Vec3::new(0.0, 0.0, 0.0));
        sys.set_pos(1, Vec3::new(d_oh * half.cos(), d_oh * half.sin(), 0.0));
        sys.set_pos(2, Vec3::new(d_oh * half.cos(), -d_oh * half.sin(), 0.0));
        sys.masses = vec![15.9994, 1.008, 1.008];
        // Radial (stretching) + tangential velocities on every bond.
        sys.set_vel(0, Vec3::new(0.3, -0.2, 0.1));
        sys.set_vel(1, Vec3::new(-0.5, 0.4, 0.2));
        sys.set_vel(2, Vec3::new(0.2, 0.1, -0.3));
        let cons = Constraints {
            bonds: vec![],
            waters: vec![RigidWater {
                first: 0,
                d_oh,
                d_hh,
            }],
        };
        cons.settle_velocities(&mut sys, 1e-12);
        for (i, j) in [(0, 1), (0, 2), (1, 2)] {
            let r = sys.pos(j) - sys.pos(i);
            let v = sys.vel(j) - sys.vel(i);
            assert!(
                r.dot(v).abs() < 1e-10,
                "radial velocity on bond ({i},{j}): {}",
                r.dot(v)
            );
        }
    }

    fn mass_weighted_com(sys: &System) -> Vec3 {
        let mut m = 0.0;
        let mut c = Vec3::zero();
        for i in 0..sys.n() {
            c += sys.pos(i) * sys.masses[i];
            m += sys.masses[i];
        }
        c * (1.0 / m)
    }
}

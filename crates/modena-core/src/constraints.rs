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

use crate::system::System;
use crate::vec3::Vec3;

/// One rigid distance constraint: atoms `i`, `j` held at fixed length `d`.
#[derive(Debug, Clone, Copy)]
pub struct Constraint {
    pub i: usize,
    pub j: usize,
    pub d: f64,
}

/// A set of distance constraints applied by the integrator.
#[derive(Debug, Clone, Default)]
pub struct Constraints {
    pub bonds: Vec<Constraint>,
}

impl Constraints {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn is_empty(&self) -> bool {
        self.bonds.is_empty()
    }

    /// SHAKE: correct the predicted positions in `system` so every constraint
    /// `|r_ij| = d` holds. `x_ref` is the flattened `3N` position snapshot at
    /// the start of the step (the reference bond vectors `r_ij(t)`).
    ///
    /// Iterates until the worst constraint residual is below `tol` (in length²)
    /// or `max_iter` sweeps are done.
    pub fn shake(&self, system: &mut System, x_ref: &[f64], tol: f64, max_iter: usize) {
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
}

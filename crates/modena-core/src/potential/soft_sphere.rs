//! Soft-sphere pair potential: `U(r) = ε·(σ/r)ⁿ`.
//!
//! A purely repulsive power-law potential (no attractive tail), useful for
//! coarse-grained and hard-sphere-like models. Shifted-force at the cutoff.
//! Multi-species via per-type-pair `(σ, ε, n)` tables.

use crate::error::{ModenaError, Result};
use crate::potential::check_cutoff;
use crate::potential::mixing::{lorentz_berthelot, MixingRule, PairParamTable};

use super::shifted_force;
use super::PairPotential;

/// Per-type soft-sphere parameters.
#[derive(Debug, Clone, Copy)]
pub struct SoftSphereTypeParams {
    pub sigma: f64,
    pub epsilon: f64,
    /// Repulsive exponent `n` (typically 12 or 18).
    pub n: i32,
}

/// Explicit per-pair soft-sphere override.
#[derive(Debug, Clone, Copy)]
pub struct SoftSphereOverride {
    pub types: [u32; 2],
    pub sigma: f64,
    pub epsilon: f64,
    pub n: i32,
}

/// Soft-sphere `ε·(σ/r)ⁿ`, shifted-force, multi-species.
#[derive(Debug, Clone)]
pub struct SoftSphere {
    sigma: PairParamTable,
    epsilon: PairParamTable,
    n: PairParamTable,
    u_rc: PairParamTable,
    f_rc: PairParamTable,
    cutoff: f64,
    cutoff2: f64,
    n_types: usize,
}

impl SoftSphere {
    /// Single-species constructor.
    pub fn new(sigma: f64, epsilon: f64, n: i32, cutoff: f64) -> Result<Self> {
        #![allow(clippy::neg_cmp_op_on_partial_ord)]
        if !(sigma > 0.0) {
            return Err(ModenaError::Config("soft-sphere sigma must be > 0".into()));
        }
        if !(epsilon >= 0.0) {
            return Err(ModenaError::Config(
                "soft-sphere epsilon must be >= 0".into(),
            ));
        }
        if n <= 0 {
            return Err(ModenaError::Config("soft-sphere n must be > 0".into()));
        }
        check_cutoff(cutoff, "soft-sphere")?;
        Self::from_types_inner(&[SoftSphereTypeParams { sigma, epsilon, n }], &[], cutoff)
    }

    /// Multi-species constructor. Cross-pair σ uses Lorentz–Berthelot (arithmetic
    /// mean), ε uses the geometric mean; `n` is taken from the diagonal. Overrides
    /// take precedence.
    pub fn from_types(
        types: &[SoftSphereTypeParams],
        _mixing: MixingRule,
        overrides: &[SoftSphereOverride],
        cutoff: f64,
    ) -> Result<Self> {
        for (i, t) in types.iter().enumerate() {
            if !(t.sigma > 0.0) {
                return Err(ModenaError::Config(format!(
                    "soft-sphere sigma for type {i} must be > 0"
                )));
            }
            if !(t.epsilon >= 0.0) {
                return Err(ModenaError::Config(format!(
                    "soft-sphere epsilon for type {i} must be >= 0"
                )));
            }
            if t.n <= 0 {
                return Err(ModenaError::Config(format!(
                    "soft-sphere n for type {i} must be > 0"
                )));
            }
        }
        check_cutoff(cutoff, "soft-sphere")?;
        Self::from_types_inner(types, overrides, cutoff)
    }

    fn from_types_inner(
        types: &[SoftSphereTypeParams],
        overrides: &[SoftSphereOverride],
        cutoff: f64,
    ) -> Result<Self> {
        let n = types.len();
        let cutoff2 = cutoff * cutoff;
        let mut sigma = PairParamTable::new(n);
        let mut epsilon = PairParamTable::new(n);
        let mut n_table = PairParamTable::new(n);
        for i in 0..n {
            for j in 0..n {
                let s = lorentz_berthelot::sigma(types[i].sigma, types[j].sigma);
                let e = lorentz_berthelot::epsilon(types[i].epsilon, types[j].epsilon);
                sigma.set_symmetric(i, j, s);
                epsilon.set_symmetric(i, j, e);
                // Exponent: take the larger (more repulsive) of the two.
                n_table.set_symmetric(i, j, types[i].n.max(types[j].n) as f64);
            }
        }
        for ov in overrides {
            let [x, y] = ov.types;
            let (x, y) = (x as usize, y as usize);
            if x >= n || y >= n {
                return Err(ModenaError::Config(format!(
                    "soft-sphere override type index {x}/{y} out of range (n_types={n})"
                )));
            }
            sigma.set_symmetric(x, y, ov.sigma);
            epsilon.set_symmetric(x, y, ov.epsilon);
            n_table.set_symmetric(x, y, ov.n as f64);
        }
        let mut u_rc = PairParamTable::new(n);
        let mut f_rc = PairParamTable::new(n);
        for i in 0..n {
            for j in 0..n {
                let (u, f) = raw_energy_force(
                    cutoff,
                    sigma.get(i as u32, j as u32),
                    epsilon.get(i as u32, j as u32),
                    n_table.get(i as u32, j as u32) as i32,
                );
                u_rc.set_symmetric(i, j, u);
                f_rc.set_symmetric(i, j, f);
            }
        }
        Ok(Self {
            sigma,
            epsilon,
            n: n_table,
            u_rc,
            f_rc,
            cutoff,
            cutoff2,
            n_types: n,
        })
    }

    pub fn n_types(&self) -> usize {
        self.n_types
    }
}

/// Raw (unshifted) soft-sphere energy and radial force magnitude at `r`.
#[inline]
fn raw_energy_force(r: f64, sigma: f64, epsilon: f64, n: i32) -> (f64, f64) {
    let s_over_r = sigma / r;
    let u = epsilon * s_over_r.powi(n);
    // F_radial = −dU/dr = n·ε·σⁿ / rⁿ⁺¹  (positive = repulsive).
    let f = (n as f64) * epsilon * sigma.powi(n) / r.powi(n + 1);
    (u, f)
}

impl PairPotential for SoftSphere {
    #[inline]
    fn cutoff(&self) -> f64 {
        self.cutoff
    }

    #[inline]
    fn cutoff2(&self) -> f64 {
        self.cutoff2
    }

    #[inline]
    fn energy_force(&self, r2: f64, ti: u32, tj: u32) -> (f64, f64) {
        let r = r2.sqrt();
        let (u_raw, f_raw) = raw_energy_force(
            r,
            self.sigma.get(ti, tj),
            self.epsilon.get(ti, tj),
            self.n.get(ti, tj) as i32,
        );
        let (u, f_sf) = shifted_force::shift(
            u_raw,
            f_raw,
            r,
            self.cutoff,
            self.u_rc.get(ti, tj),
            self.f_rc.get(ti, tj),
        );
        (u, f_sf / r)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn ss() -> SoftSphere {
        SoftSphere::new(1.0, 1.0, 12, 2.5).unwrap()
    }

    #[test]
    fn energy_and_force_zero_at_cutoff() {
        let p = ss();
        let (u, f_over_r) = p.energy_force(p.cutoff2, 0, 0);
        assert!(u.abs() < 1e-9, "shifted U at cutoff = {u}");
        assert!((f_over_r * p.cutoff()).abs() < 1e-9, "shifted F at cutoff");
    }

    #[test]
    fn purely_repulsive() {
        let p = ss();
        for &r in &[0.9, 1.0, 1.5, 2.0] {
            let (_, f_over_r) = p.energy_force(r * r, 0, 0);
            assert!(
                (f_over_r * r) > 0.0,
                "soft-sphere must be repulsive at r={r}"
            );
        }
    }
}

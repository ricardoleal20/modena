//! Buckingham pair potential: `U(r) = A·exp(−r/ρ) − B/r⁶`.
//!
//! A short-range repulsive exponential plus a dispersive `−r⁻⁶` attraction,
//! common in ionic/biomolecular force fields. Shifted-force at the cutoff so
//! both `U` and `F` are continuous. Multi-species via per-type-pair `(A, B, ρ)`
//! tables.

use crate::error::{ModenaError, Result};
use crate::potential::check_cutoff;
use crate::potential::mixing::{MixingRule, PairParamTable};

use super::shifted_force;
use super::PairPotential;

/// Per-type Buckingham parameters.
#[derive(Debug, Clone, Copy)]
pub struct BuckinghamTypeParams {
    /// Repulsive prefactor `A` (energy units).
    pub a: f64,
    /// Dispersive prefactor `B` (energy·length⁶ units).
    pub b: f64,
    /// Repulsive decay length `ρ` (length units).
    pub rho: f64,
}

/// Explicit per-pair Buckingham override.
#[derive(Debug, Clone, Copy)]
pub struct BuckinghamOverride {
    pub types: [u32; 2],
    pub a: f64,
    pub b: f64,
    pub rho: f64,
}

/// Buckingham `A·exp(−r/ρ) − B/r⁶`, shifted-force, multi-species.
#[derive(Debug, Clone)]
pub struct Buckingham {
    a: PairParamTable,
    b: PairParamTable,
    rho: PairParamTable,
    u_rc: PairParamTable,
    f_rc: PairParamTable,
    cutoff: f64,
    cutoff2: f64,
    n_types: usize,
}

impl Buckingham {
    /// Single-species constructor.
    pub fn new(a: f64, b: f64, rho: f64, cutoff: f64) -> Result<Self> {
        #![allow(clippy::neg_cmp_op_on_partial_ord)]
        if !(rho > 0.0) {
            return Err(ModenaError::Config("Buckingham rho must be > 0".into()));
        }
        if !(a >= 0.0) {
            return Err(ModenaError::Config("Buckingham A must be >= 0".into()));
        }
        if !(b >= 0.0) {
            return Err(ModenaError::Config("Buckingham B must be >= 0".into()));
        }
        check_cutoff(cutoff, "Buckingham")?;
        Self::from_types_inner(&[BuckinghamTypeParams { a, b, rho }], &[], cutoff)
    }

    /// Multi-species constructor. Cross-pair parameters default to the
    /// single-type values of type 0 (Buckingham has no universal combining
    /// rule); supply explicit overrides for real cross interactions.
    pub fn from_types(
        types: &[BuckinghamTypeParams],
        _mixing: MixingRule,
        overrides: &[BuckinghamOverride],
        cutoff: f64,
    ) -> Result<Self> {
        for (i, t) in types.iter().enumerate() {
            if !(t.rho > 0.0) {
                return Err(ModenaError::Config(format!(
                    "Buckingham rho for type {i} must be > 0"
                )));
            }
            if !(t.a >= 0.0) {
                return Err(ModenaError::Config(format!(
                    "Buckingham A for type {i} must be >= 0"
                )));
            }
            if !(t.b >= 0.0) {
                return Err(ModenaError::Config(format!(
                    "Buckingham B for type {i} must be >= 0"
                )));
            }
        }
        check_cutoff(cutoff, "Buckingham")?;
        Self::from_types_inner(types, overrides, cutoff)
    }

    fn from_types_inner(
        types: &[BuckinghamTypeParams],
        overrides: &[BuckinghamOverride],
        cutoff: f64,
    ) -> Result<Self> {
        let n = types.len();
        let cutoff2 = cutoff * cutoff;
        let mut a = PairParamTable::new(n);
        let mut b = PairParamTable::new(n);
        let mut rho = PairParamTable::new(n);
        for i in 0..n {
            for j in 0..n {
                // No universal combining rule: use the diagonal (i==i) value for
                // cross pairs as a placeholder until an override is supplied.
                let p = &types[i.min(j)];
                a.set_symmetric(i, j, p.a);
                b.set_symmetric(i, j, p.b);
                rho.set_symmetric(i, j, p.rho);
            }
        }
        for ov in overrides {
            let [x, y] = ov.types;
            let (x, y) = (x as usize, y as usize);
            if x >= n || y >= n {
                return Err(ModenaError::Config(format!(
                    "Buckingham override type index {x}/{y} out of range (n_types={n})"
                )));
            }
            a.set_symmetric(x, y, ov.a);
            b.set_symmetric(x, y, ov.b);
            rho.set_symmetric(x, y, ov.rho);
        }
        let mut u_rc = PairParamTable::new(n);
        let mut f_rc = PairParamTable::new(n);
        for i in 0..n {
            for j in 0..n {
                let (u, f) = raw_energy_force(
                    cutoff,
                    a.get(i as u32, j as u32),
                    b.get(i as u32, j as u32),
                    rho.get(i as u32, j as u32),
                );
                u_rc.set_symmetric(i, j, u);
                f_rc.set_symmetric(i, j, f);
            }
        }
        Ok(Self {
            a,
            b,
            rho,
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

/// Raw (unshifted) Buckingham energy and radial force magnitude at distance `r`.
#[inline]
fn raw_energy_force(r: f64, a: f64, b: f64, rho: f64) -> (f64, f64) {
    let r6 = r.powi(6);
    let u = a * (-r / rho).exp() - b / r6;
    // F_radial = −dU/dr = (A/ρ)·exp(−r/ρ) − 6B/r⁷  (positive = repulsive).
    let f = (a / rho) * (-r / rho).exp() - 6.0 * b / (r6 * r);
    (u, f)
}

impl PairPotential for Buckingham {
    #[inline]
    fn cutoff(&self) -> f64 {
        self.cutoff
    }
    #[inline]
    fn n_types(&self) -> usize {
        self.n_types
    }

    #[inline]
    fn cutoff2(&self) -> f64 {
        self.cutoff2
    }

    #[inline]
    fn energy_force(&self, r2: f64, ti: u32, tj: u32) -> (f64, f64) {
        let r = r2.sqrt();
        let a = self.a.get(ti, tj);
        let b = self.b.get(ti, tj);
        let rho = self.rho.get(ti, tj);
        let (u_raw, f_raw) = raw_energy_force(r, a, b, rho);
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

    fn buck() -> Buckingham {
        // A representative Buckingham form.
        Buckingham::new(1000.0, 1.0, 0.3, 2.5).unwrap()
    }

    #[test]
    fn energy_and_force_zero_at_cutoff() {
        let p = buck();
        let (u, f_over_r) = p.energy_force(p.cutoff2, 0, 0);
        assert!(u.abs() < 1e-9, "shifted U at cutoff = {u}");
        assert!((f_over_r * p.cutoff()).abs() < 1e-9, "shifted F at cutoff");
    }

    #[test]
    fn no_nan_inside_cutoff() {
        let p = buck();
        let r = 0.5;
        let (u, f) = p.energy_force(r * r, 0, 0);
        assert!(u.is_finite());
        assert!(f.is_finite());
    }
}

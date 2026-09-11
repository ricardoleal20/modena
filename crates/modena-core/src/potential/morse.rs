//! Morse pair potential: `U(r) = D·[1 − exp(−α(r−r₀))]²`.
//!
//! A bond-like pair potential with a single minimum at `r = r₀` of depth `D`
//! and curvature set by `α`. Shifted-force at the cutoff. Multi-species via
//! per-type-pair `(D, α, r₀)` tables.

use crate::error::{ModenaError, Result};
use crate::potential::check_cutoff;
use crate::potential::mixing::{MixingRule, PairParamTable};

use super::shifted_force;
use super::PairPotential;

/// Per-type Morse parameters.
#[derive(Debug, Clone, Copy)]
pub struct MorseTypeParams {
    /// Well depth `D` (energy units).
    pub d: f64,
    /// Curvature `α` (1/length).
    pub alpha: f64,
    /// Equilibrium distance `r₀` (length).
    pub r0: f64,
}

/// Explicit per-pair Morse override.
#[derive(Debug, Clone, Copy)]
pub struct MorseOverride {
    pub types: [u32; 2],
    pub d: f64,
    pub alpha: f64,
    pub r0: f64,
}

/// Morse `D·[1 − exp(−α(r−r₀))]²`, shifted-force, multi-species.
#[derive(Debug, Clone)]
pub struct Morse {
    d: PairParamTable,
    alpha: PairParamTable,
    r0: PairParamTable,
    u_rc: PairParamTable,
    f_rc: PairParamTable,
    cutoff: f64,
    cutoff2: f64,
    n_types: usize,
}

impl Morse {
    /// Single-species constructor.
    pub fn new(d: f64, alpha: f64, r0: f64, cutoff: f64) -> Result<Self> {
        #![allow(clippy::neg_cmp_op_on_partial_ord)]
        if !(d >= 0.0) {
            return Err(ModenaError::Config("Morse D must be >= 0".into()));
        }
        if !(alpha > 0.0) {
            return Err(ModenaError::Config("Morse alpha must be > 0".into()));
        }
        if !(r0 > 0.0) {
            return Err(ModenaError::Config("Morse r0 must be > 0".into()));
        }
        if !(cutoff > r0) {
            return Err(ModenaError::Config("Morse cutoff must be > r0".into()));
        }
        check_cutoff(cutoff, "Morse")?;
        Self::from_types_inner(&[MorseTypeParams { d, alpha, r0 }], &[], cutoff)
    }

    /// Multi-species constructor. Cross pairs default to the type-0 diagonal;
    /// supply overrides for real cross interactions.
    pub fn from_types(
        types: &[MorseTypeParams],
        _mixing: MixingRule,
        overrides: &[MorseOverride],
        cutoff: f64,
    ) -> Result<Self> {
        for (i, t) in types.iter().enumerate() {
            if !(t.d >= 0.0) {
                return Err(ModenaError::Config(format!(
                    "Morse D for type {i} must be >= 0"
                )));
            }
            if !(t.alpha > 0.0) {
                return Err(ModenaError::Config(format!(
                    "Morse alpha for type {i} must be > 0"
                )));
            }
            if !(t.r0 > 0.0) {
                return Err(ModenaError::Config(format!(
                    "Morse r0 for type {i} must be > 0"
                )));
            }
        }
        check_cutoff(cutoff, "Morse")?;
        Self::from_types_inner(types, overrides, cutoff)
    }

    fn from_types_inner(
        types: &[MorseTypeParams],
        overrides: &[MorseOverride],
        cutoff: f64,
    ) -> Result<Self> {
        let n = types.len();
        let cutoff2 = cutoff * cutoff;
        let mut d = PairParamTable::new(n);
        let mut alpha = PairParamTable::new(n);
        let mut r0 = PairParamTable::new(n);
        for i in 0..n {
            for j in 0..n {
                let p = &types[i.min(j)];
                d.set_symmetric(i, j, p.d);
                alpha.set_symmetric(i, j, p.alpha);
                r0.set_symmetric(i, j, p.r0);
            }
        }
        for ov in overrides {
            let [x, y] = ov.types;
            let (x, y) = (x as usize, y as usize);
            if x >= n || y >= n {
                return Err(ModenaError::Config(format!(
                    "Morse override type index {x}/{y} out of range (n_types={n})"
                )));
            }
            d.set_symmetric(x, y, ov.d);
            alpha.set_symmetric(x, y, ov.alpha);
            r0.set_symmetric(x, y, ov.r0);
        }
        let mut u_rc = PairParamTable::new(n);
        let mut f_rc = PairParamTable::new(n);
        for i in 0..n {
            for j in 0..n {
                let (u, f) = raw_energy_force(
                    cutoff,
                    d.get(i as u32, j as u32),
                    alpha.get(i as u32, j as u32),
                    r0.get(i as u32, j as u32),
                );
                u_rc.set_symmetric(i, j, u);
                f_rc.set_symmetric(i, j, f);
            }
        }
        Ok(Self {
            d,
            alpha,
            r0,
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

/// Raw (unshifted) Morse energy and radial force magnitude at `r`.
#[inline]
fn raw_energy_force(r: f64, d: f64, alpha: f64, r0: f64) -> (f64, f64) {
    let e = (-alpha * (r - r0)).exp();
    let u = d * (1.0 - e).powi(2);
    // dU/dr = 2Dα(1 − e)·e, so F_radial = −dU/dr = −2Dα(1 − e)·e.
    // Repulsive (F > 0) when r < r0 (compressed), attractive when r > r0.
    let f = -2.0 * d * alpha * (1.0 - e) * e;
    (u, f)
}

impl PairPotential for Morse {
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
            self.d.get(ti, tj),
            self.alpha.get(ti, tj),
            self.r0.get(ti, tj),
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

    fn morse() -> Morse {
        Morse::new(1.0, 2.0, 1.0, 2.5).unwrap()
    }

    #[test]
    fn minimum_at_r0_has_zero_raw_force() {
        let p = morse();
        // At r = r0 the raw force is zero; shifted force = −F(rc) (non-zero).
        let d = p.d.get(0, 0);
        let alpha = p.alpha.get(0, 0);
        let r0 = p.r0.get(0, 0);
        let (_, f_raw) = raw_energy_force(r0, d, alpha, r0);
        assert!(f_raw.abs() < 1e-12);
        // Raw energy at r0 is 0.
        let (u_raw, _) = raw_energy_force(r0, d, alpha, r0);
        assert!(u_raw.abs() < 1e-12);
    }

    #[test]
    fn energy_and_force_zero_at_cutoff() {
        let p = morse();
        let (u, f_over_r) = p.energy_force(p.cutoff2, 0, 0);
        assert!(u.abs() < 1e-9, "shifted U at cutoff = {u}");
        assert!((f_over_r * p.cutoff()).abs() < 1e-9, "shifted F at cutoff");
    }
}

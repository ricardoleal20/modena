//! Lennard-Jones 12-6 with a shifted-force-neutral cutoff.
//!
//! `U(r) = 4ε[(σ/r)¹² − (σ/r)⁶]`, energy- and force-shifted so both `U` and `F`
//! are continuous (→0) at the cutoff — the property that gives velocity-Verlet
//! its good NVE energy conservation under a truncated potential.
//!
//! Supports multi-species systems via per-type-pair `(σ, ε)` tables, combined
//! with the Lorentz–Berthelot rule by default (`σ_ij = (σ_i+σ_j)/2`,
//! `ε_ij = √(ε_i·ε_j)`), with optional explicit per-pair overrides. Single-
//! component systems use the 1-type [`LennardJones::new`] constructor, which is
//! numerically identical to the original Phase-1 kernel.

use crate::error::{ModenaError, Result};
use crate::potential::check_cutoff;
use crate::potential::mixing::{lorentz_berthelot, MixingRule, PairParamTable};

use super::shifted_force;
use super::PairPotential;

/// Per-type Lennard-Jones parameters.
#[derive(Debug, Clone, Copy)]
pub struct LjTypeParams {
    pub sigma: f64,
    pub epsilon: f64,
}

/// An explicit per-pair override of the mixing-rule default.
#[derive(Debug, Clone, Copy)]
pub struct LjOverride {
    pub types: [u32; 2],
    pub sigma: f64,
    pub epsilon: f64,
}

/// Lennard-Jones 12-6, shifted-force, multi-species.
///
/// The 1-type constructor [`new`](Self::new) reproduces the original Phase-1
/// numerics exactly (single `(σ, ε)` pair); [`from_types`](Self::from_types)
/// builds a per-type-pair table. The per-pair cutoff is uniform; the neighbor
/// list is sized to that single cutoff.
#[derive(Debug, Clone)]
pub struct LennardJones {
    /// Per-type-pair σ², indexed `ti*n_types+tj` (symmetric).
    sigma2: PairParamTable,
    /// Per-type-pair ε, indexed `ti*n_types+tj` (symmetric).
    epsilon: PairParamTable,
    /// Per-type-pair `U(rc)`, indexed `ti*n_types+tj`.
    u_rc: PairParamTable,
    /// Per-type-pair `F(rc)` (raw radial force magnitude at the cutoff).
    f_rc: PairParamTable,
    cutoff: f64,
    cutoff2: f64,
    n_types: usize,
}

impl LennardJones {
    /// Single-species constructor: one `(σ, ε)` pair. Numerically identical to
    /// the Phase-1 kernel — the 1-type table lookup is a direct index into a
    /// length-1 table.
    pub fn new(sigma: f64, epsilon: f64, cutoff: f64) -> Result<Self> {
        #![allow(clippy::neg_cmp_op_on_partial_ord)]
        if !(sigma > 0.0) {
            return Err(ModenaError::Config("LJ sigma must be > 0".into()));
        }
        if !(epsilon >= 0.0) {
            return Err(ModenaError::Config("LJ epsilon must be >= 0".into()));
        }
        check_cutoff(cutoff, "LJ")?;
        Self::from_types_inner(
            &[LjTypeParams { sigma, epsilon }],
            MixingRule::LorentzBerthelot,
            &[],
            cutoff,
        )
    }

    /// Reduced LJ benchmark preset: σ = ε = m = 1, cutoff = `cutoff_sigma`·σ.
    pub fn reduced(cutoff_sigma: f64) -> Result<Self> {
        Self::new(1.0, 1.0, cutoff_sigma)
    }

    /// Constructor that is also used via reflection-free config deserialization.
    pub fn from_params(sigma: f64, epsilon: f64, cutoff: f64) -> Result<Self> {
        Self::new(sigma, epsilon, cutoff)
    }

    /// Multi-species constructor. Builds the per-type-pair `(σ, ε)` table using
    /// the given combining rule, then applies explicit overrides. A uniform
    /// `cutoff` is used for all pairs (the standard assumption for cell-list
    /// neighbor search).
    pub fn from_types(
        types: &[LjTypeParams],
        mixing: MixingRule,
        overrides: &[LjOverride],
        cutoff: f64,
    ) -> Result<Self> {
        for (i, t) in types.iter().enumerate() {
            if !(t.sigma > 0.0) {
                return Err(ModenaError::Config(format!(
                    "LJ sigma for type {i} must be > 0"
                )));
            }
            if !(t.epsilon >= 0.0) {
                return Err(ModenaError::Config(format!(
                    "LJ epsilon for type {i} must be >= 0"
                )));
            }
        }
        check_cutoff(cutoff, "LJ")?;
        Self::from_types_inner(types, mixing, overrides, cutoff)
    }

    fn from_types_inner(
        types: &[LjTypeParams],
        _mixing: MixingRule,
        overrides: &[LjOverride],
        cutoff: f64,
    ) -> Result<Self> {
        let n = types.len();
        let cutoff2 = cutoff * cutoff;

        // Build per-pair σ and ε tables via Lorentz–Berthelot, then apply
        // explicit overrides.
        let mut sigma2 = PairParamTable::new(n);
        let mut epsilon = PairParamTable::new(n);
        for i in 0..n {
            for j in 0..n {
                let s = lorentz_berthelot::sigma(types[i].sigma, types[j].sigma);
                let e = lorentz_berthelot::epsilon(types[i].epsilon, types[j].epsilon);
                sigma2.set_symmetric(i, j, s * s);
                epsilon.set_symmetric(i, j, e);
            }
        }
        for ov in overrides {
            let [a, b] = ov.types;
            let (a, b) = (a as usize, b as usize);
            if a >= n || b >= n {
                return Err(ModenaError::Config(format!(
                    "LJ override type index {a}/{b} out of range (n_types={n})"
                )));
            }
            sigma2.set_symmetric(a, b, ov.sigma * ov.sigma);
            epsilon.set_symmetric(a, b, ov.epsilon);
        }

        // Precompute U(rc) and F(rc) per pair for the shifted-force form.
        let mut u_rc = PairParamTable::new(n);
        let mut f_rc = PairParamTable::new(n);
        for i in 0..n {
            for j in 0..n {
                let s2 = sigma2.get(i as u32, j as u32);
                let e = epsilon.get(i as u32, j as u32);
                let s = (s2 / cutoff.powi(2)).sqrt(); // σ/rc
                let s6 = s.powi(6);
                let s12 = s6 * s6;
                u_rc.set_symmetric(i, j, 4.0 * e * (s12 - s6));
                f_rc.set_symmetric(i, j, (24.0 * e / cutoff) * (2.0 * s12 - s6));
            }
        }

        Ok(Self {
            sigma2,
            epsilon,
            u_rc,
            f_rc,
            cutoff,
            cutoff2,
            n_types: n,
        })
    }

    /// Number of species types in the parameter table.
    pub fn n_types(&self) -> usize {
        self.n_types
    }
}

impl PairPotential for LennardJones {
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
        let s2 = self.sigma2.get(ti, tj) / r2;
        let s6 = s2 * s2 * s2;
        let s12 = s6 * s6;
        let e = self.epsilon.get(ti, tj);
        let u_raw = 4.0 * e * (s12 - s6);
        let f_raw = (24.0 * e / r) * (2.0 * s12 - s6);
        let u_rc = self.u_rc.get(ti, tj);
        let f_rc = self.f_rc.get(ti, tj);
        let (u, f_sf) = shifted_force::shift(u_raw, f_raw, r, self.cutoff, u_rc, f_rc);
        (u, f_sf / r)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn lj() -> LennardJones {
        LennardJones::new(1.0, 1.0, 2.5).unwrap()
    }

    #[test]
    fn energy_and_force_at_minimum() {
        // LJ minimum at r = 2^(1/6) σ: raw U = −ε, raw F = 0. The shifted-force
        // form offsets both; recompute the expected values from σ, ε, rc.
        let lj = lj();
        let sigma = 1.0;
        let epsilon = 1.0;
        let r = 2f64.powf(1.0 / 6.0) * sigma;
        let (u, f_over_r) = lj.energy_force(r * r, 0, 0);
        let rc = lj.cutoff();
        let s = sigma / rc;
        let s6c = s.powi(6);
        let s12c = s6c * s6c;
        let u_rc = 4.0 * epsilon * (s12c - s6c);
        let f_rc = (24.0 * epsilon / rc) * (2.0 * s12c - s6c);
        let sm = sigma / r;
        let s6m = sm.powi(6);
        let s12m = s6m * s6m;
        let u_expected = 4.0 * epsilon * (s12m - s6m) - u_rc + (r - rc) * f_rc;
        assert!(
            (u - u_expected).abs() < 1e-9,
            "U at min = {u}, expected {u_expected}"
        );
        // Shifted force at the minimum = −F(rc) (small, non-zero).
        assert!((f_over_r * r + f_rc).abs() < 1e-9, "F_sf at min wrong");
    }

    #[test]
    fn energy_and_force_zero_at_cutoff() {
        let lj = lj();
        let (u, f_over_r) = lj.energy_force(lj.cutoff2, 0, 0);
        assert!(u.abs() < 1e-12, "shifted U at cutoff = {u}, expected 0");
        assert!(
            (f_over_r * lj.cutoff()).abs() < 1e-12,
            "shifted F at cutoff != 0"
        );
    }

    #[test]
    fn repulsive_regime_has_positive_radial_force() {
        let lj = lj();
        let r = 0.9; // well inside σ, strongly repulsive
        let (_, f_over_r) = lj.energy_force(r * r, 0, 0);
        assert!((f_over_r * r) > 0.0);
    }

    #[test]
    fn single_type_matches_phase1_numerics() {
        // The 1-type table must produce exactly the same (U, f_over_r) as the
        // original closed-form kernel for several distances.
        let lj = LennardJones::new(1.0, 1.0, 2.5).unwrap();
        for &r in &[0.95_f64, 1.122462, 1.5, 2.0, 2.49] {
            let (u, f) = lj.energy_force(r * r, 0, 0);
            // Recompute the reference closed-form shifted-force values.
            let rc: f64 = 2.5;
            let s2 = 1.0 / (r * r);
            let s6 = s2.powi(3);
            let s12 = s6 * s6;
            let u_raw = 4.0 * (s12 - s6);
            let f_raw = (24.0 / r) * (2.0 * s12 - s6);
            let sc = 1.0 / rc;
            let s6c = sc.powi(6);
            let s12c = s6c * s6c;
            let u_rc = 4.0 * (s12c - s6c);
            let f_rc = (24.0 / rc) * (2.0 * s12c - s6c);
            let u_ref = u_raw - u_rc + (r - rc) * f_rc;
            let f_ref = (f_raw - f_rc) / r;
            assert!((u - u_ref).abs() < 1e-12, "U mismatch at r={r}");
            assert!((f - f_ref).abs() < 1e-12, "f mismatch at r={r}");
        }
    }

    #[test]
    fn multi_species_uses_cross_interaction_params() {
        // Two types: Ar (σ=1, ε=1) and a larger species (σ=2, ε=4).
        // LB cross: σ_01 = 1.5, ε_01 = 2.0.
        let lj = LennardJones::from_types(
            &[
                LjTypeParams {
                    sigma: 1.0,
                    epsilon: 1.0,
                },
                LjTypeParams {
                    sigma: 2.0,
                    epsilon: 4.0,
                },
            ],
            MixingRule::LorentzBerthelot,
            &[],
            2.5,
        )
        .unwrap();
        assert_eq!(lj.n_types(), 2);

        // Build a single-pair reference LJ with the cross params directly.
        let cross = LennardJones::new(1.5, 2.0, 2.5).unwrap();
        for &r in &[1.0, 1.5, 2.0, 2.4] {
            let (u_multi, f_multi) = lj.energy_force(r * r, 0, 1);
            let (u_cross, f_cross) = cross.energy_force(r * r, 0, 0);
            assert!((u_multi - u_cross).abs() < 1e-12, "cross U at r={r}");
            assert!((f_multi - f_cross).abs() < 1e-12, "cross f at r={r}");
        }
        // Symmetry: (0,1) == (1,0).
        let (u01, f01) = lj.energy_force(2.25, 0, 1);
        let (u10, f10) = lj.energy_force(2.25, 1, 0);
        assert!((u01 - u10).abs() < 1e-12);
        assert!((f01 - f10).abs() < 1e-12);
    }

    #[test]
    fn override_replaces_lb_cross_params() {
        let lj = LennardJones::from_types(
            &[
                LjTypeParams {
                    sigma: 1.0,
                    epsilon: 1.0,
                },
                LjTypeParams {
                    sigma: 2.0,
                    epsilon: 4.0,
                },
            ],
            MixingRule::LorentzBerthelot,
            &[LjOverride {
                types: [0, 1],
                sigma: 1.8,
                epsilon: 0.5,
            }],
            2.5,
        )
        .unwrap();
        let direct = LennardJones::new(1.8, 0.5, 2.5).unwrap();
        let (u_ov, f_ov) = lj.energy_force(2.0 * 2.0, 0, 1);
        let (u_dir, f_dir) = direct.energy_force(2.0 * 2.0, 0, 0);
        assert!((u_ov - u_dir).abs() < 1e-12);
        assert!((f_ov - f_dir).abs() < 1e-12);
    }
}

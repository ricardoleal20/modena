//! Pair potentials.
//!
//! A `PairPotential` is a short-ranged, pairwise-additive interaction. It knows
//! its cutoff and, given `r²` and the two species indices, returns the pair
//! energy `U` and the scalar `f_over_r = F_radial / r`, where `F_radial` is the
//! magnitude of the force on atom *i* along the pair vector `r_ij = r_i - r_j`.
//! The caller forms the force vector as `f_over_r * r_ij` (and its negative on
//! *j*), satisfying Newton's third law exactly.
//!
//! This is the Phase 1 interaction surface. Bonded (2/3/4-body) and
//! long-range (Coulomb/Ewald) interactions live behind the broader `ForceModel`
//! trait in `force_model.rs` and are staged in the roadmap (M2–M4).

use crate::error::{ModenaError, Result};

/// Short-ranged pairwise-additive potential.
///
/// # Conservation contract
///
/// Implementations MUST satisfy **force/energy consistency**: the energy `U(r)`
/// returned by [`energy_force`](Self::energy_force) must be a function whose
/// radial derivative is the negative of the reported force,
///
/// `dU/dr = −F(r)`, where `F(r) = f_over_r · r`.
///
/// Equivalently, the reported potential energy is the *exact integral* of the
/// reported force. This makes the NVE Hamiltonian `H = KE + PE` the true
/// conserved quantity under the symplectic integrator, so energy drift is a
/// genuine numerical signal rather than a truncation artifact. The invariant is
/// verified for every registered potential by the property-based tests in
/// `tests/property_consistency.rs` (force = −∇U by finite difference, no NaN/Inf,
/// cutoff continuity). See `conservation::ConservationReport`.
pub trait PairPotential: Send + Sync {
    /// Cutoff distance (length units).
    fn cutoff(&self) -> f64;

    /// Square of the cutoff (cached to avoid recompute in tight loops).
    fn cutoff2(&self) -> f64;

    /// Pair interaction for squared distance `r2` and species `(ti, tj)`.
    /// Returns `(U, f_over_r)` where `U` is the (shifted) pair energy and
    /// `f_over_r = F_radial / r`. Caller must skip pairs with `r2 > cutoff2`.
    ///
    /// **Contract:** `dU/dr = −(f_over_r · r)` (see the trait-level docs).
    fn energy_force(&self, r2: f64, ti: u32, tj: u32) -> (f64, f64);
}

/// Lennard-Jones 12-6 with a shifted-force-neutral cutoff.
///
/// `U(r) = 4ε[(σ/r)¹² − (σ/r)⁶]`, energy-shifted to zero at the cutoff so the
/// potential is continuous (good for NVE energy conservation). The force is the
/// exact analytical derivative:
/// `F = (24ε/r²)(2(σ/r)¹² − (σ/r)⁶) · r_vec`.
#[derive(Debug, Clone)]
pub struct LennardJones {
    pub sigma: f64,
    pub epsilon: f64,
    cutoff: f64,
    cutoff2: f64,
    sigma2: f64,
    /// `U(rc)` — potential energy at the cutoff.
    u_rc: f64,
    /// `F(rc)` — radial force magnitude at the cutoff (for shifted-force).
    f_rc: f64,
}

impl LennardJones {
    pub fn new(sigma: f64, epsilon: f64, cutoff: f64) -> Result<Self> {
        // `!(x > 0.0)` deliberately rejects NaN (NaN > 0.0 is false, so !false = true).
        // This is the intended guard — do not "simplify" to x <= 0.0 (which lets NaN through).
        #![allow(clippy::neg_cmp_op_on_partial_ord)]
        if !(sigma > 0.0) {
            return Err(ModenaError::Config("LJ sigma must be > 0".into()));
        }
        if !(cutoff > 0.0) {
            return Err(ModenaError::Config("LJ cutoff must be > 0".into()));
        }
        let cutoff2 = cutoff * cutoff;
        let sigma2 = sigma * sigma;
        // Precompute cutoff energy and force for the shifted-force form.
        let s = sigma / cutoff;
        let s6 = s.powi(6);
        let s12 = s6 * s6;
        let u_rc = 4.0 * epsilon * (s12 - s6);
        let f_rc = (24.0 * epsilon / cutoff) * (2.0 * s12 - s6);
        Ok(Self {
            sigma,
            epsilon,
            cutoff,
            cutoff2,
            sigma2,
            u_rc,
            f_rc,
        })
    }

    /// Reduced LJ benchmark preset: σ = ε = m = 1, cutoff = 2.5σ.
    pub fn reduced(cutoff_sigma: f64) -> Result<Self> {
        Self::new(1.0, 1.0, cutoff_sigma)
    }

    /// Constructor that is also used via reflection-free config deserialization.
    pub fn from_params(sigma: f64, epsilon: f64, cutoff: f64) -> Result<Self> {
        Self::new(sigma, epsilon, cutoff)
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
    fn energy_force(&self, r2: f64, _ti: u32, _tj: u32) -> (f64, f64) {
        let r = r2.sqrt();
        // (σ/r)², then ^3 and ^6.
        let s2 = self.sigma2 / r2;
        let s6 = s2 * s2 * s2;
        let s12 = s6 * s6;
        // Raw potential and radial force magnitude (positive = repulsive).
        let u_raw = 4.0 * self.epsilon * (s12 - s6);
        let f_raw = (24.0 * self.epsilon / r) * (2.0 * s12 - s6);
        // Shifted-force form: U_sf(r) = U(r) − U(rc) − (r − rc)·U'(rc),
        // with U'(rc) = −F(rc), so U_sf = U − U(rc) + (r − rc)·F(rc); and
        // F_sf(r) = F(r) − F(rc). This makes both U and F continuous (→0) at
        // the cutoff, which is what gives velocity-Verlet good energy
        // conservation under a truncated potential.
        let u = u_raw - self.u_rc + (r - self.cutoff) * self.f_rc;
        let f_sf = f_raw - self.f_rc;
        let f_over_r = f_sf / r;
        (u, f_over_r)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn lj() -> LennardJones {
        // reduced units, but with explicit shift off-by-construction it's fine
        LennardJones::new(1.0, 1.0, 2.5).unwrap()
    }

    #[test]
    fn energy_and_force_at_minimum() {
        // LJ minimum at r = 2^(1/6) σ: raw U = −ε, raw F = 0. The shifted-force
        // form offsets both; recompute the expected values from σ, ε, rc.
        let lj = lj();
        let r = 2f64.powf(1.0 / 6.0) * lj.sigma;
        let (u, f_over_r) = lj.energy_force(r * r, 0, 0);
        let rc = lj.cutoff();
        let s = lj.sigma / rc;
        let s6c = s.powi(6);
        let s12c = s6c * s6c;
        let u_rc = 4.0 * lj.epsilon * (s12c - s6c);
        let f_rc = (24.0 * lj.epsilon / rc) * (2.0 * s12c - s6c);
        let sm = lj.sigma / r;
        let s6m = sm.powi(6);
        let s12m = s6m * s6m;
        let u_expected = 4.0 * lj.epsilon * (s12m - s6m) - u_rc + (r - rc) * f_rc;
        assert!(
            (u - u_expected).abs() < 1e-9,
            "U at min = {u}, expected {u_expected}"
        );
        // Shifted force at the minimum = −F(rc) (small, non-zero).
        assert!((f_over_r * r + f_rc).abs() < 1e-9, "F_sf at min wrong");
    }

    #[test]
    fn energy_and_force_zero_at_cutoff() {
        // Shifted-force: both U and F vanish at the cutoff.
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
        // F_radial = f_over_r * r > 0 (pushes i away from j)
        assert![(f_over_r * r) > 0.0];
    }
}

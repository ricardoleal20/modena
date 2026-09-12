//! Pair potentials.
//!
//! A `PairPotential` is a short-ranged, pairwise-additive interaction. It knows
//! its cutoff and, given `r²` and the two species indices, returns the pair
//! energy `U` and the scalar `f_over_r = F_radial / r`, where `F_radial` is the
//! magnitude of the force on atom *i* along the pair vector `r_ij = r_i - r_j`.
//! The caller forms the force vector as `f_over_r * r_ij` (and its negative on
//! *j*), satisfying Newton's third law exactly.
//!
//! M2 ships four analytic pair potentials — Lennard-Jones, Buckingham, Morse,
//! and soft-sphere — plus a tabulated potential, all in the shifted-force form
//! so both `U` and `F` are continuous (→0) at the cutoff. Multi-species systems
//! use per-type-pair parameter tables (Lorentz–Berthelot mixing by default).

pub mod buckingham;
pub mod lennard_jones;
pub mod mixing;
pub mod morse;
pub mod soft_sphere;
pub mod tabulated;

pub use buckingham::{Buckingham, BuckinghamOverride, BuckinghamTypeParams};
pub use lennard_jones::{LennardJones, LjOverride, LjTypeParams};
pub use mixing::{MixingRule, PairParamTable};
pub use morse::{Morse, MorseOverride, MorseTypeParams};
pub use soft_sphere::{SoftSphere, SoftSphereOverride, SoftSphereTypeParams};
pub use tabulated::{InterpKind, Tabulated};

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
///
/// # Species indices
///
/// `ti` and `tj` index per-type parameters in the potential's parameter table
/// (see [`PairParamTable`]). For single-component systems both are `0`. The
/// pair index is symmetric: `(ti, tj)` and `(tj, ti)` return identical values.
pub trait PairPotential: Send + Sync {
    /// Cutoff distance (length units). This is the *maximum* cutoff across all
    /// type pairs, used to size the neighbor list.
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

/// Shared shifted-force helpers.
///
/// The shifted-force form makes both `U` and `F` continuous at the cutoff:
///   `U_sf(r) = U(r) − U(rc) + (r − rc)·F(rc)`
///   `F_sf(r) = F(r) − F(rc)`
/// where `F(rc)` is the *raw* radial force magnitude at the cutoff. With
/// `f_over_r = F_sf / r`, this gives `dU_sf/dr = −F_sf` exactly.
pub(crate) mod shifted_force {
    /// Build the shifted pair `(U, F_radial)` from the raw values and the
    /// cutoff quantities. Returns `(U_sf, f_over_r)` where `F_radial_sf = F_raw − F_rc`.
    #[inline]
    pub fn shift(u_raw: f64, f_raw: f64, r: f64, cutoff: f64, u_rc: f64, f_rc: f64) -> (f64, f64) {
        let u = u_raw - u_rc + (r - cutoff) * f_rc;
        let f_sf = f_raw - f_rc;
        (u, f_sf)
    }
}

/// Validate that a cutoff is positive (NaN-safe, matching the LJ guard style).
pub(crate) fn check_cutoff(cutoff: f64, name: &str) -> Result<f64> {
    #![allow(clippy::neg_cmp_op_on_partial_ord)]
    if !(cutoff > 0.0) {
        return Err(ModenaError::Config(format!("{name} cutoff must be > 0")));
    }
    Ok(cutoff)
}

#[cfg(test)]
mod tests {
    // The per-potential analytic tests live alongside each impl. Cross-cutting
    // property tests (force = −∇U, cutoff continuity, no NaN/Inf) for every
    // registered potential live in tests/property_consistency.rs.
}

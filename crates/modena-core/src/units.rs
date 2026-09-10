//! Unit handling.
//!
//! **Phase 1** runs in *reduced (LJ) units*: σ = ε = m = 1, time in √(mσ²/ε),
//! temperature in ε/kB. This keeps the PoC verifiable against the canonical LJ
//! benchmarks and removes unit-conversion bugs from the correctness path.
//!
//! Roadmap (M2): a `UnitSystem` enum (`Reduced`, `Real { length, energy, mass,
//! temperature }`) with conversion factors applied at the I/O boundary, so the
//! core always computes in reduced units internally.

/// The unit regime the simulation runs in.
#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
#[serde(rename_all = "snake_case")]
#[derive(Default)]
pub enum UnitSystem {
    /// Reduced Lennard-Jones units (σ = ε = m = 1).
    #[default]
    Reduced,
    /// Real units (Å, kcal/mol, amu, K). Conversion applied at I/O boundary —
    /// not yet wired in Phase 1; specifying `Real` returns a config error today.
    Real,
}

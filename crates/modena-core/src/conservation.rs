//! Energy-conservation contract.
//!
//! MODENA's headline correctness guarantee is that the **reported potential
//! energy is the exact integral of the reported force**: for every
//! `PairPotential`, the energy `U(r)` returned by `energy_force` satisfies
//! `dU/dr = −F(r)`, where `F(r) = f_over_r · r` is the radial force magnitude.
//! Because the integrator (velocity-Verlet) is symplectic, this makes the NVE
//! Hamiltonian `H = KE + PE` the true conserved quantity — so a rising energy
//! drift is a *signal* (numerical breakdown, bad timestep, phase transition),
//! not an expected artifact of the potential form.
//!
//! This invariant is *already* satisfied by `LennardJones` (shifted-force,
//! analytic). The contract is made *enforced* rather than incidental by the
//! property-based tests in `tests/property_consistency.rs`, and surfaced to
//! callers (the agent interface, the GUI) as a structured [`ConservationReport`].

/// Verdict on an NVE run's energy conservation, derived from the relative drift.
///
/// Thresholds are calibrated for the LJ liquid at `dt = 0.005` (reduced units).
/// They are advisory labels, not hard gates — the caller decides how to react.
#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize)]
#[serde(rename_all = "snake_case")]
pub enum ConservationVerdict {
    /// Relative drift `< 1e-4`. The run is numerically sound.
    Healthy,
    /// Relative drift `< 1e-3`. Likely a too-large timestep or a stiff
    /// configuration; investigate before trusting derived quantities.
    Warning,
    /// Relative drift `>= 1e-3`. The integrator is no longer conserving the
    /// Hamiltonian — results from this run are not reliable.
    Drifting,
}

impl ConservationVerdict {
    /// Classify a relative energy drift into a verdict.
    pub fn from_drift(relative_drift: f64) -> Self {
        let d = relative_drift.abs();
        if d < 1e-4 {
            Self::Healthy
        } else if d < 1e-3 {
            Self::Warning
        } else {
            Self::Drifting
        }
    }
}

/// A structured energy-conservation report for a finished (or in-flight) run.
///
/// This is the typed signal an LLM agent reasons over ("drift 4.7e-6 → healthy,
/// transition evidence is trustworthy") and that the GUI's conservation gauge
/// displays. Drift is *data*, not an error: a degrading verdict mid-sweep can
/// flag a phase transition or numerical breakdown.
#[derive(Debug, Clone, Copy, PartialEq, serde::Serialize)]
pub struct ConservationReport {
    /// Total energy at step 0.
    pub e0: f64,
    /// Total energy at the last sampled step.
    pub e_final: f64,
    /// `(E_final − E0) / |E0|` (signed).
    pub relative_drift: f64,
    /// Classified verdict.
    pub verdict: ConservationVerdict,
}

impl ConservationReport {
    /// Build a report from the initial and final total energies.
    pub fn new(e0: f64, e_final: f64) -> Self {
        let denom = e0.abs().max(1e-12);
        let relative_drift = (e_final - e0) / denom;
        let verdict = ConservationVerdict::from_drift(relative_drift);
        Self {
            e0,
            e_final,
            relative_drift,
            verdict,
        }
    }

    /// A report for a run that has not yet produced a final energy (e.g. a
    /// freshly started async run). Verdict is `Healthy` by convention.
    pub fn pending(e0: f64) -> Self {
        Self::new(e0, e0)
    }

    /// Update an in-flight report with a newly observed total energy.
    pub fn observe(&mut self, e: f64) {
        self.e_final = e;
        let denom = self.e0.abs().max(1e-12);
        self.relative_drift = (e - self.e0) / denom;
        self.verdict = ConservationVerdict::from_drift(self.relative_drift);
    }
}

impl Default for ConservationReport {
    fn default() -> Self {
        Self::pending(0.0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn verdict_thresholds() {
        assert_eq!(
            ConservationVerdict::from_drift(4.7e-6),
            ConservationVerdict::Healthy
        );
        assert_eq!(
            ConservationVerdict::from_drift(-3e-5),
            ConservationVerdict::Healthy
        );
        assert_eq!(
            ConservationVerdict::from_drift(5e-4),
            ConservationVerdict::Warning
        );
        assert_eq!(
            ConservationVerdict::from_drift(2e-3),
            ConservationVerdict::Drifting
        );
        assert_eq!(
            ConservationVerdict::from_drift(-1e-2),
            ConservationVerdict::Drifting
        );
    }

    #[test]
    fn report_signs_and_observe() {
        let mut r = ConservationReport::new(-3622.0, -3622.0);
        assert!(r.relative_drift.abs() < 1e-9);
        assert_eq!(r.verdict, ConservationVerdict::Healthy);

        r.observe(-3622.5); // ~1.4e-4 drift
        assert_eq!(r.verdict, ConservationVerdict::Warning);
        assert!(r.relative_drift < 0.0);
    }

    #[test]
    fn report_handles_near_zero_e0() {
        // Must not panic or produce inf when E0 ≈ 0.
        let r = ConservationReport::new(0.0, 1.0);
        assert!(r.relative_drift.is_finite());
    }
}

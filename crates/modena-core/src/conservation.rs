//! Energy diagnostics, distinct from validation of forces, sampling or materials.
//!
//! The runner uses [`EnergyMonitor`] to sample every integration step in reduced
//! time, independently of output cadence. It normalizes by the initial kinetic
//! energy, which is invariant under adding a constant to potential energy.
//! Zero initial kinetic energy leaves normalized quantities unavailable.
//!
//! Energy error, its bounded fluctuations and its regression slope are separate
//! signals. A symplectic scheme has discretization error; neither a small error
//! nor a force-energy gradient identity certifies the sampled ensemble. Physical
//! total energy is not conserved under thermal or pressure coupling.

use crate::error::{ModenaError, Result};

/// Advisory classification of an energy-error metric, not a validity certificate.
///
/// The legacy endpoint classifier retains 1e-4/1e-3. Runner diagnostics require
/// explicit thresholds for maximum excursion / initial kinetic energy: the old
/// cut points are not transferred to this different metric. Inspect timestep
/// convergence, the band and slope rather than interpreting a label as proof.
#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize)]
#[serde(rename_all = "snake_case")]
pub enum ConservationVerdict {
    /// The supplied error metric is below the selected warning threshold.
    Healthy,
    /// The supplied error metric is between the selected thresholds.
    Warning,
    /// Error at or above the selected drifting threshold. This legacy name does not establish secular drift;
    /// even a bounded oscillation can have a large excursion.
    Drifting,
    /// Energy exchange is expected; no NVE verdict is applicable.
    NotApplicable,
    /// A physical reference scale or finite error metric is unavailable.
    Unavailable,
    /// No integration interval has yet been observed.
    InsufficientData,
    /// A normalized metric exists, but the caller has not selected tolerances.
    Uncalibrated,
}

impl ConservationVerdict {
    /// Classify a relative energy drift into a verdict.
    pub fn from_drift(relative_drift: f64) -> Self {
        if !relative_drift.is_finite() {
            return Self::Unavailable;
        }
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

/// Whether the physical total energy can be assessed as an NVE invariant.
#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize)]
#[serde(rename_all = "snake_case")]
pub enum ConservationApplicability {
    Nve,
    EnergyExchange,
    /// Legacy endpoint constructors have no ensemble information.
    Unspecified,
}

/// Normalization used by the runner, independent of the arbitrary zero of U.
#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize)]
#[serde(rename_all = "snake_case")]
pub enum EnergyNormalization {
    InitialKineticEnergy,
}

/// Caller-selected dimensionless advisory thresholds for max|E(t)-E0| / KE0.
/// These are policy inputs, not proof that a material or timestep is validated.
#[derive(Debug, Clone, Copy, PartialEq, serde::Serialize, serde::Deserialize)]
#[serde(deny_unknown_fields)]
pub struct EnergyErrorThresholds {
    pub warning: f64,
    pub drifting: f64,
}

impl EnergyErrorThresholds {
    pub fn validate(&self) -> Result<()> {
        if !self.warning.is_finite()
            || !self.drifting.is_finite()
            || self.warning <= 0.0
            || self.drifting <= self.warning
        {
            return Err(ModenaError::Config(
                "energy_error_thresholds require finite 0 < warning < drifting".into(),
            ));
        }
        Ok(())
    }

    fn classify(&self, metric: f64) -> ConservationVerdict {
        if metric < self.warning {
            ConservationVerdict::Healthy
        } else if metric < self.drifting {
            ConservationVerdict::Warning
        } else {
            ConservationVerdict::Drifting
        }
    }
}

/// Statistics over the initial state and every step of this run invocation.
/// Energies and time are reduced, even when thermodynamic output uses real units.
#[derive(Debug, Clone, Copy, PartialEq, serde::Serialize)]
pub struct EnergyDiagnostics {
    pub samples: usize,
    pub duration: f64,
    pub energy_units: &'static str,
    pub time_units: &'static str,
    pub delta_energy: f64,
    pub delta_per_atom: f64,
    pub max_abs_deviation: f64,
    pub energy_range: f64,
    /// Population RMS around the sample mean, not endpoint drift or RMS from E0.
    pub rms_fluctuation: f64,
    /// OLS slope in reduced energy / reduced time. Not a significance test.
    pub slope: Option<f64>,
    pub normalization: EnergyNormalization,
    pub reference_energy: Option<f64>,
    pub normalized_change: Option<f64>,
    pub normalized_max_abs_deviation: Option<f64>,
    pub advisory_thresholds: Option<EnergyErrorThresholds>,
}

/// Energy changes plus an explicitly scoped advisory and optional time series.
#[derive(Debug, Clone, Copy, PartialEq, serde::Serialize)]
pub struct ConservationReport {
    /// Total energy at step 0.
    pub e0: f64,
    /// Total energy at the last sampled step.
    pub e_final: f64,
    /// Legacy endpoint change divided by `max(|E0|,1e-12)`. Retained for numeric
    /// compatibility; depends on the zero of U and is not the runner's verdict basis.
    pub relative_drift: f64,
    /// Legacy constructors classify relative_drift; the runner classifies the
    /// normalized maximum excursion in diagnostics, only for NVE.
    pub verdict: ConservationVerdict,
    pub applicability: ConservationApplicability,
    pub diagnostics: Option<EnergyDiagnostics>,
}

impl ConservationReport {
    /// Legacy endpoint report without an ensemble, time or physical scale.
    /// Prefer [`EnergyMonitor`] for scientifically scoped diagnostics.
    pub fn new(e0: f64, e_final: f64) -> Self {
        let denom = e0.abs().max(1e-12);
        let relative_drift = (e_final - e0) / denom;
        let verdict = ConservationVerdict::from_drift(relative_drift);
        Self {
            e0,
            e_final,
            relative_drift,
            verdict,
            applicability: ConservationApplicability::Unspecified,
            diagnostics: None,
        }
    }

    /// A report for a run that has not yet produced a final energy (e.g. a
    /// freshly started legacy calculation). No ensemble or scale is inferred.
    pub fn pending(e0: f64) -> Self {
        Self::new(e0, e0)
    }

    /// Update a legacy endpoint report. This clears contextual statistics;
    /// runners and live consumers should use EnergyMonitor / runner events.
    pub fn observe(&mut self, e: f64) {
        *self = Self::new(self.e0, e);
    }
}

impl Default for ConservationReport {
    fn default() -> Self {
        Self {
            verdict: ConservationVerdict::InsufficientData,
            ..Self::pending(0.0)
        }
    }
}

/// Constant-memory, centered online energy/time statistics. A new monitor starts
/// at t=0 for each invocation, including continuation from a checkpoint.
#[derive(Debug)]
pub struct EnergyMonitor {
    e0: f64,
    e_final: f64,
    atoms: usize,
    reference_energy: Option<f64>,
    applicability: ConservationApplicability,
    samples: usize,
    last_time: f64,
    mean_delta: f64,
    m2_delta: f64,
    min_delta: f64,
    max_delta: f64,
    mean_time: f64,
    m2_time: f64,
    covariance: f64,
    thresholds: Option<EnergyErrorThresholds>,
}

impl EnergyMonitor {
    /// Construct from reduced initial total/kinetic energy and a positive atom count.
    pub fn new(
        e0: f64,
        atoms: usize,
        initial_kinetic_energy: f64,
        applicability: ConservationApplicability,
    ) -> Result<Self> {
        if atoms == 0
            || !e0.is_finite()
            || !initial_kinetic_energy.is_finite()
            || initial_kinetic_energy < 0.0
        {
            return Err(diagnostic_error(
                "initial energy/kinetic scale must be finite, kinetic >= 0 and atoms > 0",
            ));
        }
        Ok(Self {
            e0,
            e_final: e0,
            atoms,
            reference_energy: (initial_kinetic_energy > 0.0).then_some(initial_kinetic_energy),
            applicability,
            samples: 1,
            last_time: 0.0,
            mean_delta: 0.0,
            m2_delta: 0.0,
            min_delta: 0.0,
            max_delta: 0.0,
            mean_time: 0.0,
            m2_time: 0.0,
            covariance: 0.0,
            thresholds: None,
        })
    }

    /// Enable advisory classification with explicitly selected tolerances.
    pub fn with_thresholds(mut self, thresholds: EnergyErrorThresholds) -> Result<Self> {
        thresholds.validate()?;
        self.thresholds = Some(thresholds);
        Ok(self)
    }

    /// Observe at a strictly increasing reduced time relative to this invocation.
    /// Invalid or overflowing samples leave the previous statistics intact.
    pub fn observe(&mut self, time: f64, energy: f64) -> Result<()> {
        if !energy.is_finite() {
            return Err(diagnostic_error(&format!(
                "non-finite energy at reduced time {time}: {energy}"
            )));
        }
        if !time.is_finite() || time <= self.last_time {
            return Err(diagnostic_error(
                "samples need finite energy and increasing finite time",
            ));
        }
        let samples = self
            .samples
            .checked_add(1)
            .ok_or_else(|| diagnostic_error("sample count overflow"))?;
        // Center at E0 before Welford updates: an energy offset must not dominate
        // variance or the time covariance through cancellation of large sums.
        let delta = energy - self.e0;
        let de = delta - self.mean_delta;
        let dt = time - self.mean_time;
        let mean_delta = self.mean_delta + de / samples as f64;
        let mean_time = self.mean_time + dt / samples as f64;
        let m2_delta = self.m2_delta + de * (delta - mean_delta);
        let m2_time = self.m2_time + dt * (time - mean_time);
        let covariance = self.covariance + dt * (delta - mean_delta);
        let min_delta = self.min_delta.min(delta);
        let max_delta = self.max_delta.max(delta);
        if ![
            delta,
            mean_delta,
            mean_time,
            m2_delta,
            m2_time,
            covariance,
            max_delta - min_delta,
        ]
        .iter()
        .all(|x| x.is_finite())
        {
            return Err(diagnostic_error("energy/time statistics overflow"));
        }
        self.e_final = energy;
        self.samples = samples;
        self.last_time = time;
        self.mean_delta = mean_delta;
        self.mean_time = mean_time;
        self.m2_delta = m2_delta;
        self.m2_time = m2_time;
        self.covariance = covariance;
        self.min_delta = min_delta;
        self.max_delta = max_delta;
        Ok(())
    }

    /// Snapshot without discarding intermediate energy errors or changing the monitor.
    pub fn report(&self) -> ConservationReport {
        let delta_energy = self.e_final - self.e0;
        let max_abs_deviation = self.min_delta.abs().max(self.max_delta.abs());
        let normalize = |value: f64| {
            self.reference_energy
                .and_then(|scale| finite(value / scale))
        };
        let normalized_max_abs_deviation = normalize(max_abs_deviation);
        let verdict = match self.applicability {
            ConservationApplicability::EnergyExchange => ConservationVerdict::NotApplicable,
            ConservationApplicability::Unspecified => ConservationVerdict::Unavailable,
            ConservationApplicability::Nve if self.samples < 2 => {
                ConservationVerdict::InsufficientData
            }
            ConservationApplicability::Nve => match (normalized_max_abs_deviation, self.thresholds)
            {
                (None, _) => ConservationVerdict::Unavailable,
                (Some(_), None) => ConservationVerdict::Uncalibrated,
                (Some(metric), Some(thresholds)) => thresholds.classify(metric),
            },
        };
        ConservationReport {
            verdict,
            applicability: self.applicability,
            diagnostics: Some(EnergyDiagnostics {
                samples: self.samples,
                duration: self.last_time,
                energy_units: "reduced",
                time_units: "reduced",
                delta_energy,
                delta_per_atom: delta_energy / self.atoms as f64,
                max_abs_deviation,
                energy_range: self.max_delta - self.min_delta,
                rms_fluctuation: (self.m2_delta.max(0.0) / self.samples as f64).sqrt(),
                slope: if self.m2_time > 0.0 {
                    finite(self.covariance / self.m2_time)
                } else {
                    None
                },
                normalization: EnergyNormalization::InitialKineticEnergy,
                reference_energy: self.reference_energy,
                normalized_change: normalize(delta_energy),
                normalized_max_abs_deviation,
                advisory_thresholds: self.thresholds,
            }),
            ..ConservationReport::new(self.e0, self.e_final)
        }
    }
}

fn finite(value: f64) -> Option<f64> {
    value.is_finite().then_some(value)
}

fn diagnostic_error(message: &str) -> ModenaError {
    ModenaError::Physics(format!("energy diagnostics: {message}"))
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

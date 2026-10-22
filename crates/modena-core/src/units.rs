//! Unit handling.
//!
//! The core always computes in *reduced (LJ) units*: σ = ε = m = 1, time in
//! √(mσ²/ε), temperature in ε/kB. This keeps the physics verifiable against the
//! canonical LJ benchmarks and removes unit-conversion bugs from the
//! correctness path.
//!
//! M2 adds a `UnitSystem::Real` regime: the user writes `units = "real"` and
//! supplies physical parameters (Å, kcal/mol, amu, K). A [`UnitFactors`] built
//! from a reference (σ_ref, ε_ref, m_ref) converts every numeric field to
//! reduced units at the I/O boundary (`InputConfig::from_path`), so the core
//! never sees real units. Output (CSV logs, trajectory, the CLI summary) is
//! converted back to real units for the user.

/// The unit regime the simulation runs in.
#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
#[serde(rename_all = "snake_case")]
#[derive(Default)]
pub enum UnitSystem {
    /// Reduced Lennard-Jones units (σ = ε = m = 1).
    #[default]
    Reduced,
    /// Real units (Å, kcal/mol, amu, K). Conversion applied at the I/O
    /// boundary via [`UnitFactors`]; the core stays reduced internally.
    Real,
}

/// Conversion factors between real and reduced units, derived from a reference
/// length (σ_ref, Å), energy (ε_ref, kcal/mol), and mass (m_ref, amu).
///
/// Reduced = real / factor for each dimension. To go back to real units,
/// multiply the reduced value by the factor.
#[derive(Debug, Clone, Copy)]
pub struct UnitFactors {
    /// Reference length σ (Å). `length_reduced = length_Å / sigma`.
    pub sigma: f64,
    /// Reference energy ε (kcal/mol). `energy_reduced = energy_kcal / epsilon`.
    pub epsilon: f64,
    /// Reference mass m (amu). `mass_reduced = mass_amu / mass_ref`.
    pub mass_ref: f64,
    /// Time factor: `t_reduced = t_ps / t_factor` where
    /// `t_factor = √(m_ref·σ²/ε)` in ps.
    pub time: f64,
    /// Temperature factor: `T_reduced = T_K / T_factor` where
    /// `T_factor = ε / kB` and `kB = 0.001987204259 kcal/(mol·K)`.
    pub temperature: f64,
    /// Pressure factor: `P_reduced = P_atm / P_factor` where
    /// `P_factor = ε/σ³` converted to atm.
    pub pressure: f64,
    /// Charge factor: `q_reduced = q_e · charge` where
    /// `charge = √(COULOMB_KCAL_ANGSTROM / (σ · ε))` so that, in reduced units,
    /// `U = q_reduced_i · q_reduced_j / r_reduced`. Charges are input in units
    /// of the elementary charge `e`.
    pub charge: f64,
}

/// Boltzmann constant in kcal/(mol·K).
pub const KB_KCAL: f64 = 0.001987204259;

/// Conversion: 1 kcal/mol per Å³ → atm. `1 kcal/(mol·Å³) = 69477.8 atm`.
pub const KCAL_PER_A3_TO_ATM: f64 = 69477.8;

impl UnitFactors {
    /// Build the conversion factors from reference σ (Å), ε (kcal/mol), m (amu).
    pub fn new(sigma: f64, epsilon: f64, mass_ref: f64) -> Self {
        // Reduced time unit τ = √(m_ref·σ²/ε_ref) expressed in picoseconds.
        // Working in amu/Å/(kcal·mol⁻¹):
        //   τ = √(amu·Å²/(kcal/mol)) · √(1.66054e-27 kg/amu · 1e-20 m²/Å² · N_A / (4184 J/kcal))
        //     = √(m σ²/ε) · 4.8888e-14 s  =  √(m σ²/ε) · 0.048888 ps.
        // For argon (σ=3.40, ε=0.238, m=39.948) this gives τ ≈ 2.15 ps, the
        // canonical LJ reduced-time value.
        const TIME_UNIT_PS: f64 = 0.0488882;
        let time = TIME_UNIT_PS * (mass_ref * sigma * sigma / epsilon).sqrt();
        let temperature = epsilon / KB_KCAL; // K per reduced unit
        let pressure = epsilon / (sigma * sigma * sigma) * KCAL_PER_A3_TO_ATM; // atm per reduced unit
                                                                               // Reduced-charge factor: with real Coulomb U = COULOMB · q_i q_j / r
                                                                               // (kcal/mol, r in Å, q in e), writing U/ε = (q·charge)² / (r/σ) requires
                                                                               // q_reduced² = q² · COULOMB / (σ · ε), i.e. charge = √(COULOMB / (σ · ε)).
        const COULOMB_KCAL_ANGSTROM: f64 = 332.06371;
        let charge = (COULOMB_KCAL_ANGSTROM / (sigma * epsilon)).sqrt();
        Self {
            sigma,
            epsilon,
            mass_ref,
            time,
            temperature,
            pressure,
            charge,
        }
    }

    /// Convert a real length (Å) to reduced.
    #[inline]
    pub fn to_reduced_length(&self, l: f64) -> f64 {
        l / self.sigma
    }

    /// Convert a real energy (kcal/mol) to reduced.
    #[inline]
    pub fn to_reduced_energy(&self, e: f64) -> f64 {
        e / self.epsilon
    }

    /// Convert a real mass (amu) to reduced.
    #[inline]
    pub fn to_reduced_mass(&self, m: f64) -> f64 {
        m / self.mass_ref
    }

    /// Convert a real temperature (K) to reduced.
    #[inline]
    pub fn to_reduced_temperature(&self, t: f64) -> f64 {
        t / self.temperature
    }

    /// Convert a real time (ps) to reduced.
    #[inline]
    pub fn to_reduced_time(&self, t: f64) -> f64 {
        t / self.time
    }

    /// Convert a reduced energy back to real (kcal/mol).
    #[inline]
    pub fn to_real_energy(&self, e: f64) -> f64 {
        e * self.epsilon
    }

    /// Convert a reduced temperature back to real (K).
    #[inline]
    pub fn to_real_temperature(&self, t: f64) -> f64 {
        t * self.temperature
    }

    /// Convert a reduced pressure back to real (atm).
    #[inline]
    pub fn to_real_pressure(&self, p: f64) -> f64 {
        p * self.pressure
    }

    /// Convert a reduced time back to real (ps).
    #[inline]
    pub fn to_real_time(&self, t: f64) -> f64 {
        t * self.time
    }

    /// Convert a real pressure (atm) to reduced.
    #[inline]
    pub fn to_reduced_pressure(&self, p: f64) -> f64 {
        p / self.pressure
    }

    /// Convert a real charge (in units of `e`) to reduced. The reduced Coulomb
    /// energy is then `U = q_red_i · q_red_j / r_red` (multiplicative, since
    /// charge is a composite unit: `q_red = q_e · charge`).
    #[inline]
    pub fn to_reduced_charge(&self, q: f64) -> f64 {
        q * self.charge
    }

    /// Convert a reduced length back to real (Å).
    #[inline]
    pub fn to_real_length(&self, l: f64) -> f64 {
        l * self.sigma
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn argon_factors_sensible() {
        // Argon: σ = 3.40 Å, ε = 0.238 kcal/mol, m = 39.948 amu.
        let f = UnitFactors::new(3.40, 0.238, 39.948);
        // T_factor = ε/kB ≈ 0.238 / 0.001987 ≈ 119.8 K. The textbook ε/kB for
        // argon is ~120 K (critical T ≈ 151 K).
        assert!(
            (f.temperature - 119.8).abs() < 1.0,
            "T_factor = {}",
            f.temperature
        );
        // Reduced time unit τ = √(mσ²/ε) ≈ 2.15 ps — the canonical LJ value.
        assert!((f.time - 2.15).abs() < 0.1, "t_factor = {}", f.time);
    }

    #[test]
    fn reduced_round_trip_is_identity() {
        let f = UnitFactors::new(3.40, 0.238, 39.948);
        let l = 5.26; // Å
        let r = f.to_reduced_length(l);
        assert!((f.to_real_length(r) - l).abs() < 1e-12);
        let e = 0.5; // kcal/mol
        assert!((f.to_real_energy(f.to_reduced_energy(e)) - e).abs() < 1e-12);
    }
}

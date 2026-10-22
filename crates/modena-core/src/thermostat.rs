//! Thermostats — velocity-only temperature control that composes with the
//! `VelocityVerlet` integrator.
//!
//! For stochastic / genuinely integrator-coupled thermostats (Langevin BAOAB),
//! see `integrator::LangevinBAOAB`, which is a self-contained integrator rather
//! than a thermostat. The split: a `Thermostat` rescales velocities given the
//! current KE; a Langevin *integrator* interleaves noise/friction with the
//! force evaluations, so it cannot be expressed as a velocity rescale.

use chacha20::ChaCha12Rng;
use rand::{RngExt, SeedableRng};
use rand_distr::{Distribution, StandardNormal};
use serde::{Deserialize, Serialize};
use serde_json::{json, Value};

use crate::error::{ModenaError, Result};
use crate::system::System;

/// Velocity rescaling thermostat, applied once per step.
pub trait Thermostat: Send {
    fn apply(&mut self, system: &mut System, dt: f64);
    /// Stateless thermostats encode null; stateful implementations override.
    fn checkpoint_state(&self) -> Result<Value> {
        Ok(Value::Null)
    }
    fn restore_state(&mut self, state: &Value) -> Result<()> {
        restore_stateless(state)
    }
}

pub(crate) fn state_error(message: &str) -> ModenaError {
    ModenaError::Input(format!("checkpoint state: {message}"))
}
pub(crate) fn restore_stateless(state: &Value) -> Result<()> {
    if state.is_null() {
        Ok(())
    } else {
        Err(state_error("stateless component requires null"))
    }
}

#[derive(Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct RngState {
    algorithm: String,
    seed: [u8; 32],
    stream: u64,
    // A ChaCha word offset is 68 bits; a decimal string avoids JSON's u64
    // limit and prevents conversion through an imprecise floating number.
    word_pos: String,
}

pub(crate) fn checkpoint_rng(rng: &ChaCha12Rng) -> Result<Value> {
    serde_json::to_value(RngState {
        algorithm: "chacha12-v1".into(),
        seed: rng.get_seed(),
        stream: rng.get_stream(),
        word_pos: rng.get_word_pos().to_string(),
    })
    .map_err(|err| state_error(&err.to_string()))
}
pub(crate) fn restore_rng(value: &Value, expected_seed: [u8; 32]) -> Result<ChaCha12Rng> {
    let state: RngState =
        serde_json::from_value(value.clone()).map_err(|err| state_error(&err.to_string()))?;
    let pos = state
        .word_pos
        .parse::<u128>()
        .map_err(|_| state_error("invalid RNG word position"))?;
    if state.algorithm != "chacha12-v1" || pos >= (1_u128 << 68) || state.seed != expected_seed {
        return Err(state_error("RNG kind, seed or word position is invalid"));
    }
    let mut rng = ChaCha12Rng::from_seed(state.seed);
    rng.set_stream(state.stream);
    rng.set_word_pos(pos);
    Ok(rng)
}
pub(crate) fn validate_state_parameters(
    state: &Value,
    expected: &Value,
    keys: &[&str],
) -> Result<()> {
    if !state.is_object() || keys.iter().any(|key| state.get(*key) != expected.get(*key)) {
        Err(state_error("component kind/parameters differ"))
    } else {
        Ok(())
    }
}

/// No thermostat — pure NVE.
pub struct NoneThermostat;

impl Thermostat for NoneThermostat {
    fn apply(&mut self, _system: &mut System, _dt: f64) {}
}

/// Berendsen weak-coupling thermostat. Velocity-rescales every step toward
/// `target_temperature` with time constant `tau`. Deterministic, simple, good
/// for equilibration; does *not* sample the exact canonical ensemble (use
/// Langevin or Nosé–Hoover for correct NVT sampling — roadmap M5).
pub struct Berendsen {
    pub target_temperature: f64,
    pub tau: f64,
}

impl Berendsen {
    pub fn new(target_temperature: f64, tau: f64) -> Self {
        Self {
            target_temperature,
            tau,
        }
    }
}

impl Thermostat for Berendsen {
    fn checkpoint_state(&self) -> Result<Value> {
        if !self.target_temperature.is_finite() || !self.tau.is_finite() {
            return Err(state_error("nonfinite Berendsen parameters"));
        }
        Ok(
            json!({"kind":"berendsen","version":1,"target_temperature":self.target_temperature,"tau":self.tau}),
        )
    }
    fn restore_state(&mut self, state: &Value) -> Result<()> {
        if state == &self.checkpoint_state()? {
            Ok(())
        } else {
            Err(state_error("Berendsen kind/parameters differ"))
        }
    }
    fn apply(&mut self, system: &mut System, dt: f64) {
        let t = system.temperature();
        if t <= 0.0 || self.tau <= 0.0 {
            return;
        }
        let lambda = (1.0 + (dt / self.tau) * (self.target_temperature / t - 1.0)).sqrt();
        for v in system.velocities.iter_mut() {
            *v *= lambda;
        }
    }
}

/// Stochastic velocity rescaling (Bussi–Parrinello 2007). Samples the exact
/// canonical kinetic-energy distribution; deterministic given the seed.
pub struct VelocityRescaling {
    pub target_temperature: f64,
    pub tau: f64,
    rng: ChaCha12Rng,
    seed: u64,
}

impl VelocityRescaling {
    pub fn new(target_temperature: f64, tau: f64, seed: u64) -> Self {
        Self {
            target_temperature,
            tau,
            rng: ChaCha12Rng::seed_from_u64(seed),
            seed,
        }
    }
}

impl Thermostat for VelocityRescaling {
    fn checkpoint_state(&self) -> Result<Value> {
        if !self.target_temperature.is_finite() || !self.tau.is_finite() {
            return Err(state_error("nonfinite CSVR parameters"));
        }
        Ok(
            json!({"kind":"csvr","version":1,"target_temperature":self.target_temperature,
            "tau":self.tau,"seed":self.seed,"rng":checkpoint_rng(&self.rng)?}),
        )
    }
    fn restore_state(&mut self, state: &Value) -> Result<()> {
        validate_state_parameters(
            state,
            &self.checkpoint_state()?,
            &["kind", "version", "target_temperature", "tau", "seed"],
        )?;
        self.rng = restore_rng(&state["rng"], self.rng.get_seed())?;
        Ok(())
    }
    fn apply(&mut self, system: &mut System, dt: f64) {
        let dof = system.dof() as f64;
        if dof < 1.0 || self.tau <= 0.0 || dt <= 0.0 {
            return;
        }
        let ke = system.kinetic_energy();
        if ke <= 0.0 {
            return;
        }
        let lambda = rescaling_factor(
            ke,
            0.5 * dof * self.target_temperature,
            dof,
            dt,
            self.tau,
            &mut self.rng,
        );
        for v in system.velocities.iter_mut() {
            *v *= lambda;
        }
    }
}

/// Exact Bussi–Parrinello–Nosé–Hoover velocity-rescaling factor (Bussi,
/// Donadio & Parrinello, J. Chem. Phys. 126, 014101 (2007)).
///
/// Given the current kinetic energy `K`, the target `K*` and
/// `a = exp(−dt/τ)`, the post-rescaling kinetic energy is sampled from
///
/// ```text
/// K' = a·K + (1−a)·(K*/Nf)·(R1 + R2²) + 2·R2·√(a·(1−a)·K·K*/Nf)
/// ```
///
/// with `R1 ~ χ²(Nf−1) = 2 Gamma((Nf−1)/2,1)` and `R2 ~ N(0,1)`.
/// The sign follows the Gaussian component of the transition, as required
/// for a reversible rescaling (including Nf=1). At K=0 no direction is
/// defined; callers must initialize a nonzero velocity before using CSVR.
fn rescaling_factor(
    ke: f64,
    target_ke: f64,
    ndof: f64,
    dt: f64,
    tau: f64,
    rng: &mut ChaCha12Rng,
) -> f64 {
    let a = (-dt / tau).exp();
    let r1 = if ndof > 1.0 {
        2.0 * sample_gamma(rng, (ndof - 1.0) / 2.0)
    } else {
        0.0
    };
    let r2: f64 = StandardNormal.sample(rng);
    // Sum of squares avoids cancellation at small Nf and near-zero K'.
    let noise = ((1.0 - a) * target_ke / ndof).sqrt();
    let radial = (a * ke).sqrt() + noise * r2;
    let new_ke = radial * radial + noise * noise * r1;
    (new_ke / ke).sqrt().copysign(radial)
}

/// Marsaglia–Tsang gamma sampler (shape > 0). Used for the KE distribution.
fn sample_gamma(rng: &mut ChaCha12Rng, shape: f64) -> f64 {
    if shape < 1.0 {
        let u: f64 = rng.random_range(0.0..1.0);
        return sample_gamma(rng, shape + 1.0) * u.powf(1.0 / shape);
    }
    let d = shape - 1.0 / 3.0;
    let c = (9.0 * d).sqrt().recip();
    loop {
        let x: f64 = StandardNormal.sample(rng);
        let v = 1.0 + c * x;
        if v <= 0.0 {
            continue;
        }
        let v = v * v * v;
        let u: f64 = rng.random_range(0.0..1.0);
        if u.ln() < 0.5 * x * x + d * (1.0 - v + v.ln()) {
            return d * v;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bussi_rescaling_reproduces_theoretic_moments() {
        // Canonical transition: G is chi-square(Nf), E[G]=Nf, Var[G]=2Nf.
        let k = 2.0_f64;
        let kstar = 1.0_f64;
        let ndof = 6.0_f64;
        let dt = 0.01_f64;
        let tau = 1.0_f64;
        let a = (-dt / tau).exp();
        let c1 = (1.0 - a) * kstar / ndof;
        let c2 = 2.0 * (a * (1.0 - a) * k * kstar / ndof).sqrt();
        let e_mean = a * k + c1 * ndof;
        let e_g_var = 2.0 * ndof;
        let e_mean2 = e_mean * e_mean + c1 * c1 * e_g_var + c2 * c2;

        const N: usize = 400_000;
        let mut rng = ChaCha12Rng::seed_from_u64(7);
        let mut sum = 0.0_f64;
        let mut sum2 = 0.0_f64;
        for _ in 0..N {
            let lam = rescaling_factor(k, kstar, ndof, dt, tau, &mut rng);
            let kp = k * lam * lam;
            sum += kp;
            sum2 += kp * kp;
        }
        let mean = sum / N as f64;
        let mean2 = sum2 / N as f64;
        // Monte-Carlo error on a sample of 400k with variance O(0.01) is
        // ~5e-5 for the mean and ~1e-3 for the second moment; 2e-3 bounds the
        // old implementation's ~10% mean bias by two orders of magnitude.
        assert!(
            (mean - e_mean).abs() < 2e-3,
            "mean K' = {mean}, expected {e_mean}"
        );
        assert!(
            (mean2 - e_mean2).abs() < 5e-3,
            "mean K'² = {mean2}, expected {e_mean2}"
        );
    }
}

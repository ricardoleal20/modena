//! Thermostats — velocity-only temperature control that composes with the
//! `VelocityVerlet` integrator.
//!
//! For stochastic / genuinely integrator-coupled thermostats (Langevin BAOAB),
//! see `integrator::LangevinBAOAB`, which is a self-contained integrator rather
//! than a thermostat. The split: a `Thermostat` rescales velocities given the
//! current KE; a Langevin *integrator* interleaves noise/friction with the
//! force evaluations, so it cannot be expressed as a velocity rescale.

use rand::rngs::StdRng;
use rand::{RngExt, SeedableRng};
use rand_distr::{Distribution, StandardNormal};

use crate::system::System;

/// Velocity rescaling thermostat, applied once per step.
pub trait Thermostat: Send {
    fn apply(&mut self, system: &mut System, dt: f64);
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
    rng: StdRng,
}

impl VelocityRescaling {
    pub fn new(target_temperature: f64, tau: f64, seed: u64) -> Self {
        Self {
            target_temperature,
            tau,
            rng: StdRng::seed_from_u64(seed),
        }
    }
}

impl Thermostat for VelocityRescaling {
    fn apply(&mut self, system: &mut System, dt: f64) {
        let dof = system.dof() as f64;
        if dof < 3.0 || self.tau <= 0.0 {
            return;
        }
        let ke = system.kinetic_energy();
        if ke <= 0.0 {
            return;
        }
        // Draw the rescaling factor from the exact canonical distribution of KE
        // (Bussi–Parrinello), then relax toward it with time constant tau.
        let ndof = dof;
        let target_ke = 0.5 * ndof * self.target_temperature;
        let r1 = sample_gamma(&mut self.rng, (ndof - 2.0) / 2.0);
        let r2: f64 = self.rng.random_range(0.0..1.0);
        let scale = ((r1 + r2 * r2) * target_ke / ke).sqrt();
        // Exponential relaxation toward the stochastic target scale.
        let alpha = (dt / self.tau).min(1.0);
        let lambda = (1.0 - alpha) + alpha * scale;
        for v in system.velocities.iter_mut() {
            *v *= lambda;
        }
    }
}

/// Marsaglia–Tsang gamma sampler (shape > 0). Used for the KE distribution.
fn sample_gamma(rng: &mut StdRng, shape: f64) -> f64 {
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

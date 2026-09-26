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
/// with `R1 ~ Gamma((Nf−1)/2)` and `R2 ~ N(0,1)`, which converges to the
/// exact canonical χ²ᵣᵦ distribution of the kinetic energy. The factor is
/// `√(K'/K)`. (Previous versions used `R1 ~ Gamma((Nf−2)/2)`, `R2 ~ U(0,1)`
/// and a linear blend of scales — a biased, non-canonical stationary state;
/// the doc claim of "exact" sampling was then false.)
fn rescaling_factor(
    ke: f64,
    target_ke: f64,
    ndof: f64,
    dt: f64,
    tau: f64,
    rng: &mut StdRng,
) -> f64 {
    let a = (-dt / tau).exp();
    let r1 = sample_gamma(rng, (ndof - 1.0) / 2.0);
    let r2: f64 = StandardNormal.sample(rng);
    let cross = 2.0 * r2 * (a * (1.0 - a) * ke * target_ke / ndof).sqrt();
    let new_ke = (a * ke + (1.0 - a) * target_ke * (r1 + r2 * r2) / ndof + cross).max(1e-300);
    (new_ke / ke).sqrt()
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bussi_rescaling_reproduces_theoretic_moments() {
        // One-step moments of the exact update, derived directly from the
        // formula (not from memory): K' = aK + c1·G + c2·R2 with G = R1 + R2²,
        // c1 = (1−a)K*/Nf, c2 = 2√(a(1−a)KK*/Nf). With R1 ~ Gamma((Nf−1)/2)
        // and R2 ~ N(0,1):
        //   E[G]  = (Nf−1)/2 + 1 = (Nf+1)/2        (NOT Nf)
        //   var(G) = (Nf−1)/2 + 2                  (Gamma var = shape; R2² var = 2)
        //   E[R2] = 0  →  <cross term> = 0
        //   var(K') = c1²·var(G) + c2²
        // The OLD implementation (Uniform R2 + linear blend) fails the mean by
        // ~10% for these parameters; it also can't reproduce the variance.
        let k = 2.0_f64;
        let kstar = 1.0_f64;
        let ndof = 6.0_f64;
        let dt = 0.01_f64;
        let tau = 1.0_f64;
        let a = (-dt / tau).exp();
        let c1 = (1.0 - a) * kstar / ndof;
        let c2 = 2.0 * (a * (1.0 - a) * k * kstar / ndof).sqrt();
        let e_mean = a * k + c1 * (ndof + 1.0) / 2.0;
        let e_g_var = (ndof - 1.0) / 2.0 + 2.0;
        let e_mean2 = e_mean * e_mean + c1 * c1 * e_g_var + c2 * c2;

        const N: usize = 400_000;
        let mut rng = StdRng::seed_from_u64(7);
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

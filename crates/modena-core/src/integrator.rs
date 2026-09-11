//! Integrators.
//!
//! An `Integrator` owns the time-stepping scheme and drives force evaluation
//! through a `ForceModel`. Phase 1 ships:
//! - `VelocityVerlet` — symplectic, NVE.
//! - `LangevinBAOAB`  — stochastic, NVT (Leimkuhler–Matthews BAOAB).
//!
//! Convention: `system.forces` holds F(t) on entry and F(t+dt) on exit, so the
//! runner computes forces once before the loop and the first half-kick uses the
//! cached F(t). This avoids a redundant force evaluation per step.
//!
//! M7: the kick/drift loops operate on the flat `3N` SoA buffers directly.

use rand::rngs::StdRng;
use rand::SeedableRng;
use rand_distr::{Distribution, StandardNormal};

use crate::error::Result;
use crate::force_model::{ForceModel, Virial};
use crate::neighbor::NeighborList;
use crate::system::System;
use crate::thermo::Energies;

/// Time-stepping scheme. `step` advances the system by `dt` and returns the
/// post-step thermodynamic observables.
pub trait Integrator: Send {
    fn step(
        &mut self,
        system: &mut System,
        model: &dyn ForceModel,
        neighbors: &mut NeighborList,
        dt: f64,
    ) -> Result<Energies>;
}

/// Velocity-Verlet (NVE). Expects F(t) already in `system.forces`.
pub struct VelocityVerlet;

impl VelocityVerlet {
    pub fn new() -> Self {
        Self
    }
}

impl Default for VelocityVerlet {
    fn default() -> Self {
        Self::new()
    }
}

impl Integrator for VelocityVerlet {
    fn step(
        &mut self,
        system: &mut System,
        model: &dyn ForceModel,
        neighbors: &mut NeighborList,
        dt: f64,
    ) -> Result<Energies> {
        let n = system.n();
        let half = 0.5 * dt;

        // First half-kick: v(t+dt/2) = v(t) + (dt/2) a(t), using cached F(t).
        kick(system, n, half);

        // Drift: x(t+dt) = x(t) + dt v(t+dt/2), then wrap.
        drift(system, n, dt);
        system.wrap_positions();

        // New forces F(t+dt), potential energy, virial.
        let mut virial = Virial::default();
        let pe = model.compute(system, &mut virial, neighbors)?;

        // Second half-kick: v(t+dt) = v(t+dt/2) + (dt/2) a(t+dt).
        kick(system, n, half);

        Ok(Energies::from(system, pe, &virial))
    }
}

/// BAOAB Langevin velocity-Verlet (Leimkuhler & Matthews, 2013). NVT.
///
/// Ordering per step (with F(t) cached): B½ A½ O A½ B½(recompute F).
/// The O step applies exact Ornstein–Uhlenbeck relaxation to velocities:
///   v ← c1·v + sqrt((1−c1²)·kT/m)·N(0,1),  c1 = exp(−γ·dt).
/// Seeded for reproducibility — important for research-grade runs.
pub struct LangevinBAOAB {
    pub gamma: f64,
    pub target_temperature: f64,
    rng: StdRng,
}

impl LangevinBAOAB {
    pub fn new(gamma: f64, target_temperature: f64, seed: u64) -> Self {
        Self {
            gamma,
            target_temperature,
            rng: StdRng::seed_from_u64(seed),
        }
    }
}

impl Integrator for LangevinBAOAB {
    fn step(
        &mut self,
        system: &mut System,
        model: &dyn ForceModel,
        neighbors: &mut NeighborList,
        dt: f64,
    ) -> Result<Energies> {
        let n = system.n();
        let half = 0.5 * dt;
        let c1 = (-self.gamma * dt).exp();
        let kt = self.target_temperature;
        let c1sq = c1 * c1;

        // B½: v += (dt/2) F/m   (using cached F(t))
        kick(system, n, half);
        // A½: x += (dt/2) v
        drift(system, n, half);
        // O: exact OU velocity update.
        {
            let v = &mut system.velocities;
            let m = &system.masses;
            for (i, &mi) in m.iter().enumerate() {
                let s = ((1.0 - c1sq) * kt / mi).sqrt();
                for k in 0..3 {
                    let gauss: f64 = StandardNormal.sample(&mut self.rng);
                    let idx = 3 * i + k;
                    v[idx] = c1 * v[idx] + s * gauss;
                }
            }
        }
        // A½: x += (dt/2) v
        drift(system, n, half);
        system.wrap_positions();

        // Recompute F(t+dt).
        let mut virial = Virial::default();
        let pe = model.compute(system, &mut virial, neighbors)?;

        // B½: v += (dt/2) F/m
        kick(system, n, half);

        Ok(Energies::from(system, pe, &virial))
    }
}

/// Half-kick: `v += factor * F / m` over all atoms (flat SoA).
#[inline]
fn kick(system: &mut System, n: usize, factor: f64) {
    let v = &mut system.velocities;
    let f = &system.forces;
    let m = &system.masses;
    for i in 0..n {
        let inv_m = factor / m[i];
        for k in 0..3 {
            v[3 * i + k] += f[3 * i + k] * inv_m;
        }
    }
}

/// Drift: `x += factor * v` over all atoms (flat SoA).
#[inline]
fn drift(system: &mut System, n: usize, factor: f64) {
    let x = &mut system.positions;
    let v = &system.velocities;
    for k in 0..3 * n {
        x[k] += v[k] * factor;
    }
}

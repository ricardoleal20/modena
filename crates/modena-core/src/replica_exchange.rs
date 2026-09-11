//! Replica-exchange molecular dynamics (parallel tempering).
//!
//! Runs N replicas in parallel, each at a different temperature, and
//! periodically attempts configuration swaps between adjacent replicas. A swap
//! of replicas `i` (β_i, U_i) and `j` (β_j, U_j) is accepted with probability
//!
//! ```text
//!   min(1, exp((β_j − β_i)(U_i − U_j)))
//! ```
//!
//! (potential-energy only — the velocities are rescaled to the new temperature
//! on a successful swap, so the kinetic contribution cancels). This drives
//! ergodic exploration of rough energy landscapes: a replica trapped in a local
//! minimum at low T can hop to high T, equilibrate, and return.
//!
//! Each replica runs the self-contained [`LangevinBAOAB`](crate::integrator::
//! LangevinBAOAB) NVT integrator at its own temperature, so the driver is
//! independent of the thermostat/integrator wiring in the runner.

use rand::rngs::StdRng;
use rand::RngExt;
use rand::SeedableRng;

use crate::error::Result;
use crate::force_model::{ForceModel, Virial};
use crate::integrator::{Integrator, LangevinBAOAB};
use crate::neighbor::NeighborList;
use crate::system::System;

struct Replica {
    system: System,
    neighbors: NeighborList,
    integ: LangevinBAOAB,
    temperature: f64,
    beta: f64,
    pe: f64,
}

/// Statistics from a replica-exchange run.
#[derive(Debug, Clone)]
pub struct ReplicaExchangeStats {
    pub swaps_attempted: usize,
    pub swaps_accepted: usize,
    /// Final potential energy per replica (indexed by temperature slot).
    pub final_energies: Vec<f64>,
}

/// Parallel-tempering driver.
pub struct ReplicaExchange {
    replicas: Vec<Replica>,
    rng: StdRng,
}

impl ReplicaExchange {
    /// Build `temperatures.len()` replicas, each a clone of `system` thermalized
    /// at its own temperature. `dt`/`gamma` configure the Langevin integrator.
    pub fn new(system: &System, temperatures: Vec<f64>, dt: f64, gamma: f64, seed: u64) -> Self {
        let mut rng = StdRng::seed_from_u64(seed);
        let replicas = temperatures
            .iter()
            .enumerate()
            .map(|(i, &t)| {
                let mut sys = system.clone();
                // Independent Maxwellian init at this temperature.
                init_velocities(&mut sys, t, &mut rng);
                let mut integ = LangevinBAOAB::new(gamma, t, seed.wrapping_add(i as u64));
                let _ = dt; // dt is used per-step in `run`, not at construction.
                let _ = &mut integ;
                Replica {
                    system: sys,
                    neighbors: NeighborList::new(2.5),
                    integ: LangevinBAOAB::new(gamma, t, seed.wrapping_add(i as u64)),
                    temperature: t,
                    beta: 1.0 / t,
                    pe: 0.0,
                }
            })
            .collect();
        Self { replicas, rng }
    }

    /// Run `steps` of MD with a swap attempt every `swap_every` steps.
    pub fn run(
        &mut self,
        model: &dyn ForceModel,
        steps: usize,
        swap_every: usize,
        dt: f64,
    ) -> Result<ReplicaExchangeStats> {
        // Initial PE for each replica.
        for r in &mut self.replicas {
            r.system.forces.fill(0.0);
            let mut v = Virial::default();
            r.pe = model.compute(&mut r.system, &mut v, &mut r.neighbors)?;
        }

        let mut attempted = 0usize;
        let mut accepted = 0usize;
        for step in 1..=steps {
            for r in &mut self.replicas {
                let _ = r.integ.step(&mut r.system, model, &mut r.neighbors, dt)?;
                // Track current PE (Langevin step returns Energies; recompute PE
                // cheaply via the model is avoided — the step already evaluated
                // forces; we read PE from a fresh compute only at swap time).
            }
            if swap_every > 0 && step % swap_every == 0 {
                // Refresh PE for the swap decision.
                for r in &mut self.replicas {
                    r.system.forces.fill(0.0);
                    let mut v = Virial::default();
                    r.pe = model.compute(&mut r.system, &mut v, &mut r.neighbors)?;
                }
                // Alternating even/odd swap pairs.
                let offset = (step / swap_every) % 2;
                let mut i = offset;
                while i + 1 < self.replicas.len() {
                    attempted += 1;
                    if self.try_swap(i, i + 1) {
                        accepted += 1;
                    }
                    i += 2;
                }
            }
        }

        let final_energies = self.replicas.iter().map(|r| r.pe).collect();
        Ok(ReplicaExchangeStats {
            swaps_attempted: attempted,
            swaps_accepted: accepted,
            final_energies,
        })
    }

    /// Attempt a swap between slots `i` and `j`. On accept, exchanges positions
    /// and rescales velocities to the new temperatures. Returns `true` if
    /// accepted.
    fn try_swap(&mut self, i: usize, j: usize) -> bool {
        let (pe_i, pe_j, beta_i, beta_j, t_i, t_j) = {
            let ri = &self.replicas[i];
            let rj = &self.replicas[j];
            (
                ri.pe,
                rj.pe,
                ri.beta,
                rj.beta,
                ri.temperature,
                rj.temperature,
            )
        };
        let delta = (beta_j - beta_i) * (pe_i - pe_j);
        let accept = delta <= 0.0 || {
            let r: f64 = self.rng.random_range(0.0..1.0);
            (-delta).exp() > r
        };
        if accept {
            // Swap positions; rescale velocities to the new temperature
            // (v ← v · sqrt(T_new / T_old)) so each config thermalizes correctly
            // at its (unchanged) slot temperature. `split_at_mut` gives two
            // disjoint mutable borrows (j > i).
            let scale_i = (t_j / t_i).sqrt();
            let scale_j = (t_i / t_j).sqrt();
            let (left, right) = self.replicas.split_at_mut(j);
            let ri = &mut left[i];
            let rj = &mut right[0];
            std::mem::swap(&mut ri.system.positions, &mut rj.system.positions);
            // The config that moved into slot i came from temperature t_j → rescale.
            for v in ri.system.velocities.iter_mut() {
                *v *= scale_i;
            }
            for v in rj.system.velocities.iter_mut() {
                *v *= scale_j;
            }
            std::mem::swap(&mut ri.pe, &mut rj.pe);
        }
        accept
    }

    /// Borrow the replica systems (e.g. to read final configurations).
    pub fn systems(&self) -> Vec<&System> {
        self.replicas.iter().map(|r| &r.system).collect()
    }
}

fn init_velocities(system: &mut System, temperature: f64, rng: &mut StdRng) {
    use rand_distr::{Distribution, StandardNormal};
    let n = system.n();
    let mut p_com = crate::vec3::Vec3::zero();
    let mut total_m = 0.0;
    for i in 0..n {
        let m = system.masses[i];
        let s = (temperature / m).sqrt();
        let gx: f64 = StandardNormal.sample(rng);
        let gy: f64 = StandardNormal.sample(rng);
        let gz: f64 = StandardNormal.sample(rng);
        let v = crate::vec3::Vec3::new(s * gx, s * gy, s * gz);
        system.set_vel(i, v);
        p_com += v * m;
        total_m += m;
    }
    let v_com = p_com * (1.0 / total_m);
    for i in 0..n {
        let v = system.vel(i) - v_com;
        system.set_vel(i, v);
    }
}

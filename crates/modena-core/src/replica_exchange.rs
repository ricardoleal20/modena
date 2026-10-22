//! Replica-exchange molecular dynamics (parallel tempering).
//!
//! Runs N replicas in parallel, each at a different temperature, and
//! periodically attempts configuration swaps between adjacent replicas. A swap
//! of replicas `i` (β_i, U_i) and `j` (β_j, U_j) is accepted with probability
//!
//! ```text
//!   min(1, exp((β_i − β_j)(U_i − U_j)))
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

use crate::error::{ModenaError, Result};
use crate::force_model::{ForceModel, Virial};
use crate::integrator::{Integrator, LangevinBAOAB};
use crate::neighbor::NeighborList;
use crate::system::System;

struct Replica {
    system: System,
    neighbors: Option<NeighborList>,
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
    dt: f64,
}

impl ReplicaExchange {
    /// Build `temperatures.len()` replicas, each a clone of `system` thermalized
    /// at its own temperature. `dt`/`gamma` configure the Langevin integrator.
    pub fn new(
        system: &System,
        temperatures: Vec<f64>,
        dt: f64,
        gamma: f64,
        seed: u64,
    ) -> Result<Self> {
        if temperatures.is_empty()
            || temperatures.iter().any(|t| !t.is_finite() || *t <= 0.0)
            || !dt.is_finite()
            || dt <= 0.0
            || !gamma.is_finite()
            || gamma < 0.0
        {
            return Err(ModenaError::Input("replica exchange requires positive finite temperatures/dt and nonnegative finite gamma".into()));
        }
        if system.n() == 0
            || system.masses.iter().any(|m| !m.is_finite() || *m <= 0.0)
            || system.positions.len() != 3 * system.n()
            || system.velocities.len() != 3 * system.n()
            || system.forces.len() != 3 * system.n()
            || system.positions.iter().any(|x| !x.is_finite())
        {
            return Err(ModenaError::Input(
                "invalid replica mechanical state or masses".into(),
            ));
        }
        if system.constraint_count() != 0 {
            return Err(ModenaError::Input(
                "replica exchange currently supports unconstrained systems only".into(),
            ));
        }
        let mut rng = StdRng::seed_from_u64(seed);
        let replicas = temperatures
            .iter()
            .enumerate()
            .map(|(i, &t)| {
                let mut sys = system.clone();
                // Independent Maxwellian init at this temperature.
                init_velocities(&mut sys, t, &mut rng);
                sys.set_conserved_com_dof(0)?;
                Ok(Replica {
                    system: sys,
                    neighbors: None,
                    integ: LangevinBAOAB::new(gamma, t, seed.wrapping_add(i as u64)),
                    temperature: t,
                    beta: 1.0 / t,
                    pe: 0.0,
                })
            })
            .collect::<Result<Vec<_>>>()?;
        Ok(Self { replicas, rng, dt })
    }

    /// Run `steps` of MD with a swap attempt every `swap_every` steps.
    pub fn run(
        &mut self,
        model: &dyn ForceModel,
        steps: usize,
        swap_every: usize,
        dt: f64,
    ) -> Result<ReplicaExchangeStats> {
        if !dt.is_finite() || dt <= 0.0 || dt != self.dt {
            return Err(ModenaError::Input(
                "replica exchange timestep must match its finite positive configured dt".into(),
            ));
        }
        let cutoff = model
            .neighbor_cutoff()
            .filter(|c| c.is_finite() && *c > 0.0)
            .ok_or_else(|| {
                ModenaError::Input(
                    "replica exchange requires an explicit finite positive model neighbor cutoff"
                        .into(),
                )
            })?;
        // The model defines the list; never infer a hard-coded LJ radius.
        for r in &mut self.replicas {
            r.neighbors = Some(NeighborList::new(cutoff));
            r.system.forces.fill(0.0);
            let mut v = Virial::default();
            let neighbors = r
                .neighbors
                .as_mut()
                .ok_or_else(|| ModenaError::Physics("replica neighbor list missing".into()))?;
            r.pe = model.compute(&mut r.system, &mut v, neighbors)?;
        }
        let mut attempted = 0usize;
        let mut accepted = 0usize;
        for step in 1..=steps {
            for r in &mut self.replicas {
                let neighbors = r
                    .neighbors
                    .as_mut()
                    .ok_or_else(|| ModenaError::Physics("replica neighbor list missing".into()))?;
                r.pe = r.integ.step(&mut r.system, model, neighbors, dt)?.potential;
            }
            if swap_every > 0 && step % swap_every == 0 {
                let offset = (step / swap_every - 1) % 2;
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

    /// Attempt a swap between slots `i` and `j`. On accept, exchanges the full
    /// mechanical states and their neighbor caches, then rescales velocities. Returns `true` if
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
        // log(weight_after/weight_before), with temperature slots fixed.
        let log_accept = (beta_i - beta_j) * (pe_i - pe_j);
        let accept = log_accept >= 0.0 || {
            let r: f64 = self.rng.random_range(0.0..1.0);
            r.ln() < log_accept
        };
        if accept {
            let scale_i = (t_i / t_j).sqrt();
            let scale_j = (t_j / t_i).sqrt();
            let (left, right) = self.replicas.split_at_mut(j);
            let ri = &mut left[i];
            let rj = &mut right[0];
            std::mem::swap(&mut ri.system, &mut rj.system);
            std::mem::swap(&mut ri.neighbors, &mut rj.neighbors);
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
    for i in 0..n {
        let m = system.masses[i];
        let s = (temperature / m).sqrt();
        let gx: f64 = StandardNormal.sample(rng);
        let gy: f64 = StandardNormal.sample(rng);
        let gz: f64 = StandardNormal.sample(rng);
        let v = crate::vec3::Vec3::new(s * gx, s * gy, s * gz);
        system.set_vel(i, v);
    }
}

#[cfg(test)]
mod hardening_tests {
    use super::*;
    use crate::{SimBox, Vec3};

    #[test]
    fn accepted_swap_moves_cached_forces_box_and_velocities() {
        let mut system = System::new(SimBox::cubic(20.0));
        system.resize(1);
        let mut rex = ReplicaExchange::new(&system, vec![1.0, 4.0], 0.01, 1.0, 17).unwrap();
        rex.replicas[0].system.set_pos(0, Vec3::new(2.0, 0.0, 0.0));
        rex.replicas[1].system.set_pos(0, Vec3::new(8.0, 0.0, 0.0));
        rex.replicas[0].system.set_vel(0, Vec3::new(1.0, 0.0, 0.0));
        rex.replicas[1].system.set_vel(0, Vec3::new(4.0, 0.0, 0.0));
        rex.replicas[0]
            .system
            .set_force(0, Vec3::new(-2.0, 0.0, 0.0));
        rex.replicas[1]
            .system
            .set_force(0, Vec3::new(-8.0, 0.0, 0.0));
        rex.replicas[1].system.box_size = SimBox::cubic(30.0);
        // Favorable ratio exp((beta_i-beta_j)*(Ui-Uj))=exp(15).
        rex.replicas[0].pe = 22.0;
        rex.replicas[1].pe = 2.0;
        assert!(rex.try_swap(0, 1));
        assert_eq!(
            rex.replicas[0].system.force(0),
            Vec3::new(-8.0, 0.0, 0.0),
            "force cache must move with configuration"
        );
        assert_eq!(rex.replicas[0].system.box_size, SimBox::cubic(30.0));
        assert_eq!(rex.replicas[0].system.vel(0), Vec3::new(2.0, 0.0, 0.0));
        assert_eq!(rex.replicas[1].system.vel(0), Vec3::new(2.0, 0.0, 0.0));
    }

    #[test]
    fn boltzmann_swap_rejects_strongly_unfavorable_energy_exchange() {
        let mut system = System::new(SimBox::cubic(20.0));
        system.resize(1);
        let mut rex = ReplicaExchange::new(&system, vec![1.0, 4.0], 0.01, 1.0, 17).unwrap();
        rex.replicas[0].pe = 0.0;
        rex.replicas[1].pe = 1000.0;
        assert!(!rex.try_swap(0, 1));
    }

    struct Harmonic;
    impl ForceModel for Harmonic {
        fn neighbor_cutoff(&self) -> Option<f64> {
            Some(3.7)
        }
        fn compute(
            &self,
            system: &mut System,
            _: &mut Virial,
            _: &mut NeighborList,
        ) -> Result<f64> {
            let x = system.pos(0);
            system.add_force(0, x * (-1.0));
            Ok(0.5 * x.norm2())
        }
    }

    #[test]
    fn final_energies_track_steps_after_last_swap() {
        let mut system = System::new(SimBox::cubic(20.0));
        system.resize(1);
        system.set_pos(0, Vec3::new(2.0, 1.0, 0.0));
        let mut rex = ReplicaExchange::new(&system, vec![1.0], 0.01, 1.0, 17).unwrap();
        let stats = rex.run(&Harmonic, 5, 0, 0.01).unwrap();
        let pe = 0.5 * rex.replicas[0].system.pos(0).norm2();
        assert_eq!(
            stats.final_energies[0], pe,
            "final energy must match final geometry even with no swap attempts"
        );
    }

    #[test]
    fn run_rejects_unknown_model_cutoff() {
        struct MissingCutoff;
        impl ForceModel for MissingCutoff {
            fn compute(&self, _: &mut System, _: &mut Virial, _: &mut NeighborList) -> Result<f64> {
                Ok(0.0)
            }
        }
        let mut s = System::new(SimBox::cubic(20.0));
        s.resize(1);
        let mut rex = ReplicaExchange::new(&s, vec![1.0], 0.01, 1.0, 17).unwrap();
        assert!(
            rex.run(&MissingCutoff, 1, 1, 0.01).is_err(),
            "unknown cutoff was silently guessed"
        );
    }

    #[test]
    fn run_rejects_invalid_timestep() {
        let mut s = System::new(SimBox::cubic(20.0));
        s.resize(1);
        let mut rex = ReplicaExchange::new(&s, vec![1.0], 0.01, 1.0, 17).unwrap();
        assert!(
            rex.run(&Harmonic, 1, 1, f64::NAN).is_err(),
            "nonfinite dt silently integrated"
        );
    }

    #[test]
    fn first_two_replica_swap_interval_attempts_the_pair() {
        let mut s = System::new(SimBox::cubic(20.0));
        s.resize(1);
        s.set_pos(0, Vec3::new(2.0, 1.0, 0.0));
        let mut rex = ReplicaExchange::new(&s, vec![1.0, 2.0], 0.01, 1.0, 17).unwrap();
        assert_eq!(rex.run(&Harmonic, 1, 1, 0.01).unwrap().swaps_attempted, 1);
    }

    #[test]
    fn constructor_rejects_invalid_temperatures_dt_and_gamma() {
        let mut s = System::new(SimBox::cubic(20.0));
        s.resize(1);
        for temperatures in [
            vec![],
            vec![0.0],
            vec![-1.0],
            vec![f64::NAN],
            vec![f64::INFINITY],
        ] {
            assert!(ReplicaExchange::new(&s, temperatures, 0.01, 1.0, 17).is_err());
        }
        for dt in [0.0, -0.1, f64::NAN, f64::INFINITY] {
            assert!(ReplicaExchange::new(&s, vec![1.0], dt, 1.0, 17).is_err());
        }
        for gamma in [-1.0, f64::NAN, f64::INFINITY] {
            assert!(ReplicaExchange::new(&s, vec![1.0], 0.01, gamma, 17).is_err());
        }
    }
}

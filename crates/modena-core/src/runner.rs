//! The simulation runner: builds engine objects from an `InputConfig`, runs the
//! time-stepping loop, and streams trajectory/energy output.

use std::path::Path;

use rand::rngs::StdRng;
use rand::SeedableRng;
use rand_distr::{Distribution, StandardNormal};

use crate::error::{ModenaError, Result};
use crate::force_model::{ForceModel, PairForceModel, Virial};
use crate::integrator::{Integrator, LangevinBAOAB, VelocityVerlet};
use crate::io::input::{
    InputConfig, IntegratorSpec, ModelSpec, RunSpec, StructureSpec, ThermostatSpec,
};
use crate::io::log::LogWriter;
use crate::io::xyz::{read_xyz, TrajectoryWriter};
use crate::neighbor::NeighborList;
use crate::potential::{LennardJones, PairPotential};
use crate::system::{SimBox, System};
use crate::thermo::Energies;
use crate::thermostat::{Berendsen, NoneThermostat, Thermostat, VelocityRescaling};
use crate::vec3::Vec3;

pub struct Simulation {
    system: System,
    model: Box<dyn ForceModel>,
    integrator: Box<dyn Integrator>,
    thermostat: Box<dyn Thermostat>,
    neighbors: NeighborList,
    dt: f64,
    steps: usize,
    trajectory_every: usize,
    trajectory_writer: Option<TrajectoryWriter>,
    log_every: usize,
    log_writer: Option<LogWriter>,
    current_pe: f64,
    current_virial: Virial,
}

#[derive(Debug, Clone)]
pub struct RunSummary {
    pub steps: usize,
    pub final_energies: Energies,
}

impl Simulation {
    /// Build a `Simulation` from a parsed config. `base_dir` resolves relative
    /// structure-file paths.
    pub fn from_config(config: &InputConfig, base_dir: &Path) -> Result<Self> {
        let system = build_system(&config.structure, base_dir)?;
        let mut system = system;

        init_velocities(&mut system, config.run.initial_temperature, config.run.seed);

        let (model, cutoff) = build_model(&config.model)?;
        let mut neighbors = NeighborList::new(cutoff);
        let integrator = build_integrator(&config.integrator);
        let thermostat = build_thermostat(&config.thermostat, &config.integrator);

        // Initial force evaluation so F(t) is cached for the first integrator
        // half-kick.
        let mut virial = Virial::default();
        let pe = model.compute(&mut system, &mut virial, &mut neighbors)?;

        let dt = match &config.integrator {
            IntegratorSpec::VelocityVerlet { dt } => *dt,
            IntegratorSpec::LangevinBaoab { dt, .. } => *dt,
        };

        let trajectory_writer = open_trajectory(&config.run, base_dir)?;
        let log_writer = open_log(&config.run, base_dir)?;

        Ok(Self {
            system,
            model,
            integrator,
            thermostat,
            neighbors,
            dt,
            steps: config.run.steps,
            trajectory_every: config.run.trajectory_every,
            trajectory_writer,
            log_every: config.run.log_every,
            log_writer,
            current_pe: pe,
            current_virial: virial,
        })
    }

    pub fn system(&self) -> &System {
        &self.system
    }

    pub fn run(&mut self) -> Result<RunSummary> {
        // Log initial state (step 0).
        let mut e = Energies::from(&self.system, self.current_pe, &self.current_virial);
        self.write_sample(0, 0.0, &e)?;

        let dt = self.dt;
        for step in 1..=self.steps {
            e = self
                .integrator
                .step(&mut self.system, &*self.model, &mut self.neighbors, dt)?;
            self.current_pe = e.potential;
            self.current_virial = Virial { w: e.virial };

            // Velocity-only thermostats compose with velocity-Verlet. For
            // Langevin the thermostat is `None` and the integrator already
            // handled temperature control.
            self.thermostat.apply(&mut self.system, dt);
            if !is_none_thermostat(&*self.thermostat) {
                // Velocities changed → recompute kinetic observables.
                e = Energies::from(&self.system, self.current_pe, &self.current_virial);
            }

            let time = step as f64 * dt;
            if step % self.log_every == 0 {
                self.log_sample(step, time, &e)?;
            }
            if step % self.trajectory_every == 0 {
                self.frame(step, time, &e)?;
            }
        }

        if let Some(w) = self.log_writer.as_mut() {
            w.flush()?;
        }
        Ok(RunSummary {
            steps: self.steps,
            final_energies: e,
        })
    }

    fn log_sample(&mut self, step: usize, time: f64, e: &Energies) -> Result<()> {
        if let Some(w) = self.log_writer.as_mut() {
            w.log(
                step,
                time,
                e.kinetic,
                e.potential,
                e.total,
                e.temperature,
                e.pressure,
                e.virial,
            )?;
        }
        Ok(())
    }

    fn frame(&mut self, step: usize, time: f64, e: &Energies) -> Result<()> {
        if let Some(w) = self.trajectory_writer.as_mut() {
            w.write_frame(
                &self.system,
                step,
                time,
                &format!("E={:.6} T={:.6}", e.total, e.temperature),
            )?;
        }
        Ok(())
    }

    fn write_sample(&mut self, step: usize, time: f64, e: &Energies) -> Result<()> {
        self.log_sample(step, time, e)?;
        self.frame(step, time, e)
    }
}

fn is_none_thermostat(t: &dyn Thermostat) -> bool {
    // `NoneThermostat` is the only no-op; we detect it by type name to avoid
    // adding a `is_noop()` method to the trait. Cheap, called once per step.
    let type_name = std::any::type_name_of_val(t);
    type_name.contains("NoneThermostat")
}

fn build_system(spec: &StructureSpec, base_dir: &Path) -> Result<System> {
    match spec {
        StructureSpec::Fcc {
            lattice_constant,
            replicas,
            element,
        } => Ok(crate::lattice::fcc(*lattice_constant, *replicas, element)),
        StructureSpec::Xyz {
            file,
            box_size,
            element,
        } => {
            let path = base_dir.join(file);
            let mut sys = read_xyz(&path, SimBox::from_array(*box_size))?;
            // Ensure element labels are populated.
            for el in sys.elements.iter_mut() {
                if el == "X" {
                    *el = element.clone();
                }
            }
            Ok(sys)
        }
    }
}

fn init_velocities(system: &mut System, temperature: f64, seed: u64) {
    let mut rng = StdRng::seed_from_u64(seed);
    let n = system.n();
    let normal = StandardNormal;
    let mut p_com = Vec3::zero();
    let mut total_m = 0.0;
    for i in 0..n {
        let m = system.masses[i];
        let s = (temperature / m).sqrt();
        let gx: f64 = normal.sample(&mut rng);
        let gy: f64 = normal.sample(&mut rng);
        let gz: f64 = normal.sample(&mut rng);
        let v = Vec3::new(s * gx, s * gy, s * gz);
        system.velocities[i] = v;
        p_com += v * m;
        total_m += m;
    }
    // Remove center-of-mass drift.
    let v_com = p_com * (1.0 / total_m);
    for i in 0..n {
        system.velocities[i] -= v_com;
    }
    // Rescale to the exact target temperature (compensate for COM removal).
    let ke = system.kinetic_energy();
    let dof = system.dof() as f64;
    if ke > 0.0 && dof > 0.0 {
        let current_t = 2.0 * ke / dof;
        if current_t > 0.0 {
            let lambda = (temperature / current_t).sqrt();
            for i in 0..n {
                system.velocities[i] *= lambda;
            }
        }
    }
}

fn build_model(spec: &ModelSpec) -> Result<(Box<dyn ForceModel>, f64)> {
    match spec {
        ModelSpec::LennardJones {
            sigma,
            epsilon,
            cutoff,
        } => {
            let lj = LennardJones::new(*sigma, *epsilon, *cutoff)?;
            let cutoff = lj.cutoff();
            Ok((Box::new(PairForceModel::new(lj)), cutoff))
        }
    }
}

fn build_integrator(spec: &IntegratorSpec) -> Box<dyn Integrator> {
    match spec {
        IntegratorSpec::VelocityVerlet { .. } => Box::new(VelocityVerlet::new()),
        IntegratorSpec::LangevinBaoab {
            gamma,
            target_temperature,
            seed,
            ..
        } => Box::new(LangevinBAOAB::new(*gamma, *target_temperature, *seed)),
    }
}

fn build_thermostat(
    spec: &Option<ThermostatSpec>,
    integrator: &IntegratorSpec,
) -> Box<dyn Thermostat> {
    // Langevin handles its own temperature control; ignore any thermostat block.
    if matches!(integrator, IntegratorSpec::LangevinBaoab { .. }) {
        return Box::new(NoneThermostat);
    }
    match spec {
        Some(ThermostatSpec::Berendsen {
            target_temperature,
            tau,
        }) => Box::new(Berendsen::new(*target_temperature, *tau)),
        Some(ThermostatSpec::VelocityRescaling {
            target_temperature,
            tau,
            seed,
        }) => Box::new(VelocityRescaling::new(*target_temperature, *tau, *seed)),
        _ => Box::new(NoneThermostat),
    }
}

fn open_trajectory(run: &RunSpec, base_dir: &Path) -> Result<Option<TrajectoryWriter>> {
    match &run.trajectory_file {
        Some(f) => {
            let path = base_dir.join(f);
            Ok(Some(TrajectoryWriter::create(path)?))
        }
        None => Ok(None),
    }
}

fn open_log(run: &RunSpec, base_dir: &Path) -> Result<Option<LogWriter>> {
    match &run.log_file {
        Some(f) => {
            let path = base_dir.join(f);
            Ok(Some(LogWriter::create(path)?))
        }
        None => Ok(None),
    }
}

/// Public entry point: parse + build + run from a single input file path.
pub fn run_from_file<P: AsRef<Path>>(path: P) -> Result<RunSummary> {
    let path = path.as_ref();
    let config = InputConfig::from_path(path)?;
    let base_dir = path
        .parent()
        .ok_or_else(|| ModenaError::Input("input file has no parent dir".into()))?;
    let mut sim = Simulation::from_config(&config, base_dir)?;
    sim.run()
}

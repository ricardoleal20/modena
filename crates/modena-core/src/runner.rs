//! The simulation runner: builds engine objects from an `InputConfig`, runs the
//! time-stepping loop, and streams trajectory/energy output.

use std::path::Path;
use std::sync::atomic::{AtomicBool, Ordering};

use rand::rngs::StdRng;
use rand::SeedableRng;
use rand_distr::{Distribution, StandardNormal};

use crate::conservation::ConservationReport;
use crate::error::{ModenaError, Result};
use crate::force_model::{ForceModel, PairForceModel, Virial};
use crate::integrator::{Integrator, LangevinBAOAB, VelocityVerlet};
use crate::io::input::{
    InputConfig, IntegratorSpec, InterpKindSpec, ModelSpec, RunSpec, StructureSpec, ThermostatSpec,
};
use crate::io::log::LogWriter;
use crate::io::xyz::{read_xyz, TrajectoryWriter};
use crate::lattice::{fcc, random_alloy, AlloySpecies};
use crate::neighbor::NeighborList;
use crate::potential::{
    Buckingham, BuckinghamOverride, BuckinghamTypeParams, InterpKind, LennardJones, LjOverride,
    LjTypeParams, Morse, MorseOverride, MorseTypeParams, PairPotential, SoftSphere,
    SoftSphereOverride, SoftSphereTypeParams, Tabulated,
};
use crate::system::{SimBox, System};
use crate::thermo::Energies;
use crate::thermostat::{Berendsen, NoneThermostat, Thermostat, VelocityRescaling};
use crate::units::UnitFactors;
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
    /// When `Some`, output (CSV log, trajectory, events) is converted from the
    /// internal reduced units back to real units for display. The physics stays
    /// reduced; only the rendered observables are scaled.
    unit_factors: Option<UnitFactors>,
}

/// Events emitted by [`Simulation::run_streaming`]. One consumer is the MCP
/// agent layer (`modena-mcp`); another is the future GUI. The blocking
/// [`Simulation::run`] path simply discards them.
#[derive(Debug, Clone)]
pub enum RunEvent {
    /// Emitted once before the first integration step.
    Started { atoms: usize, dt: f64, steps: usize },
    /// Thermodynamic sample at a logged step (`log_every` cadence, plus step 0).
    Thermo {
        step: usize,
        time: f64,
        energies: Energies,
    },
    /// Trajectory frame at a frame step (`trajectory_every` cadence). Positions
    /// are flattened `x,y,z` per atom as `f32` for cheap cross-boundary
    /// transfer (the agent/GUI do not need full `f64` precision for display).
    Frame {
        step: usize,
        time: f64,
        positions: Vec<f32>,
    },
    /// Run reached its step count (or was cancelled) and stopped cleanly.
    Finished { summary: RunSummary },
    /// Run aborted with an error.
    Failed { error: String },
}

#[derive(Debug, Clone, Copy)]
pub struct RunSummary {
    pub steps: usize,
    pub final_energies: Energies,
    pub conservation: ConservationReport,
}

impl Simulation {
    /// Build a `Simulation` from a parsed config. `base_dir` resolves relative
    /// structure-file paths.
    pub fn from_config(config: &InputConfig, base_dir: &Path) -> Result<Self> {
        let mut system = build_system(&config.structure, base_dir)?;

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
            unit_factors: config.unit_factors,
        })
    }

    pub fn system(&self) -> &System {
        &self.system
    }

    /// Run to completion, discarding streamed events. Equivalent to
    /// [`run_streaming`](Self::run_streaming) with no stop flag and a no-op
    /// sink. Kept as the simple blocking entry point for the CLI and tests.
    pub fn run(&mut self) -> Result<RunSummary> {
        self.run_streaming(None, |_| {})
    }

    /// Run the time-stepping loop, emitting [`RunEvent`]s to `on_event` and
    /// stopping early when `stop` (if given) is set to `true`.
    ///
    /// This is the shared enabler for the agent interface and the GUI: both
    /// consume the same step-level event stream over the same core. File output
    /// (CSV log, XYZ trajectory) is written exactly as in the blocking path, in
    /// addition to the callback — so a run observed by an agent produces the
    /// same on-disk artifacts as a CLI run.
    pub fn run_streaming(
        &mut self,
        stop: Option<&AtomicBool>,
        mut on_event: impl FnMut(&RunEvent),
    ) -> Result<RunSummary> {
        let dt = self.dt;
        let mut e = Energies::from(&self.system, self.current_pe, &self.current_virial);
        let e0 = e.total;
        let mut conservation = ConservationReport::pending(e0);

        // Step-0 sample (file + events). Internal `e` stays in reduced units
        // (conservation tracking); the displayed time/energies are scaled to
        // real units when a `unit_factors` is configured.
        self.write_sample(0, 0.0, &e)?;
        on_event(&RunEvent::Started {
            atoms: self.system.n(),
            dt: self.display_time(dt),
            steps: self.steps,
        });
        on_event(&RunEvent::Thermo {
            step: 0,
            time: 0.0,
            energies: self.display_energies(&e),
        });

        let mut last_step = 0usize;
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
            conservation.observe(e.total);

            let time_reduced = step as f64 * dt;
            let time = self.display_time(time_reduced);
            if step % self.log_every == 0 {
                self.log_sample(step, time, &e)?;
                on_event(&RunEvent::Thermo {
                    step,
                    time,
                    energies: self.display_energies(&e),
                });
            }
            if step % self.trajectory_every == 0 {
                self.frame(step, time, &e)?;
                on_event(&RunEvent::Frame {
                    step,
                    time,
                    positions: self.snapshot_positions(),
                });
            }
            last_step = step;

            if let Some(flag) = stop {
                if flag.load(Ordering::Relaxed) {
                    break;
                }
            }
        }

        if let Some(w) = self.log_writer.as_mut() {
            w.flush()?;
        }
        let summary = RunSummary {
            steps: last_step,
            final_energies: self.display_energies(&e),
            conservation,
        };
        on_event(&RunEvent::Finished { summary });
        Ok(summary)
    }

    /// Scale a reduced time to real units (ps) when factors are configured.
    fn display_time(&self, t: f64) -> f64 {
        match self.unit_factors {
            Some(f) => f.to_real_time(t),
            None => t,
        }
    }

    /// Scale reduced energies/temperature/pressure to real units when factors
    /// are configured. The internal physics is untouched.
    fn display_energies(&self, e: &Energies) -> Energies {
        match self.unit_factors {
            Some(f) => Energies {
                kinetic: f.to_real_energy(e.kinetic),
                potential: f.to_real_energy(e.potential),
                total: f.to_real_energy(e.total),
                temperature: f.to_real_temperature(e.temperature),
                pressure: f.to_real_pressure(e.pressure),
                virial: f.to_real_energy(e.virial),
            },
            None => *e,
        }
    }

    /// Flattened positions (`x,y,z` per atom) as `f32`, for `RunEvent::Frame`.
    /// Scaled to real units (Å) when factors are configured.
    fn snapshot_positions(&self) -> Vec<f32> {
        let scale = self.unit_factors.map(|f| f.sigma).unwrap_or(1.0) as f32;
        let mut s = Vec::with_capacity(self.system.n() * 3);
        for p in &self.system.positions {
            s.push(p.0[0] as f32 * scale);
            s.push(p.0[1] as f32 * scale);
            s.push(p.0[2] as f32 * scale);
        }
        s
    }

    fn log_sample(&mut self, step: usize, time: f64, e: &Energies) -> Result<()> {
        let d = self.display_energies(e);
        if let Some(w) = self.log_writer.as_mut() {
            w.log(
                step,
                time,
                d.kinetic,
                d.potential,
                d.total,
                d.temperature,
                d.pressure,
                d.virial,
            )?;
        }
        Ok(())
    }

    fn frame(&mut self, step: usize, time: f64, e: &Energies) -> Result<()> {
        let scale = self.unit_factors.map(|f| f.sigma).unwrap_or(1.0);
        let d = self.display_energies(e);
        let comment = format!("E={:.6} T={:.6}", d.total, d.temperature);
        if let Some(w) = self.trajectory_writer.as_mut() {
            w.write_frame_scaled(&self.system, step, time, &comment, scale)?;
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
            mass,
        } => {
            let mut sys = fcc(*lattice_constant, *replicas, element);
            if let Some(m) = mass {
                for ma in sys.masses.iter_mut() {
                    *ma = *m;
                }
            }
            Ok(sys)
        }
        StructureSpec::Xyz {
            file,
            box_size,
            element,
            mass,
        } => {
            let path = base_dir.join(file);
            let mut sys = read_xyz(&path, SimBox::from_array(*box_size))?;
            // Ensure element labels are populated.
            for el in sys.elements.iter_mut() {
                if el == "X" {
                    *el = element.clone();
                }
            }
            if let Some(m) = mass {
                for ma in sys.masses.iter_mut() {
                    *ma = *m;
                }
            }
            Ok(sys)
        }
        StructureSpec::RandomAlloy {
            lattice_constant,
            replicas,
            species,
            seed,
        } => {
            let alloy: Vec<AlloySpecies> = species
                .iter()
                .map(|s| AlloySpecies {
                    element: s.element.clone(),
                    mass: s.mass,
                    concentration: s.concentration,
                })
                .collect();
            random_alloy(*lattice_constant, *replicas, &alloy, *seed)
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
    let potential: Box<dyn PairPotential> = match spec {
        ModelSpec::LennardJones {
            sigma,
            epsilon,
            cutoff,
            types,
            overrides,
            ..
        } => {
            if !types.is_empty() {
                let params: Vec<LjTypeParams> = types
                    .iter()
                    .map(|t| LjTypeParams {
                        sigma: t.sigma,
                        epsilon: t.epsilon,
                    })
                    .collect();
                let ov: Vec<LjOverride> = overrides
                    .iter()
                    .map(|o| LjOverride {
                        types: o.types,
                        sigma: o.sigma,
                        epsilon: o.epsilon,
                    })
                    .collect();
                Box::new(LennardJones::from_types(
                    &params,
                    crate::potential::MixingRule::LorentzBerthelot,
                    &ov,
                    *cutoff,
                )?)
            } else {
                let s = sigma.ok_or_else(|| {
                    ModenaError::Config("LennardJones needs `sigma` or `types`".into())
                })?;
                let e = epsilon.ok_or_else(|| {
                    ModenaError::Config("LennardJones needs `epsilon` or `types`".into())
                })?;
                Box::new(LennardJones::new(s, e, *cutoff)?)
            }
        }
        ModelSpec::Buckingham {
            a,
            b,
            rho,
            cutoff,
            types,
            overrides,
        } => {
            if !types.is_empty() {
                let params: Vec<BuckinghamTypeParams> = types
                    .iter()
                    .map(|t| BuckinghamTypeParams {
                        a: t.a,
                        b: t.b,
                        rho: t.rho,
                    })
                    .collect();
                let ov: Vec<BuckinghamOverride> = overrides
                    .iter()
                    .map(|o| BuckinghamOverride {
                        types: o.types,
                        a: o.a,
                        b: o.b,
                        rho: o.rho,
                    })
                    .collect();
                Box::new(Buckingham::from_types(
                    &params,
                    crate::potential::MixingRule::LorentzBerthelot,
                    &ov,
                    *cutoff,
                )?)
            } else {
                let a =
                    a.ok_or_else(|| ModenaError::Config("Buckingham needs `a` or `types`".into()))?;
                let b =
                    b.ok_or_else(|| ModenaError::Config("Buckingham needs `b` or `types`".into()))?;
                let rho = rho.ok_or_else(|| {
                    ModenaError::Config("Buckingham needs `rho` or `types`".into())
                })?;
                Box::new(Buckingham::new(a, b, rho, *cutoff)?)
            }
        }
        ModelSpec::Morse {
            d,
            alpha,
            r0,
            cutoff,
            types,
            overrides,
        } => {
            if !types.is_empty() {
                let params: Vec<MorseTypeParams> = types
                    .iter()
                    .map(|t| MorseTypeParams {
                        d: t.d,
                        alpha: t.alpha,
                        r0: t.r0,
                    })
                    .collect();
                let ov: Vec<MorseOverride> = overrides
                    .iter()
                    .map(|o| MorseOverride {
                        types: o.types,
                        d: o.d,
                        alpha: o.alpha,
                        r0: o.r0,
                    })
                    .collect();
                Box::new(Morse::from_types(
                    &params,
                    crate::potential::MixingRule::LorentzBerthelot,
                    &ov,
                    *cutoff,
                )?)
            } else {
                let d =
                    d.ok_or_else(|| ModenaError::Config("Morse needs `d` or `types`".into()))?;
                let alpha = alpha
                    .ok_or_else(|| ModenaError::Config("Morse needs `alpha` or `types`".into()))?;
                let r0 =
                    r0.ok_or_else(|| ModenaError::Config("Morse needs `r0` or `types`".into()))?;
                Box::new(Morse::new(d, alpha, r0, *cutoff)?)
            }
        }
        ModelSpec::SoftSphere {
            sigma,
            epsilon,
            n,
            cutoff,
            types,
            overrides,
        } => {
            if !types.is_empty() {
                let params: Vec<SoftSphereTypeParams> = types
                    .iter()
                    .map(|t| SoftSphereTypeParams {
                        sigma: t.sigma,
                        epsilon: t.epsilon,
                        n: t.n,
                    })
                    .collect();
                let ov: Vec<SoftSphereOverride> = overrides
                    .iter()
                    .map(|o| SoftSphereOverride {
                        types: o.types,
                        sigma: o.sigma,
                        epsilon: o.epsilon,
                        n: o.n,
                    })
                    .collect();
                Box::new(SoftSphere::from_types(
                    &params,
                    crate::potential::MixingRule::LorentzBerthelot,
                    &ov,
                    *cutoff,
                )?)
            } else {
                let s = sigma.ok_or_else(|| {
                    ModenaError::Config("SoftSphere needs `sigma` or `types`".into())
                })?;
                let e = epsilon.ok_or_else(|| {
                    ModenaError::Config("SoftSphere needs `epsilon` or `types`".into())
                })?;
                let n = n.unwrap_or(12);
                Box::new(SoftSphere::new(s, e, n, *cutoff)?)
            }
        }
        ModelSpec::Tabulated {
            file,
            interpolation,
        } => {
            let text = std::fs::read_to_string(file).map_err(|e| {
                ModenaError::Input(format!("tabulated potential file `{file}`: {e}"))
            })?;
            let interp = match interpolation {
                InterpKindSpec::Linear => InterpKind::Linear,
            };
            Box::new(Tabulated::from_text(&text, interp)?)
        }
    };
    let cutoff = potential.cutoff();
    Ok((Box::new(PairForceModel::from_boxed(potential)), cutoff))
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

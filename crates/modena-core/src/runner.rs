//! The simulation runner: builds engine objects from an `InputConfig`, runs the
//! time-stepping loop, and streams trajectory/energy output.

use std::collections::HashSet;
use std::path::Path;
use std::sync::atomic::{AtomicBool, Ordering};

use rand::rngs::StdRng;
use rand::SeedableRng;
use rand_distr::{Distribution, StandardNormal};

use crate::barostat::{
    Barostat, BerendsenBarostat, MonteCarloBarostat, NoneBarostat, ParrinelloRahman,
};
use crate::bonded::{Angle, Bond, BondedForceModel, BondedTopology, Dihedral, DihedralTerm};
use crate::composite::CompositeForceModel;
use crate::composite::LongRangeForceModel;
use crate::conservation::ConservationReport;
use crate::constraints::{Constraint, Constraints, RigidWater};
use crate::error::{ModenaError, Result};
use crate::ewald::EwaldForceModel;
use crate::force_model::{ForceModel, PairForceModel, Virial};
use crate::integrator::{Integrator, LangevinBAOAB, VelocityVerlet};
use crate::io::input::{
    BarostatSpec, BondedTopologySpec, ConstraintsSpec, EwaldMethod, EwaldSpec, InputConfig,
    IntegratorSpec, InterpKindSpec, ModelSpec, RestartMode, RunSpec, StructureSpec, ThermostatSpec,
};
use crate::io::log::LogWriter;
use crate::io::xyz::{read_xyz_with_identity, TrajectoryWriter};
use crate::lattice::{fcc, random_alloy, rocksalt, water_lattice, AlloySpecies};
use crate::neighbor::{NeighborList, PairScale};
use crate::nose_hoover::NoseHooverChain;
use crate::potential::{
    Buckingham, BuckinghamOverride, BuckinghamTypeParams, InterpKind, LennardJones, LjOverride,
    LjTypeParams, Morse, MorseOverride, MorseTypeParams, PairPotential, SoftSphere,
    SoftSphereOverride, SoftSphereTypeParams, Tabulated,
};
use crate::spme::SpmeForceModel;
use crate::system::{SimBox, System};
use crate::thermo::Energies;
use crate::thermostat::{Berendsen, NoneThermostat, Thermostat, VelocityRescaling};
use crate::units::{UnitFactors, UnitSystem};
use crate::vec3::Vec3;

pub struct Simulation {
    system: System,
    model: Box<dyn ForceModel>,
    integrator: Box<dyn Integrator>,
    thermostat: Box<dyn Thermostat>,
    barostat: Box<dyn Barostat>,
    neighbors: NeighborList,
    dt: f64,
    steps: usize,
    /// Global step counter: 0 for a fresh run, or the restored step from a
    /// `restore_file` checkpoint. The run loop labels its steps (log, frames,
    /// checkpoints) from here so a resumed run continues logically.
    current_step: usize,
    restart_every: usize,
    restart_path: Option<std::path::PathBuf>,
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
    physical_signature: serde_json::Value,
    ensemble: Ensemble,
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

/// Thermodynamic interpretation shared by CLI and agent consumers.
#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize)]
#[serde(rename_all = "snake_case")]
pub enum Ensemble {
    Nve,
    Nvt,
    EquilibrationNpt,
    CanonicalNpt,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize)]
#[serde(rename_all = "snake_case")]
pub enum PressureStatus {
    UnconstrainedPairVirial,
    DiagnosticWithoutConstraintImpulseStress,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize)]
pub struct RunMetadata {
    pub ensemble: Ensemble,
    pub kinetic_dof: usize,
    pub constraint_count: usize,
    pub conserved_com_dof: usize,
    pub pressure_status: PressureStatus,
    /// Real thermo: kcal/mol, K, atm; times ps and coordinates Angstrom.
    pub output_units: UnitSystem,
}

#[derive(Debug, Clone, Copy)]
pub struct RunSummary {
    pub steps: usize,
    pub final_energies: Energies,
    pub conservation: ConservationReport,
    pub metadata: RunMetadata,
}

impl Simulation {
    /// Begin a fluent, code-typed configuration (the Rust "modena-language"):
    /// `Simulation::builder().fcc(..).lj(..).velocity_verlet(..).steps(..).build()?`
    /// produces the same validated `InputConfig` a TOML file produces.
    pub fn builder() -> crate::config::SimulationBuilder {
        crate::config::SimulationBuilder::new()
    }

    /// Build a `Simulation` from a parsed config. `base_dir` resolves relative
    /// structure-file paths.
    pub fn from_config(config: &InputConfig, base_dir: &Path) -> Result<Self> {
        Self::build_config(config, base_dir, true)
    }

    /// Perform all scientific setup and restart checks without opening outputs.
    pub fn from_config_without_outputs(config: &InputConfig, base_dir: &Path) -> Result<Self> {
        Self::build_config(config, base_dir, false)
    }

    fn build_config(raw: &InputConfig, base_dir: &Path, outputs: bool) -> Result<Self> {
        let normalized = raw.normalized()?;
        let config = &normalized;
        let mut system = build_system(&config.structure, base_dir, config.unit_factors.as_ref())?;

        crate::config::validation::validate_system(config, &system)?;
        let checkpoint = config
            .run
            .restore_file
            .as_ref()
            .map(|file| crate::io::restart::load(base_dir.join(file)))
            .transpose()?;
        if checkpoint.as_ref().is_some_and(|cp| cp.atoms != system.n()) {
            return Err(ModenaError::Config(
                "restart atom count differs from the structure".into(),
            ));
        }

        let (pair_model, _) = build_model(&config.model, base_dir)?;
        // Species-range validation: a multi-species structure (rocksalt sets
        // species {0,1}; random_alloy whatever it was asked for) fed to a
        // single-type model would panic with "index out of bounds" inside the
        // parallel force fold — an unhelpful crash at runtime instead of a
        // config error.
        if let Some(max_sp) = system.species.iter().max().copied() {
            let n_types = pair_model.potential.n_types();
            if max_sp as usize >= n_types {
                return Err(ModenaError::Config(format!(
                    "structure uses species index {max_sp}, but the pair model \
                     has only {n_types} type{s}; configure the model with \
                     matching [[model.types]] entries",
                    s = if n_types == 1 { "" } else { "s" }
                )));
            }
        }
        let mut scales = config
            .topology
            .as_ref()
            .map(|t| crate::config::topology_pair_scales(t, system.n()))
            .transpose()?
            .unwrap_or_default();
        // Rigid water: the structure implies its constraints (exact
        // per-molecule projection) and its intramolecular exclusions.
        let mut water_constraints: Option<Constraints> = None;
        if let StructureSpec::Water {
            n_waters,
            d_oh,
            hoh_angle,
            ..
        } = &config.structure
        {
            let pair_key = |i: usize, j: usize| (i.min(j), i.max(j));
            let mut policy: std::collections::BTreeMap<_, _> = scales
                .iter()
                .copied()
                .map(|p| (pair_key(p.i, p.j), p))
                .collect();
            let explicit: std::collections::HashMap<_, _> = config
                .topology
                .as_ref()
                .map(|t| {
                    t.exceptions
                        .iter()
                        .map(|e| (pair_key(e.atoms[0], e.atoms[1]), (e.lj, e.coulomb)))
                        .collect()
                })
                .unwrap_or_default();
            let d_hh = 2.0 * d_oh * (hoh_angle.to_radians() * 0.5).sin();
            let mut waters = Vec::with_capacity(*n_waters as usize);
            for w in 0..*n_waters as usize {
                let first = 3 * w;
                waters.push(RigidWater {
                    first,
                    d_oh: *d_oh,
                    d_hh,
                });
                for (i, j) in [
                    (first, first + 1),
                    (first, first + 2),
                    (first + 1, first + 2),
                ] {
                    // Rigid water's physical intramolecular policy is mandatory.
                    if explicit
                        .get(&pair_key(i, j))
                        .is_some_and(|&(lj, coulomb)| lj != 0.0 || coulomb != 0.0)
                    {
                        return Err(ModenaError::Config("explicit water intramolecular exceptions must preserve zero LJ/Coulomb exclusions".into()));
                    }
                    policy.insert(
                        pair_key(i, j),
                        PairScale {
                            i,
                            j,
                            lj: 0.0,
                            coulomb: 0.0,
                        },
                    );
                }
            }
            scales = policy.into_values().collect();
            water_constraints = Some(Constraints {
                bonds: Vec::new(),
                waters,
            });
        }
        let ewald = build_electrostatics(&config.electrostatics)?;
        // Wrap in a composite only when there is something to compose (bonded
        // topology or Ewald); otherwise the bare pair model stands alone.
        let model: Box<dyn ForceModel> = match (&config.topology, &ewald) {
            (None, None) => Box::new(pair_model),
            _ => {
                let mut comp = CompositeForceModel::new(pair_model);
                if let Some(topo_spec) = &config.topology {
                    comp = comp.with_bonded(BondedForceModel::new(build_topology(topo_spec)));
                }
                match ewald {
                    Some(LongRangeForceModel::Ewald(e)) => comp = comp.with_ewald(e),
                    Some(LongRangeForceModel::Spme(s)) => comp = comp.with_spme(s),
                    None => {}
                }
                Box::new(comp)
            }
        };
        let neighbor_cutoff = model.neighbor_cutoff().ok_or_else(|| {
            ModenaError::Config("force model must provide its physical neighbor cutoff".into())
        })?;
        // Verlet skin (config.run.neighbor_skin, default 0.3): the list is
        // rebuilt only when an atom moves > skin/2 since the last build. Each
        // sub-model still gates pairs by its own cutoff, so a list built at
        // neighbor_cutoff + skin is safe for every model.
        let skin = config.run.neighbor_skin;
        if !(skin >= 0.0) {
            return Err(ModenaError::Config("run.neighbor_skin must be >= 0".into()));
        }
        let mut neighbors = NeighborList::with_pair_scales(neighbor_cutoff, skin, &scales)?;
        neighbors.validate_pair_scales(system.n())?;
        let mut constraints = config.constraints.as_ref().map(build_constraints);
        match (&mut constraints, water_constraints) {
            (Some(c), Some(wc)) => c.waters.extend(wc.waters),
            (None, Some(wc)) => constraints = Some(wc),
            _ => {}
        }
        if let Some(c) = &constraints {
            validate_constraint_pairs(c)?;
        }
        let physical_signature = physical_signature(config, &system, base_dir)?;
        let count = constraints.as_ref().map_or(0, Constraints::count);
        system.set_constraint_count(count)?;
        let com_eligible = !matches!(config.integrator, IntegratorSpec::LangevinBaoab { .. })
            && model.conserves_total_momentum(&system);
        system.set_conserved_com_dof(if com_eligible && system.n() > 1 { 3 } else { 0 })?;
        let restored_step = if let Some(cp) = &checkpoint {
            system.positions = cp.positions.clone();
            system.velocities = cp.velocities.clone();
            system.box_size = cp.box_size;
            // Only a momentum-conserving model/controller may restore COM3.
            // Unsupported saved reductions remain COM0 and fail the full
            // identity check below. State-only assesses the new masses/model.
            if com_eligible {
                let com = match (&config.run.restart_mode, &cp.continuation) {
                    (RestartMode::Continue, Some(state)) => state.conserved_com_dof,
                    _ => conserved_zero_com_dof(&system),
                };
                system.set_conserved_com_dof(com)?;
            }
            if let Some(c) = &constraints {
                c.validate(&system)?;
            }
            cp.step
        } else {
            if let Some(c) = &constraints {
                c.validate(&system)?;
                c.settle(&mut system, 1e-10)?;
                let reference = system.positions.clone();
                c.shake(&mut system, &reference, 1e-10, 1000)?;
            }
            let remove_com = system.conserved_com_dof() != 0;
            draw_velocities(
                &mut system,
                config.run.initial_temperature,
                config.run.seed,
                remove_com,
            );
            crate::config::validation::validate_system(config, &system)?;
            if let Some(c) = &constraints {
                c.rattle(&mut system, 1e-10, 1000)?;
            }
            rescale_temperature(&mut system, config.run.initial_temperature);
            0
        };
        if matches!(config.thermostat,Some(ThermostatSpec::VelocityRescaling{target_temperature,..}) if target_temperature>0.0)
            && system.kinetic_energy() == 0.0
        {
            return Err(ModenaError::Config(
                "positive-temperature CSVR cannot heat zero kinetic energy".into(),
            ));
        }
        crate::config::validation::validate_system(config, &system)?;
        let mut integrator = build_integrator(&config.integrator, constraints.clone());
        let mut thermostat = build_thermostat(&config.thermostat, &config.integrator);
        let mut barostat = build_barostat(&config.barostat)?;

        // Initial force evaluation so F(t) is cached for the first integrator
        // half-kick.
        let mut virial = Virial::default();
        let mut pe = model.compute(&mut system, &mut virial, &mut neighbors)?;
        if let Some(cp) = &checkpoint {
            match config.run.restart_mode {
                RestartMode::Continue => {
                    let state=cp.continuation.as_ref().ok_or_else(||ModenaError::Config(
                        "legacy/mechanical checkpoint requires explicit restart_mode = \"state_only\"; full continuation state is absent".into()))?;
                    if state.physical_signature != physical_signature
                        || state.masses != system.masses
                        || state.species != system.species
                        || state.charges != system.charges
                        || state.elements != system.elements
                        || state.constraint_count != system.constraint_count()
                        || state.conserved_com_dof != system.conserved_com_dof()
                    {
                        return Err(ModenaError::Config("restart physical identity differs (model, dt, ensemble, seeds, topology, constraints or atoms); use restart_mode = \"state_only\" for an intentional new protocol".into()));
                    }
                    integrator.restore_state(&state.integrator)?;
                    thermostat.restore_state(&state.thermostat)?;
                    barostat.restore_state(&state.barostat)?;
                    system.forces = state.forces.clone();
                    pe = state.potential_energy;
                    virial.w = state.virial;
                }
                RestartMode::StateOnly => {
                    // State-only intentionally resets extended/RNG state; impose
                    // the new protocol's constraint manifold before force caching.
                    if let Some(c) = &constraints {
                        c.settle(&mut system, 1e-10)?;
                        let reference = system.positions.clone();
                        c.shake(&mut system, &reference, 1e-10, 1000)?;
                        c.rattle(&mut system, 1e-10, 1000)?;
                        system.zero_forces();
                        virial = Virial::default();
                        neighbors.force_rebuild(&system);
                        pe = model.compute(&mut system, &mut virial, &mut neighbors)?;
                    }
                    if com_eligible {
                        // Constraint projection can change the imported
                        // velocities; classify their resulting COM manifold.
                        system.set_conserved_com_dof(conserved_zero_com_dof(&system))?;
                    }
                }
            }
        }
        if !pe.is_finite() || !virial.w.is_finite() || system.forces.iter().any(|x| !x.is_finite())
        {
            return Err(ModenaError::Physics(
                "initial force evaluation produced nonfinite values".into(),
            ));
        }
        restored_step
            .checked_add(config.run.steps)
            .ok_or_else(|| ModenaError::Config("run step counter overflow".into()))?;

        let dt = match &config.integrator {
            IntegratorSpec::VelocityVerlet { dt } => *dt,
            IntegratorSpec::LangevinBaoab { dt, .. } => *dt,
            IntegratorSpec::NoseHooverChain { dt, .. } => *dt,
        };

        let trajectory_writer = if outputs {
            open_trajectory(&config.run, base_dir)?
        } else {
            None
        };
        let log_writer = if outputs {
            open_log(&config.run, base_dir)?
        } else {
            None
        };

        let ensemble = match config.barostat {
            Some(BarostatSpec::MonteCarlo { .. }) => Ensemble::CanonicalNpt,
            Some(_) => Ensemble::EquilibrationNpt,
            None => match config.integrator {
                IntegratorSpec::LangevinBaoab { gamma, .. } if gamma > 0.0 => Ensemble::Nvt,
                IntegratorSpec::NoseHooverChain { .. } => Ensemble::Nvt,
                _ => {
                    if config
                        .thermostat
                        .as_ref()
                        .is_some_and(|t| !matches!(t, ThermostatSpec::None))
                    {
                        Ensemble::Nvt
                    } else {
                        Ensemble::Nve
                    }
                }
            },
        };
        Ok(Self {
            ensemble,
            system,
            model,
            integrator,
            thermostat,
            barostat,
            neighbors,
            dt,
            steps: config.run.steps,
            current_step: restored_step,
            restart_every: config.run.restart_every,
            restart_path: if outputs {
                config.run.restart_file.as_ref().map(|f| base_dir.join(f))
            } else {
                None
            },
            trajectory_every: config.run.trajectory_every,
            trajectory_writer,
            log_every: config.run.log_every,
            log_writer,
            current_pe: pe,
            current_virial: virial,
            unit_factors: config.unit_factors,
            physical_signature,
        })
    }

    pub fn metadata(&self) -> RunMetadata {
        RunMetadata {
            ensemble: self.ensemble,
            kinetic_dof: self.system.dof(),
            constraint_count: self.system.constraint_count(),
            conserved_com_dof: self.system.conserved_com_dof(),
            output_units: if self.unit_factors.is_some() {
                UnitSystem::Real
            } else {
                UnitSystem::Reduced
            },
            pressure_status: if self.system.constraint_count() > 0 {
                PressureStatus::DiagnosticWithoutConstraintImpulseStress
            } else {
                PressureStatus::UnconstrainedPairVirial
            },
        }
    }

    pub fn system(&self) -> &System {
        &self.system
    }

    /// The global step counter (0 for a fresh run; the restored step from a
    /// `restore_file` checkpoint; advances as the run loop steps).
    pub fn current_step(&self) -> usize {
        self.current_step
    }

    /// Atomically write v2 mechanical, identity, DOF, cached force, RNG and
    /// extended state (reduced units) to `path`. Called by the run
    /// loop every `restart_every` steps and at run end; public for tests and
    /// callers that drive stepping manually.
    pub fn write_checkpoint(&self, path: &Path) -> Result<()> {
        crate::io::restart::save(
            path,
            &crate::io::restart::Checkpoint {
                atoms: self.system.n(),
                step: self.current_step,
                time: self.current_step as f64 * self.dt,
                box_size: self.system.box_size,
                positions: self.system.positions.clone(),
                velocities: self.system.velocities.clone(),
                continuation: Some(crate::io::restart::ContinuationState {
                    forces: self.system.forces.clone(),
                    masses: self.system.masses.clone(),
                    species: self.system.species.clone(),
                    charges: self.system.charges.clone(),
                    elements: self.system.elements.clone(),
                    constraint_count: self.system.constraint_count(),
                    conserved_com_dof: self.system.conserved_com_dof(),
                    potential_energy: self.current_pe,
                    virial: self.current_virial.w,
                    physical_signature: self.physical_signature.clone(),
                    integrator: self.integrator.checkpoint_state()?,
                    thermostat: self.thermostat.checkpoint_state()?,
                    barostat: self.barostat.checkpoint_state()?,
                }),
            },
        )
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

        // Start-of-run sample (file + events), labeled with the global step:
        // step 0 for a fresh run, the restored step for a resume (the seam
        // row repeats the checkpoint's state exactly). Internal `e` stays in
        // reduced units (conservation tracking); the displayed time/energies
        // are scaled to real units when a `unit_factors` is configured.
        let start_step = self.current_step;
        let start_time = self.display_time(start_step as f64 * dt);
        self.write_sample(start_step, start_time, &e)?;
        on_event(&RunEvent::Started {
            atoms: self.system.n(),
            dt: self.display_time(dt),
            steps: self.steps,
        });
        on_event(&RunEvent::Thermo {
            step: start_step,
            time: start_time,
            energies: self.display_energies(&e),
        });

        let mut executed = 0usize;
        for i in 1..=self.steps {
            // Global step label: a resumed run continues from the checkpoint.
            let step = start_step + i;
            self.current_step = step;
            e = self
                .integrator
                .step(&mut self.system, &*self.model, &mut self.neighbors, dt)?;
            self.current_pe = e.potential;
            self.current_virial = Virial { w: e.virial };

            // Velocity-only thermostats compose with velocity-Verlet. For
            // Langevin the thermostat is `None` and the integrator already
            // handled temperature control.
            self.thermostat.apply(&mut self.system, dt);
            // Pressure coupling: scales the box + positions. Only `NoneBarostat`
            // is a true no-op; any real barostat changes V (and thus the
            // pressure/KE observables), so recompute the thermo snapshot.
            self.barostat.apply_with_model(
                &mut self.system,
                &*self.model,
                &mut self.neighbors,
                &mut self.current_pe,
                &mut self.current_virial,
                dt,
            )?;
            e = Energies::from(&self.system, self.current_pe, &self.current_virial);
            conservation.observe(e.total);

            let time_reduced = step as f64 * dt;
            let time = self.display_time(time_reduced);
            if step.is_multiple_of(self.log_every) {
                self.log_sample(step, time, &e)?;
                on_event(&RunEvent::Thermo {
                    step,
                    time,
                    energies: self.display_energies(&e),
                });
            }
            if step.is_multiple_of(self.trajectory_every) {
                self.frame(step, time, &e)?;
                on_event(&RunEvent::Frame {
                    step,
                    time,
                    positions: self.snapshot_positions(),
                });
            }
            // Periodic checkpoint (crash recovery): every `restart_every`
            // global steps. The end-of-run write below covers the last step.
            if let Some(path) = self.restart_path.clone() {
                if step.is_multiple_of(self.restart_every.max(1)) {
                    self.write_checkpoint(&path)?;
                }
            }
            executed = i;

            if let Some(flag) = stop {
                if flag.load(Ordering::Relaxed) {
                    break;
                }
            }
        }

        if let Some(w) = self.log_writer.as_mut() {
            w.flush()?;
        }
        // Final checkpoint — also on early stop, so a run is always resumable
        // from exactly where it stopped.
        if let Some(path) = self.restart_path.clone() {
            self.write_checkpoint(&path)?;
        }
        let summary = RunSummary {
            steps: executed,
            final_energies: self.display_energies(&e),
            conservation,
            metadata: self.metadata(),
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
        let n = self.system.n();
        let mut s = Vec::with_capacity(n * 3);
        for i in 0..n {
            let p = self.system.pos(i);
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

fn build_system(
    spec: &StructureSpec,
    base_dir: &Path,
    factors: Option<&UnitFactors>,
) -> Result<System> {
    let to_reduced_length = factors.map(|f| f.to_reduced_length(1.0)).unwrap_or(1.0);
    match spec {
        StructureSpec::Fcc {
            lattice_constant,
            replicas,
            element,
            mass,
            charge,
        } => {
            let mut sys = fcc(*lattice_constant, *replicas, element);
            if let Some(m) = mass {
                for ma in sys.masses.iter_mut() {
                    *ma = *m;
                }
            }
            if let Some(q) = charge {
                for c in sys.charges.iter_mut() {
                    *c = *q;
                }
            }
            Ok(sys)
        }
        StructureSpec::Xyz {
            file,
            box_size,
            element,
            mass,
            charge,
        } => {
            let path = base_dir.join(file);
            let (mut sys, typed) = read_xyz_with_identity(&path, SimBox::from_array(*box_size))?;
            if typed {
                if mass.is_some() || charge.is_some() {
                    return Err(ModenaError::Config("typed XYZ carries per-atom mass/charge; global mass/charge overrides are only valid for plain XYZ (reference mass is 1 amu in real units)".into()));
                }
                if let Some(f) = factors {
                    for m in &mut sys.masses {
                        *m = f.to_reduced_mass(*m);
                    }
                    for q in &mut sys.charges {
                        *q = f.to_reduced_charge(*q);
                    }
                }
            }
            // Real units: `box_size`/masses/charges were converted to reduced
            // units by `convert_to_reduced`, but the *positions live in the
            // file* and were never converted — an Å-scale XYZ file under
            // units = "real" silently placed atoms ~1/σ outside the reduced
            // box, scattering them under wrapping. Scale them by the same
            // length factor. (Fcc/random-alloy/rocksalt build positions from
            // the already-converted `lattice_constant`, so they are fine.)
            if to_reduced_length != 1.0 {
                for p in sys.positions.iter_mut() {
                    *p *= to_reduced_length;
                }
            }
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
            if let Some(q) = charge {
                for c in sys.charges.iter_mut() {
                    *c = *q;
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
                    charge: s.charge,
                })
                .collect();
            random_alloy(*lattice_constant, *replicas, &alloy, *seed)
        }
        StructureSpec::RockSalt {
            lattice_constant,
            replicas,
            cation,
            anion,
            cation_mass,
            anion_mass,
            cation_charge,
            anion_charge,
        } => Ok(rocksalt(
            *lattice_constant,
            *replicas,
            cation,
            anion,
            cation_mass.unwrap_or(1.0),
            anion_mass.unwrap_or(1.0),
            *cation_charge,
            *anion_charge,
        )),
        StructureSpec::Water {
            n_waters,
            density,
            d_oh,
            hoh_angle,
            q_o,
            q_h,
            m_o,
            m_h,
            seed,
        } => water_lattice(
            *n_waters as usize,
            *density,
            *d_oh,
            *hoh_angle,
            *q_o,
            *q_h,
            *m_o,
            *m_h,
            *seed,
            to_reduced_length,
        ),
    }
}

/// Initialize velocities from a seeded Maxwell–Boltzmann distribution at
/// `temperature`, remove COM drift only when declared conserved, and rescale to the exact
/// target temperature. Public so out-of-runner engines (the GPU-d ladder
/// driver in `examples/gpu_ladder.rs`) reproduce deck initialization
/// bit-for-bit at the same seed and physical DOF. Manual callers must set
/// constraint count and conserved COM DOF before calling (normally three
/// conserved COM modes for multi-atom NVE). This helper preserves the caller's
/// metadata; Simulation sets it from the controller and force model. Charged
/// finite-mesh SPME and Langevin retain the full COM kinetic space.
pub fn init_velocities(system: &mut System, temperature: f64, seed: u64) {
    draw_velocities(system, temperature, seed, system.conserved_com_dof() != 0);
    rescale_temperature(system, temperature);
}

fn draw_velocities(system: &mut System, temperature: f64, seed: u64, remove_com: bool) {
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
        system.set_vel(i, v);
        p_com += v * m;
        total_m += m;
    }
    // Remove center-of-mass drift.
    if remove_com {
        let v_com = p_com * (1.0 / total_m);
        for i in 0..n {
            let v = system.vel(i) - v_com;
            system.set_vel(i, v);
        }
    }
}

fn rescale_temperature(system: &mut System, temperature: f64) {
    let n = system.n();
    let ke = system.kinetic_energy();
    let dof = system.dof() as f64;
    if ke > 0.0 && dof > 0.0 {
        let current_t = 2.0 * ke / dof;
        if current_t > 0.0 {
            let lambda = (temperature / current_t).sqrt();
            for i in 0..n {
                let v = system.vel(i) * lambda;
                system.set_vel(i, v);
            }
        }
    }
}

fn physical_signature(
    config: &InputConfig,
    system: &System,
    base_dir: &Path,
) -> Result<serde_json::Value> {
    let mut model =
        serde_json::to_value(&config.model).map_err(|e| ModenaError::Config(e.to_string()))?;
    if let ModelSpec::Tabulated { file, .. } = &config.model {
        model.as_object_mut().unwrap().remove("file");
        model["table_contents"] = std::fs::read_to_string(base_dir.join(file))?.into();
    }
    let factors = config.unit_factors.map(|f| {
        [
            f.sigma,
            f.epsilon,
            f.mass_ref,
            f.time,
            f.temperature,
            f.pressure,
            f.charge,
        ]
    });
    Ok(
        serde_json::json!({"schema":"modena-physical-v2", "model":model, "integrator":config.integrator,
        "thermostat":config.thermostat.as_ref().filter(|t|!matches!(t,ThermostatSpec::None)),
        "barostat":config.barostat,"electrostatics":config.electrostatics,
        "topology":config.topology,"constraints":config.constraints,
        "water_geometry":match &config.structure {StructureSpec::Water{n_waters,d_oh,hoh_angle,..}=>Some((*n_waters,*d_oh,*hoh_angle)),_=>None},
        "initialization_seed":config.run.seed,"initial_temperature":config.run.initial_temperature,
        "masses":system.masses,"species":system.species,"charges":system.charges,"elements":system.elements,
        "output_unit_factors":factors}),
    )
}

/// Reject repeated generated or explicit distance pairs before counting DOF.
/// General constraint-rank analysis remains out of scope.
fn validate_constraint_pairs(constraints: &Constraints) -> Result<()> {
    let mut pairs = HashSet::new();
    let distances =
        constraints
            .bonds
            .iter()
            .map(|b| (b.i, b.j))
            .chain(constraints.waters.iter().flat_map(|w| {
                [
                    (w.first, w.first + 1),
                    (w.first, w.first + 2),
                    (w.first + 1, w.first + 2),
                ]
            }));
    for (i, j) in distances {
        if !pairs.insert((i.min(j), i.max(j))) {
            return Err(ModenaError::Config(format!(
                "duplicate constraint pair ({i}, {j}) in explicit/generated distances; constraints must be independent"
            )));
        }
    }
    Ok(())
}

fn conserved_zero_com_dof(system: &System) -> usize {
    let momentum_scale: f64 = (0..system.n())
        .map(|i| system.mass(i) * system.vel(i).norm())
        .sum();
    if system.n() > 1 && system.total_momentum().norm() <= 64.0 * f64::EPSILON * momentum_scale {
        3
    } else {
        0
    }
}

fn build_model(spec: &ModelSpec, base_dir: &Path) -> Result<(PairForceModel, f64)> {
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
            // Resolve relative to the input file's directory like every other
            // path in the schema (XYZ structure, trajectory, log outputs).
            let text = std::fs::read_to_string(base_dir.join(file)).map_err(|e| {
                ModenaError::Input(format!(
                    "tabulated potential file `{}`: {e}",
                    base_dir.join(file).display()
                ))
            })?;
            let interp = match interpolation {
                InterpKindSpec::Linear => InterpKind::Linear,
            };
            Box::new(Tabulated::from_text(&text, interp)?)
        }
    };
    let cutoff = potential.cutoff();
    Ok((PairForceModel::from_boxed(potential), cutoff))
}

fn build_integrator(
    spec: &IntegratorSpec,
    constraints: Option<Constraints>,
) -> Box<dyn Integrator> {
    match spec {
        IntegratorSpec::VelocityVerlet { .. } => match constraints {
            Some(c) if !c.is_empty() => Box::new(VelocityVerlet::with_constraints(c)),
            _ => Box::new(VelocityVerlet::new()),
        },
        IntegratorSpec::LangevinBaoab {
            gamma,
            target_temperature,
            seed,
            ..
        } => match constraints {
            Some(c) if !c.is_empty() => Box::new(LangevinBAOAB::with_constraints(
                *gamma,
                *target_temperature,
                *seed,
                c,
            )),
            _ => Box::new(LangevinBAOAB::new(*gamma, *target_temperature, *seed)),
        },
        IntegratorSpec::NoseHooverChain {
            target_temperature,
            tau,
            chain_len,
            ..
        } => {
            // NHC handles its own temperature control (it *is* the thermostat),
            // so the separate `Thermostat` slot is forced to `NoneThermostat`
            // in `build_thermostat`. Constraints compose via SHAKE/RATTLE.
            let nhc = NoseHooverChain::new(*target_temperature, *tau, *chain_len);
            match constraints {
                Some(c) if !c.is_empty() => Box::new(nhc.with_constraints(c)),
                _ => Box::new(nhc),
            }
        }
    }
}

/// Translate the TOML bonded-topology spec into the core `BondedTopology`.
fn build_topology(spec: &BondedTopologySpec) -> BondedTopology {
    let mut topo = BondedTopology::new();
    topo.bonds = spec
        .bonds
        .iter()
        .map(|b| Bond {
            i: b.atoms[0],
            j: b.atoms[1],
            r0: b.r0,
            k: b.k,
        })
        .collect();
    topo.angles = spec
        .angles
        .iter()
        .map(|a| Angle {
            i: a.atoms[0],
            j: a.atoms[1],
            k: a.atoms[2],
            theta0: a.theta0,
            ktheta: a.k,
        })
        .collect();
    topo.dihedrals = spec
        .dihedrals
        .iter()
        .map(|d| Dihedral {
            i: d.atoms[0],
            j: d.atoms[1],
            k: d.atoms[2],
            l: d.atoms[3],
            terms: d
                .terms
                .iter()
                .map(|t| DihedralTerm {
                    n: t.n,
                    k: t.k,
                    delta: t.delta,
                })
                .collect(),
        })
        .collect();
    topo
}

/// Translate the TOML constraints spec into the core `Constraints`.
fn build_constraints(spec: &ConstraintsSpec) -> Constraints {
    Constraints {
        bonds: spec
            .bonds
            .iter()
            .map(|c| Constraint {
                i: c.atoms[0],
                j: c.atoms[1],
                d: c.length,
            })
            .collect(),
        waters: Vec::new(),
    }
}

/// Translate the TOML electrostatics spec into the selected long-range model
/// (`method = "spme"` → mesh SPME; anything else → direct Ewald). `None` → no
/// electrostatics (the model composes as a no-op).
fn build_electrostatics(spec: &Option<EwaldSpec>) -> Result<Option<LongRangeForceModel>> {
    match spec {
        Some(e) => {
            let model = match e.method {
                EwaldMethod::Ewald => LongRangeForceModel::Ewald(EwaldForceModel::new(
                    e.kappa,
                    e.real_cutoff,
                    e.kspace_cutoff,
                )?),
                EwaldMethod::Spme => LongRangeForceModel::Spme(SpmeForceModel::new(
                    e.kappa,
                    e.real_cutoff,
                    e.kspace_cutoff,
                )?),
            };
            Ok(Some(model))
        }
        None => Ok(None),
    }
}

/// Translate the TOML barostat spec into a `Barostat`. `None` → constant-volume
/// (`NoneBarostat`), applied as a no-op each step.
fn build_barostat(spec: &Option<BarostatSpec>) -> Result<Box<dyn Barostat>> {
    let result: Box<dyn Barostat> = match spec {
        Some(BarostatSpec::MonteCarlo {
            target_pressure,
            target_temperature,
            max_log_volume,
            interval,
            seed,
        }) => Box::new(
            MonteCarloBarostat::new(
                *target_pressure,
                *target_temperature,
                *max_log_volume,
                *seed,
            )
            .with_interval(*interval)?,
        ),
        Some(BarostatSpec::Berendsen {
            target_pressure,
            tau,
            compressibility,
            anisotropic,
        }) => {
            let mut b = BerendsenBarostat::new(*target_pressure, *tau, *compressibility);
            if *anisotropic {
                b = b.anisotropic();
            }
            Box::new(b)
        }
        Some(BarostatSpec::ParrinelloRahman {
            target_pressure,
            tau,
            target_temperature,
            anisotropic,
        }) => {
            let mut b = ParrinelloRahman::new(*target_pressure, *tau, *target_temperature);
            if *anisotropic {
                b = b.anisotropic();
            }
            Box::new(b)
        }
        None => Box::new(NoneBarostat),
    };
    Ok(result)
}

fn build_thermostat(
    spec: &Option<ThermostatSpec>,
    integrator: &IntegratorSpec,
) -> Box<dyn Thermostat> {
    // Langevin and the Nosé–Hoover chain handle their own temperature control;
    // ignore any thermostat block for them.
    if matches!(
        integrator,
        IntegratorSpec::LangevinBaoab { .. } | IntegratorSpec::NoseHooverChain { .. }
    ) {
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

//! The simulation runner: builds engine objects from an `InputConfig`, runs the
//! time-stepping loop, and streams trajectory/energy output.

use std::collections::HashSet;
use std::fs::{File, OpenOptions};
use std::io::{BufWriter, Write};
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
use crate::conservation::{
    ConservationApplicability, ConservationReport, EnergyErrorThresholds, EnergyMonitor,
};
use crate::constraints::{Constraint, Constraints, RigidWater};
use crate::error::{ModenaError, Result};
use crate::ewald::EwaldForceModel;
use crate::external::{AseEmtForceModel, ProviderIdentity};
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
use crate::stress::{
    endpoint_pressure, IntervalPressure, MissingTensor, MissingTensorReason, ModelEvaluation,
    PressureObservation, TensorAvailability,
};
use crate::system::{SimBox, System};
use crate::thermo::Energies;
use crate::thermostat::{Berendsen, NoneThermostat, Thermostat, VelocityRescaling};
use crate::units::{UnitFactors, UnitSystem};
use crate::vec3::Vec3;

pub struct Simulation {
    failed: bool,
    system: System,
    model: Box<dyn ForceModel>,
    integrator: Box<dyn Integrator>,
    thermostat: Option<Box<dyn Thermostat>>,
    constraints: Option<Constraints>,
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
    pressure_writer: Option<BufWriter<File>>,
    current_model: ModelEvaluation,
    // Only runner-owned integration and box coupling can change q/cell.
    geometry_generation: u64,
    model_generation: u64,
    endpoint: PressureObservation,
    interval: Option<IntervalPressure>,
    /// When `Some`, output (CSV log, trajectory, events) is converted from the
    /// internal reduced units back to real units for display. The physics stays
    /// reduced; only the rendered observables are scaled.
    unit_factors: Option<UnitFactors>,
    physical_signature: serde_json::Value,
    provider: Option<ProviderIdentity>,
    ensemble: Ensemble,
    energy_error_thresholds: Option<EnergyErrorThresholds>,
}

/// Events emitted by [`Simulation::run_streaming`]. One consumer is the MCP
/// agent layer (`modena-mcp`); another is the future GUI. The blocking
/// [`Simulation::run`] path simply discards them.
// Preserve the public by-value Finished payload; provider identity is emitted
// once per run, while the high-frequency thermo event remains allocation-free.
#[allow(clippy::large_enum_variant)]
#[derive(Debug, Clone)]
pub enum RunEvent {
    /// Emitted once before the first integration step.
    Started { atoms: usize, dt: f64, steps: usize },
    /// Thermodynamic sample at a logged step (`log_every` cadence, plus step 0).
    Thermo {
        step: usize,
        time: f64,
        energies: Energies,
        /// Full-rate statistics in reduced units, independent of log cadence.
        conservation: ConservationReport,
    },
    /// Pressure in display units at the thermo cadence. Endpoint is after controllers;
    /// interval describes the preceding mechanical integration step.
    Pressure {
        step: usize,
        time: f64,
        endpoint: PressureObservation,
        interval: Option<IntervalPressure>,
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
    Failed { step: usize, error: String },
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
    /// Configurational scalar virial derived from provider stress: -V trace(stress).
    UnconstrainedProviderStress,
    DiagnosticWithoutConstraintImpulseStress,
}

#[derive(Debug, Clone, PartialEq, Eq, serde::Serialize)]
pub struct RunMetadata {
    /// Validated installed identity for an external model; None for native models.
    pub provider: Option<ProviderIdentity>,
    pub ensemble: Ensemble,
    pub kinetic_dof: usize,
    pub constraint_count: usize,
    pub conserved_com_dof: usize,
    pub pressure_status: PressureStatus,
    /// Real thermo: kcal/mol, K, atm; times ps and coordinates Angstrom.
    pub output_units: UnitSystem,
}

#[derive(Debug, Clone)]
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
        Self::build_config(config, base_dir, true, None)
    }

    /// Perform all scientific setup and restart checks without opening outputs.
    pub fn from_config_without_outputs(config: &InputConfig, base_dir: &Path) -> Result<Self> {
        Self::build_config(config, base_dir, false, None)
    }

    /// Build with cancellation checks between setup operations. A single
    /// in-flight provider call finishes under its configured deadline.
    pub fn from_config_with_stop(
        config: &InputConfig,
        base_dir: &Path,
        stop: &AtomicBool,
    ) -> Result<Self> {
        Self::build_config(config, base_dir, true, Some(stop))
    }

    fn build_config(
        raw: &InputConfig,
        base_dir: &Path,
        outputs: bool,
        stop: Option<&AtomicBool>,
    ) -> Result<Self> {
        check_setup_stop(stop)?;
        let normalized = raw.normalized()?;
        let config = &normalized;
        let mut system = build_system(&config.structure, base_dir, config.unit_factors.as_ref())?;

        if let ModelSpec::AseEmt(spec) = &config.model {
            system.atomic_numbers = spec.atomic_numbers.clone();
        }
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

        check_setup_stop(stop)?;
        let (model, provider, mut neighbors, water_constraints): (Box<dyn ForceModel>, _, _, _) =
            if let ModelSpec::AseEmt(spec) = &config.model {
                let mut spec = spec.clone();
                spec.python = resolve_provider_python(&spec.python, base_dir)?;
                let model = AseEmtForceModel::new(
                    spec,
                    config.unit_factors.expect("validated real units"),
                )?;
                let identity = model.identity().clone();
                // Explicit provider-managed path: this valid compatibility list is
                // never built or consumed by ASE; it implies no physical cutoff.
                debug_assert!(model.uses_provider_neighbors());
                debug_assert_eq!(model.neighbor_cutoff(), None);
                (
                    Box::new(model),
                    Some(identity),
                    NeighborList::new(1.0),
                    None,
                )
            } else {
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
                            comp =
                                comp.with_bonded(BondedForceModel::new(build_topology(topo_spec)));
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
                    ModenaError::Config(
                        "force model must provide its physical neighbor cutoff".into(),
                    )
                })?;
                // Verlet skin (config.run.neighbor_skin, default 0.3): the list is
                // rebuilt only when an atom moves > skin/2 since the last build. Each
                // sub-model still gates pairs by its own cutoff, so a list built at
                // neighbor_cutoff + skin is safe for every model.
                let skin = config.run.neighbor_skin;
                if !(skin >= 0.0) {
                    return Err(ModenaError::Config("run.neighbor_skin must be >= 0".into()));
                }
                let neighbors = NeighborList::with_pair_scales(neighbor_cutoff, skin, &scales)?;
                neighbors.validate_pair_scales(system.n())?;
                (model, None, neighbors, water_constraints)
            };
        check_setup_stop(stop)?;
        let mut constraints = config.constraints.as_ref().map(build_constraints);
        match (&mut constraints, water_constraints) {
            (Some(c), Some(wc)) => c.waters.extend(wc.waters),
            (None, Some(wc)) => constraints = Some(wc),
            _ => {}
        }
        if let Some(c) = &constraints {
            validate_constraint_pairs(c)?;
        }
        let physical_signature = physical_signature(config, &system, base_dir, provider.as_ref())?;
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

        let mut restored_model = None;
        if let Some(cp) = &checkpoint {
            match config.run.restart_mode {
                RestartMode::Continue => {
                    let state=cp.continuation.as_ref().ok_or_else(||ModenaError::Config(
                        "legacy/mechanical checkpoint requires explicit restart_mode = \"state_only\"; full continuation state is absent".into()))?;
                    if !checkpoint_signature_matches(
                        config,
                        &state.physical_signature,
                        &physical_signature,
                    ) || state.masses != system.masses
                        || state.species != system.species
                        || !(state.atomic_numbers == system.atomic_numbers
                            || (state.atomic_numbers.is_empty()
                                && provider.is_none()
                                && system.atomic_numbers.iter().all(|z| *z == 0)))
                        || state.charges != system.charges
                        || state.elements != system.elements
                        || state.constraint_count != system.constraint_count()
                        || state.conserved_com_dof != system.conserved_com_dof()
                    {
                        return Err(ModenaError::Config("restart physical identity differs (model, dt, ensemble, seeds, topology, constraints or atoms); use restart_mode = \"state_only\" for an intentional new protocol".into()));
                    }
                    integrator.restore_state(&state.integrator)?;
                    match &mut thermostat {
                        Some(t) => t.restore_state(&state.thermostat)?,
                        None => NoneThermostat.restore_state(&state.thermostat)?,
                    }
                    barostat.restore_state(&state.barostat)?;
                    system.forces = state.forces.clone();
                    restored_model = Some(ModelEvaluation {
                        energy: state.potential_energy,
                        scalar_virial: state.virial,
                        tensor: missing("v2 checkpoint", MissingTensorReason::LegacyRestart),
                    });
                }
                RestartMode::StateOnly => {
                    // State-only intentionally resets extended/RNG state; impose
                    // the new protocol's constraint manifold before force caching.
                    if let Some(c) = &constraints {
                        c.settle(&mut system, 1e-10)?;
                        let reference = system.positions.clone();
                        c.shake(&mut system, &reference, 1e-10, 1000)?;
                        c.rattle(&mut system, 1e-10, 1000)?;
                    }
                    if com_eligible {
                        // Constraint projection can change the imported
                        // velocities; classify their resulting COM manifold.
                        system.set_conserved_com_dof(conserved_zero_com_dof(&system))?;
                    }
                }
            }
        }
        check_setup_stop(stop)?;
        // Evaluate once after all initialization projections. Continuation already owns F/E/W.
        let current_model = match restored_model {
            Some(saved) => saved,
            None => {
                system.zero_forces();
                model.evaluate(&mut system, &mut neighbors)?
            }
        };
        check_setup_stop(stop)?;
        if !current_model.energy.is_finite()
            || !current_model.scalar_virial.is_finite()
            || system.forces.iter().any(|x| !x.is_finite())
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

        check_setup_stop(stop)?;
        let endpoint = endpoint_pressure(&system, &current_model, constraints.as_ref())?;
        let pressure_writer = if outputs {
            open_pressure_output(&config.run, base_dir)?
        } else {
            None
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
            failed: false,
            provider,
            ensemble,
            energy_error_thresholds: config.run.energy_error_thresholds,
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
            current_model,
            geometry_generation: 0,
            model_generation: 0,
            constraints,
            endpoint,
            interval: None,
            pressure_writer,
            unit_factors: config.unit_factors,
            physical_signature,
        })
    }

    pub fn metadata(&self) -> RunMetadata {
        RunMetadata {
            provider: self.provider.clone(),
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
            } else if self.provider.is_some() {
                PressureStatus::UnconstrainedProviderStress
            } else {
                PressureStatus::UnconstrainedPairVirial
            },
        }
    }

    pub fn system(&self) -> &System {
        &self.system
    }

    /// Mechanical endpoint at the latest completed snapshot, in reduced units.
    pub fn pressure_observation(&self) -> Option<&PressureObservation> {
        (!self.failed).then_some(&self.endpoint)
    }

    /// Pre-controller interval in reduced units. Absent before the first measured step.
    pub fn interval_pressure(&self) -> Option<&IntervalPressure> {
        if self.failed {
            None
        } else {
            self.interval.as_ref()
        }
    }

    /// The global step counter (0 for a fresh run; the restored step from a
    /// `restore_file` checkpoint; advances as the run loop steps).
    pub fn current_step(&self) -> usize {
        self.current_step
    }

    /// Rejects failed instances without opening the destination; a failed step
    /// may leave partial buffers and requires construction from a good checkpoint.
    /// Atomically write v2 mechanical, identity, DOF, cached force, RNG and
    /// extended state (reduced units) to `path`. Called by the run
    /// loop every `restart_every` steps and at run end; public for tests and
    /// callers that drive stepping manually.
    pub fn write_checkpoint(&self, path: &Path) -> Result<()> {
        self.ensure_usable()?;
        self.ensure_current_model()?;
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
                    atomic_numbers: self.system.atomic_numbers.clone(),
                    charges: self.system.charges.clone(),
                    elements: self.system.elements.clone(),
                    constraint_count: self.system.constraint_count(),
                    conserved_com_dof: self.system.conserved_com_dof(),
                    potential_energy: self.current_model.energy,
                    virial: self.current_model.scalar_virial,
                    physical_signature: self.physical_signature.clone(),
                    integrator: self.integrator.checkpoint_state()?,
                    thermostat: match &self.thermostat {
                        Some(t) => t.checkpoint_state()?,
                        None => NoneThermostat.checkpoint_state()?,
                    },
                    barostat: self.barostat.checkpoint_state()?,
                }),
            },
        )
    }

    /// Run to completion, discarding streamed events. Equivalent to
    /// [`run_streaming`](Self::run_streaming) with no stop flag and a no-op
    /// sink. Kept as the simple blocking entry point for the CLI and tests.
    /// Any returned run error makes this instance terminal; construct a new
    /// simulation from a good checkpoint. Clean cancellation remains usable.
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
        self.ensure_usable()?;
        let result = self.run_streaming_inner(stop, &mut on_event);
        if result.is_err() {
            self.failed = true;
        }
        result
    }

    fn ensure_usable(&self) -> Result<()> {
        if self.failed {
            Err(ModenaError::Physics("simulation is unusable after a failed run; construct a new simulation from a good checkpoint".into()))
        } else {
            Ok(())
        }
    }

    fn run_streaming_inner(
        &mut self,
        stop: Option<&AtomicBool>,
        on_event: &mut impl FnMut(&RunEvent),
    ) -> Result<RunSummary> {
        self.ensure_current_model()?;
        let dt = self.dt;
        let mut e = self.energies();
        let e0 = e.total;
        let mut monitor = EnergyMonitor::new(
            e0,
            self.system.n(),
            e.kinetic,
            if self.ensemble == Ensemble::Nve {
                ConservationApplicability::Nve
            } else {
                ConservationApplicability::EnergyExchange
            },
        )?;
        if let Some(thresholds) = self.energy_error_thresholds {
            monitor = monitor.with_thresholds(thresholds)?;
        }

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
            conservation: monitor.report(),
        });

        self.emit_pressure(start_step, start_time, on_event)?;

        let mut executed = 0usize;
        for i in 1..=self.steps {
            if stop.is_some_and(|flag| flag.load(Ordering::Acquire)) {
                break;
            }
            // Global step label: a resumed run continues from the checkpoint.
            let step = start_step + i;
            e = match self.advance() {
                Ok(e) => e,
                Err(error) => {
                    on_event(&RunEvent::Failed {
                        step,
                        error: error.to_string(),
                    });
                    return Err(error);
                }
            };
            if let Err(error) = monitor.observe(i as f64 * dt, e.total) {
                on_event(&RunEvent::Failed {
                    step,
                    error: error.to_string(),
                });
                return Err(error);
            }

            self.current_step = step;
            let time_reduced = step as f64 * dt;
            let time = self.display_time(time_reduced);
            if step.is_multiple_of(self.log_every) {
                self.log_sample(step, time, &e)?;
                on_event(&RunEvent::Thermo {
                    step,
                    time,
                    energies: self.display_energies(&e),
                    conservation: monitor.report(),
                });
                self.emit_pressure(step, time, on_event)?;
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
        if let Some(w) = &mut self.pressure_writer {
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
            conservation: monitor.report(),
            metadata: self.metadata(),
        };
        on_event(&RunEvent::Finished {
            summary: summary.clone(),
        });
        Ok(summary)
    }

    fn energies(&self) -> Energies {
        Energies::from(
            &self.system,
            self.current_model.energy,
            &Virial {
                w: self.current_model.scalar_virial,
            },
        )
    }

    fn advance(&mut self) -> Result<Energies> {
        self.ensure_current_model()?;
        self.bump_geometry_generation()?;
        let out = self.integrator.step_with_stress(
            &mut self.system,
            &*self.model,
            &mut self.neighbors,
            self.dt,
            &self.current_model,
        )?;
        let (mut model, pre_controller_endpoint, mut interval) = match out.stress {
            Some(stress) => (
                stress.final_model,
                Some(stress.endpoint),
                Some(stress.interval),
            ),
            None => (
                legacy_missing(out.energies.potential, out.energies.virial),
                None,
                None,
            ),
        };
        if let Some(t) = &mut self.thermostat {
            t.apply(&mut self.system, self.dt);
        }
        let updated = self.barostat.evaluate_move(
            &mut self.system,
            &*self.model,
            &mut self.neighbors,
            model.energy,
            &Virial {
                w: model.scalar_virial,
            },
            self.dt,
        )?;
        let moved = updated.is_some();
        if let Some(updated) = updated {
            self.bump_geometry_generation()?;
            model = updated;
            if let Some(interval) = &mut interval {
                let mut reasons = match &interval.model_average {
                    TensorAvailability::Missing(reasons) => reasons.clone(),
                    TensorAvailability::Complete(_) => Vec::new(),
                };
                reasons.push(MissingTensor {
                    component: "barostat".into(),
                    reason: MissingTensorReason::PostStepTransform,
                });
                interval.model_average = TensorAvailability::Missing(reasons.clone());
                interval.total = TensorAvailability::Missing(reasons);
                interval.scalar = None;
            }
        }
        self.endpoint = match pre_controller_endpoint {
            Some(endpoint) if self.thermostat.is_none() && !moved => endpoint,
            None if moved => {
                let unavailable = legacy_missing(model.energy, model.scalar_virial);
                endpoint_pressure(&self.system, &unavailable, self.constraints.as_ref())?
            }
            _ => endpoint_pressure(&self.system, &model, self.constraints.as_ref())?,
        };
        self.current_model = model;
        self.model_generation = self.geometry_generation;
        self.interval = interval;
        Ok(self.energies())
    }

    fn ensure_current_model(&self) -> Result<()> {
        if self.geometry_generation != self.model_generation {
            return Err(ModenaError::Physics(
                "stale runner force/model cache".into(),
            ));
        }
        Ok(())
    }

    fn bump_geometry_generation(&mut self) -> Result<()> {
        self.geometry_generation = self
            .geometry_generation
            .checked_add(1)
            .ok_or_else(|| ModenaError::Physics("runner geometry generation overflow".into()))?;
        Ok(())
    }

    fn emit_pressure(
        &mut self,
        step: usize,
        time: f64,
        sink: &mut impl FnMut(&RunEvent),
    ) -> Result<()> {
        finite_display(time)?;
        let mut endpoint = self.endpoint.clone();
        let mut interval = self.interval.clone();
        if let Some(f) = self.unit_factors {
            scale_components(&mut endpoint.kinetic, f.epsilon)?;
            scale_components(&mut endpoint.constraint, f.epsilon)?;
            scale_tensor(&mut endpoint.model, f.epsilon)?;
            scale_tensor(&mut endpoint.total, f.pressure)?;
            endpoint.scalar = endpoint
                .scalar
                .map(|p| finite_display(f.to_real_pressure(p)))
                .transpose()?;
            endpoint.volume = finite_display(endpoint.volume * f.sigma.powi(3))?;
            if let Some(i) = &mut interval {
                scale_components(&mut i.kinetic_average, f.epsilon)?;
                scale_components(&mut i.constraint_average, f.epsilon)?;
                scale_tensor(&mut i.model_average, f.epsilon)?;
                scale_tensor(&mut i.total, f.pressure)?;
                i.scalar = i
                    .scalar
                    .map(|p| finite_display(f.to_real_pressure(p)))
                    .transpose()?;
                i.volume = finite_display(i.volume * f.sigma.powi(3))?;
                i.dt = finite_display(f.to_real_time(i.dt))?;
            }
        }
        if let Some(w) = &mut self.pressure_writer {
            #[derive(serde::Serialize)]
            struct Sample<'a> {
                step: usize,
                time: f64,
                endpoint: &'a PressureObservation,
                interval: Option<&'a IntervalPressure>,
            }
            serde_json::to_writer(
                &mut *w,
                &Sample {
                    step,
                    time,
                    endpoint: &endpoint,
                    interval: interval.as_ref(),
                },
            )
            .map_err(std::io::Error::from)?;
            w.write_all(b"\n")?;
        }
        sink(&RunEvent::Pressure {
            step,
            time,
            endpoint,
            interval,
        });
        Ok(())
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

fn legacy_missing(energy: f64, scalar_virial: f64) -> ModelEvaluation {
    ModelEvaluation {
        energy,
        scalar_virial,
        tensor: missing(
            "integrator without step stress",
            MissingTensorReason::LegacyForceModel,
        ),
    }
}

fn missing(component: &str, reason: MissingTensorReason) -> TensorAvailability {
    TensorAvailability::Missing(vec![MissingTensor {
        component: component.into(),
        reason,
    }])
}

fn finite_display(value: f64) -> Result<f64> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(ModenaError::Physics(
            "pressure display unit conversion overflow".into(),
        ))
    }
}

fn scale_components(value: &mut crate::stress::SymTensor6, factor: f64) -> Result<()> {
    for component in &mut value.0 {
        *component = finite_display(*component * factor)?;
    }
    Ok(())
}

fn scale_tensor(value: &mut TensorAvailability, factor: f64) -> Result<()> {
    if let TensorAvailability::Complete(t) = value {
        scale_components(t, factor)?;
    }
    Ok(())
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
    provider: Option<&ProviderIdentity>,
) -> Result<serde_json::Value> {
    let mut model =
        serde_json::to_value(&config.model).map_err(|e| ModenaError::Config(e.to_string()))?;
    if let ModelSpec::AseEmt(spec) = &config.model {
        model = serde_json::json!({"type":"ase_emt", "provider":provider, "atomic_numbers":system.atomic_numbers, "pbc":spec.pbc, "asap_cutoff":spec.asap_cutoff, "composition":spec.composition});
    }
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
        ModelSpec::AseEmt(_) => {
            return Err(ModenaError::Config(
                "external total provider cannot be built as pair model".into(),
            ))
        }
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
) -> Option<Box<dyn Thermostat>> {
    // Langevin and the Nosé–Hoover chain handle their own temperature control;
    // ignore any thermostat block for them.
    if matches!(
        integrator,
        IntegratorSpec::LangevinBaoab { .. } | IntegratorSpec::NoseHooverChain { .. }
    ) {
        return None;
    }
    match spec {
        Some(ThermostatSpec::Berendsen {
            target_temperature,
            tau,
        }) => Some(Box::new(Berendsen::new(*target_temperature, *tau))),
        Some(ThermostatSpec::VelocityRescaling {
            target_temperature,
            tau,
            seed,
        }) => Some(Box::new(VelocityRescaling::new(
            *target_temperature,
            *tau,
            *seed,
        ))),
        _ => None,
    }
}

/// Compare the actual files before truncating the retained pressure handle.
/// Failed checks can leave a newly created empty file at the configured path.
fn open_pressure_output(run: &RunSpec, base: &Path) -> Result<Option<BufWriter<File>>> {
    let Some(name) = &run.pressure_jsonl else {
        return Ok(None);
    };
    let outputs = [
        ("log_file", &run.log_file),
        ("trajectory_file", &run.trajectory_file),
        ("restart_file", &run.restart_file),
    ];
    // Each writer normally creates its parent. Do this before identity checks
    // so e.g. new-dir/../log has the same resolution as at write time.
    for output in std::iter::once(name).chain(outputs.iter().filter_map(|(_, name)| name.as_ref()))
    {
        if let Some(parent) = base.join(output).parent() {
            std::fs::create_dir_all(parent)?;
        }
    }
    // Pass the original spelling to the OS: "..", trailing "/" and "/." may
    // only traverse directories. Never unlink an externally replaceable name.
    let file = OpenOptions::new()
        .create(true)
        .truncate(false)
        .write(true)
        .open(base.join(name))?;
    for (field, name) in outputs
        .into_iter()
        .chain([("restore_file", &run.restore_file)])
    {
        let Some(other) = name else { continue };
        if path_matches_file(&base.join(other), &file)? {
            return Err(ModenaError::Config(format!(
                "run.pressure_jsonl and run.{field} resolve to the same file"
            )));
        }
    }
    file.set_len(0)?;
    Ok(Some(BufWriter::new(file)))
}

#[cfg(unix)]
fn path_matches_file(path: &Path, file: &File) -> std::io::Result<bool> {
    use std::os::unix::fs::MetadataExt;
    let other = match path.metadata() {
        Ok(m) => m,
        Err(e) if e.kind() == std::io::ErrorKind::NotFound => return Ok(false),
        Err(e) => return Err(e),
    };
    let opened = file.metadata()?;
    Ok(opened.dev() == other.dev() && opened.ino() == other.ino())
}

#[cfg(windows)]
fn path_matches_file(path: &Path, file: &File) -> std::io::Result<bool> {
    use std::ffi::c_void;
    use std::os::windows::{fs::OpenOptionsExt, io::AsRawHandle};

    // FILE_ID_INFO from GetFileInformationByHandleEx, available on Windows 8+.
    // The volume and 128-bit file ID identify both named aliases and hard links.
    #[repr(C)]
    #[derive(Default, PartialEq)]
    struct FileId {
        volume: u64,
        id: [u8; 16],
    }
    #[link(name = "kernel32")]
    unsafe extern "system" {
        fn GetFileInformationByHandleEx(
            file: *mut c_void,
            class: i32,
            info: *mut FileId,
            size: u32,
        ) -> i32;
    }
    fn identity(file: &File) -> std::io::Result<FileId> {
        let mut info = FileId::default();
        // SAFETY: file owns a live handle for this call. Class 18 writes
        // exactly FILE_ID_INFO to the valid, correctly sized output buffer.
        let ok = unsafe {
            GetFileInformationByHandleEx(
                file.as_raw_handle(),
                18,
                &mut info,
                std::mem::size_of::<FileId>() as u32,
            )
        };
        if ok == 0 {
            Err(std::io::Error::last_os_error())
        } else {
            Ok(info)
        }
    }
    // Only file attributes are needed, including for write-only output paths.
    let other = match OpenOptions::new().access_mode(0).open(path) {
        Ok(file) => file,
        Err(e) if e.kind() == std::io::ErrorKind::NotFound => return Ok(false),
        Err(e) => return Err(e),
    };
    Ok(identity(file)? == identity(&other)?)
}

#[cfg(not(any(unix, windows)))]
fn path_matches_file(path: &Path, _file: &File) -> std::io::Result<bool> {
    if path.try_exists()? {
        return Err(std::io::Error::new(
            std::io::ErrorKind::Unsupported,
            "this platform cannot verify output file identity",
        ));
    }
    Ok(false)
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

fn check_setup_stop(stop: Option<&AtomicBool>) -> Result<()> {
    if stop.is_some_and(|flag| flag.load(Ordering::Acquire)) {
        Err(ModenaError::Cancelled)
    } else {
        Ok(())
    }
}

/// Only the historical display-factor error may migrate, and only for real
/// fixed-cell protocols. All chemistry/provider/ensemble/RNG identity remains exact.
fn checkpoint_signature_matches(
    config: &InputConfig,
    saved: &serde_json::Value,
    current: &serde_json::Value,
) -> bool {
    if saved == current {
        return true;
    }
    if config.units != UnitSystem::Real || config.barostat.is_some() {
        return false;
    }
    let Some(f) = config.unit_factors else {
        return false;
    };
    let Some(old) = saved.get("output_unit_factors").and_then(|v| v.as_array()) else {
        return false;
    };
    let Some(new) = current
        .get("output_unit_factors")
        .and_then(|v| v.as_array())
    else {
        return false;
    };
    if old.len() != 7
        || new.len() != 7
        || old[5].as_f64() != Some(f.epsilon / (f.sigma * f.sigma * f.sigma) * 69477.8)
    {
        return false;
    }
    let mut migrated = saved.clone();
    migrated["output_unit_factors"][5] = new[5].clone();
    migrated == *current
}

/// Resolve explicit relative paths without changing PATH basenames or symlinks.
fn resolve_provider_python(python: &str, base_dir: &Path) -> Result<String> {
    let path = Path::new(python);
    if path.is_relative() && python.chars().any(std::path::is_separator) {
        let resolved = std::path::absolute(base_dir.join(path))?;
        resolved.into_os_string().into_string().map_err(|_| {
            ModenaError::Config("resolved model.python path is not valid UTF-8; select a UTF-8 executable path or a bare PATH name".into())
        })
    } else {
        Ok(python.to_owned())
    }
}

#[cfg(all(test, unix))]
mod interpreter_path_tests {
    use super::resolve_provider_python;
    use crate::ModenaError;
    use std::path::Path;

    #[test]
    fn bare_names_remain_on_path_and_explicit_paths_resolve_lexically() {
        for name in ["python3", ".emt-python"] {
            assert_eq!(
                resolve_provider_python(name, Path::new("/config")).unwrap(),
                name
            );
        }
        assert_eq!(
            resolve_provider_python("./venv/bin/python", Path::new("/config")).unwrap(),
            "/config/venv/bin/python"
        );
        assert_eq!(
            resolve_provider_python("/venv/bin/python", Path::new("/config")).unwrap(),
            "/venv/bin/python"
        );
    }

    #[cfg(unix)]
    #[test]
    fn non_utf8_base_is_explicitly_rejected_only_when_a_path_is_resolved() {
        use std::os::unix::ffi::OsStringExt;
        let base =
            std::path::PathBuf::from(std::ffi::OsString::from_vec(b"/tmp/modena-\xff".to_vec()));
        let error = resolve_provider_python("./venv/bin/python", &base).unwrap_err();
        assert!(
            matches!(error,ModenaError::Config(ref message) if message.contains("UTF-8")),
            "{error}"
        );
        assert_eq!(
            resolve_provider_python(".emt-python", &base).unwrap(),
            ".emt-python"
        );
    }
}

#[cfg(test)]
mod pressure_cache_tests {
    use super::*;
    use crate::stress::SymTensor6;
    use std::sync::{Arc, Mutex};

    const ONE_STEP: &str = r#"
[structure]
type = "fcc"
lattice_constant = 4.0
replicas = [1,1,1]
[model]
type = "lennard_jones"
sigma = 1.0
epsilon = 1.0
cutoff = 1.5
[integrator]
type = "velocity_verlet"
dt = 0.001
[run]
steps = 1
log_every = 1
"#;

    struct CountingTensor(Arc<Mutex<Vec<ModelEvaluation>>>);
    impl ForceModel for CountingTensor {
        fn compute(&self, _: &mut System, _: &mut Virial, _: &mut NeighborList) -> Result<f64> {
            panic!("runner lost the full model evaluation");
        }
        fn evaluate(&self, s: &mut System, _: &mut NeighborList) -> Result<ModelEvaluation> {
            let v = s.box_size.volume();
            let tensor = SymTensor6([v, v * 2., v * 3., v + 1., v + 2., v + 3.]);
            let m = ModelEvaluation {
                energy: v,
                scalar_virial: tensor.trace(),
                tensor: TensorAvailability::Complete(tensor),
            };
            self.0.lock().unwrap().push(m.clone());
            Ok(m)
        }
    }
    struct Scalar;
    impl ForceModel for Scalar {
        fn compute(&self, _: &mut System, w: &mut Virial, _: &mut NeighborList) -> Result<f64> {
            w.w += 3.;
            Ok(2.)
        }
    }
    struct LegacyIntegrator;
    impl Integrator for LegacyIntegrator {
        fn step(
            &mut self,
            s: &mut System,
            model: &dyn ForceModel,
            n: &mut NeighborList,
            dt: f64,
        ) -> Result<Energies> {
            VelocityVerlet::new().step(s, model, n, dt)
        }
    }

    fn install(sim: &mut Simulation, model: Box<dyn ForceModel>) {
        sim.model = model;
        sim.system.zero_forces();
        sim.current_model = sim
            .model
            .evaluate(&mut sim.system, &mut sim.neighbors)
            .unwrap();
        sim.endpoint = endpoint_pressure(&sim.system, &sim.current_model, None).unwrap();
    }

    #[test]
    fn moved_cell_refreshes_full_cache_once_and_noop_reuses_it() {
        let dir = tempfile::tempdir().unwrap();
        let config = InputConfig::from_str(ONE_STEP).unwrap();
        let mut sim = Simulation::from_config(&config, dir.path()).unwrap();
        let calls = Arc::new(Mutex::new(Vec::new()));
        install(&mut sim, Box::new(CountingTensor(calls.clone())));
        sim.run().unwrap();
        assert_eq!(calls.lock().unwrap().len(), 2);
        // Configured barostat at its exact target is also a no-op.
        let unchanged = sim.system.box_size;
        let target = sim.energies().pressure;
        sim.barostat = Box::new(BerendsenBarostat::new(target, 1., 1.));
        sim.run().unwrap();
        assert_eq!(sim.system.box_size, unchanged);
        assert_eq!(calls.lock().unwrap().len(), 3);
        sim.barostat = Box::new(BerendsenBarostat::new(-2., 1., 1.));
        sim.run().unwrap();
        let sample = sim.pressure_observation().unwrap();
        assert_ne!(sim.system.box_size, unchanged);
        let records = calls.lock().unwrap().clone();
        assert_eq!(records.len(), 5);
        assert_eq!(sample.model, records[4].tensor);
        assert_ne!(sample.model, records[3].tensor);
        let i = sim.interval_pressure().unwrap();
        assert_eq!(i.scalar, None);
        assert!(
            matches!(&i.total,TensorAvailability::Missing(m) if m.iter().any(|r|r.reason==MissingTensorReason::PostStepTransform))
        );
        sim.barostat = Box::new(NoneBarostat);
        sim.run().unwrap();
        assert_eq!(calls.lock().unwrap().len(), 6);
        assert_eq!(
            sim.interval.as_ref().unwrap().model_average,
            records[4].tensor
        );
        println!("counted calls after initialization, no barostat, exact-target no-op, moved box, next fixed step: 1, 2, 3, 5, 6; cached moved tensor {:?}",records[4].tensor);
        // Simulate an internal unpaired geometry edit, then check rejection before F/trajectory mutation.
        sim.bump_geometry_generation().unwrap();
        let q = sim.system.positions.clone();
        let mut emitted = 0;
        assert!(sim
            .run_streaming(None, |_| emitted += 1)
            .unwrap_err()
            .to_string()
            .contains("stale runner"));
        assert_eq!(emitted, 0);
        assert_eq!(q, sim.system.positions);
        assert_eq!(calls.lock().unwrap().len(), 6);
        assert!(sim.pressure_observation().is_none());
    }

    #[test]
    fn identity_mc_proposal_preserves_interval_cache_and_controller_draws() {
        use rand::RngExt;
        const WIDTH: f64 = 1e-18;
        const SEED: u64 = 310;
        let dir = tempfile::tempdir().unwrap();
        let cfg = InputConfig::from_str(ONE_STEP).unwrap();
        let mut sim = Simulation::from_config(&cfg, dir.path()).unwrap();
        let mut fixed = Simulation::from_config(&cfg, dir.path()).unwrap();
        let calls = Arc::new(Mutex::new(Vec::new()));
        install(&mut sim, Box::new(CountingTensor(calls.clone())));
        install(
            &mut fixed,
            Box::new(CountingTensor(Arc::new(Mutex::new(Vec::new())))),
        );
        let controller = || {
            MonteCarloBarostat::new(1., 1., WIDTH, SEED)
                .with_interval(2)
                .unwrap()
        };
        sim.barostat = Box::new(controller());
        let mut expected_rng = chacha20::ChaCha12Rng::seed_from_u64(SEED);
        let mut accepted = 0;
        for step in 1..=6 {
            sim.run().unwrap();
            fixed.run().unwrap();
            if step % 2 == 0 {
                let delta = expected_rng.random_range(-WIDTH..WIDTH);
                assert_eq!((delta / 3.).exp(), 1.);
                let log_draw = expected_rng.random::<f64>().ln();
                let log_ratio = (sim.system.n() + 1) as f64 * delta;
                accepted += usize::from(log_ratio >= 0. || log_draw < log_ratio);
            }
            assert_eq!(calls.lock().unwrap().len(), step + 1);
            assert_eq!(sim.interval_pressure(), fixed.interval_pressure());
            assert!(sim.interval_pressure().unwrap().scalar.is_some());
            assert_eq!(sim.system.positions, fixed.system.positions);
            assert_eq!(sim.system.forces, fixed.system.forces);
            let state = sim.barostat.checkpoint_state().unwrap();
            assert_eq!(state["attempted"], step / 2);
            assert_eq!(state["elapsed_steps"], step);
            assert_eq!(state["accepted"], accepted);
            assert_eq!(
                state["rng"],
                crate::thermostat::checkpoint_rng(&expected_rng).unwrap()
            );
            if step == 3 {
                let mut resumed = controller();
                resumed.restore_state(&state).unwrap();
                sim.barostat = Box::new(resumed);
            }
        }
        assert_eq!(accepted, 3);
        println!("identity MC steps=6, calls={}, attempts=3 accepted={accepted}, RNG final word {}, restored after step 3", calls.lock().unwrap().len(), expected_rng.get_word_pos());
    }

    #[test]
    fn nonfinite_display_tensor_cannot_be_written_as_null() {
        let dir = tempfile::tempdir().unwrap();
        let mut config = InputConfig::from_str(ONE_STEP).unwrap();
        config.run.pressure_jsonl = Some("pressure.jsonl".into());
        let mut sim = Simulation::from_config(&config, dir.path()).unwrap();
        install(
            &mut sim,
            Box::new(CountingTensor(Arc::new(Mutex::new(Vec::new())))),
        );
        let mut factors = UnitFactors::new(1., 1., 1.);
        factors.pressure = f64::MAX;
        sim.unit_factors = Some(factors);
        assert!(sim
            .run()
            .unwrap_err()
            .to_string()
            .contains("display unit conversion overflow"));
        assert!(sim.pressure_observation().is_none());
        assert!(std::fs::read(dir.path().join("pressure.jsonl"))
            .unwrap()
            .is_empty());
    }

    #[test]
    fn legacy_provider_and_custom_integrator_report_explicit_missing_diagnostics() {
        let dir = tempfile::tempdir().unwrap();
        for custom in [false, true] {
            let path = format!("pressure-{custom}.jsonl");
            let mut config = InputConfig::from_str(ONE_STEP).unwrap();
            config.run.pressure_jsonl = Some(path.clone());
            let mut sim = Simulation::from_config(&config, dir.path()).unwrap();
            if custom {
                sim.integrator = Box::new(LegacyIntegrator);
            } else {
                install(&mut sim, Box::new(Scalar));
            }
            let e = sim.run().unwrap().final_energies;
            let raw = std::fs::read_to_string(dir.path().join(path)).unwrap();
            let row: serde_json::Value = serde_json::from_str(raw.lines().last().unwrap()).unwrap();
            let endpoint = &row["endpoint"];
            assert!(endpoint["scalar"].is_null());
            assert_eq!(endpoint["total"]["kind"], "missing");
            assert_eq!(
                endpoint["model"]["value"][0]["reason"]["kind"],
                "legacy_force_model"
            );
            assert_eq!(endpoint["validation"]["kind"], "diagnostic");
            if custom {
                assert!(row["interval"].is_null());
                assert_eq!(
                    endpoint["model"]["value"][0]["component"],
                    "integrator without step stress"
                );
            } else {
                assert_eq!(e.virial, 3.);
                assert!(row["interval"]["scalar"].is_null());
                assert_eq!(row["interval"]["total"]["kind"], "missing");
            }
            println!(
                "custom={custom}, legacy W={}, endpoint={}, interval={}",
                e.virial, endpoint["model"], row["interval"]["total"]["kind"]
            );
        }
    }
}

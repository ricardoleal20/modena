//! The Rust "modena-language": a typed, fluent builder that produces the same
//! validated `InputConfig` a TOML file produces (docs/PLAN_MODENA_PY.md,
//! M-RL-0). One validator (`InputConfig::validate`) serves TOML, Python, and
//! Rust.
//!
//! ```rust,ignore
//! let sim = Simulation::builder()
//!     .fcc(1.6796, [4, 4, 4], "Ar")
//!     .lj(1.0, 1.0, 2.5)
//!     .velocity_verlet(0.005)
//!     .steps(2000).temperature(0.8).seed(42)
//!     .build()?;
//! ```

use crate::error::{ModenaError, Result};
use crate::io::input::{
    BarostatSpec, BondedTopologySpec, ConstraintsSpec, EwaldSpec, InputConfig, IntegratorSpec,
    ModelSpec, RunSpec, StructureSpec, ThermostatSpec,
};
use crate::units::UnitSystem;

/// Fluent builder over [`InputConfig`]. Omitted optional sections stay `None`;
/// omitted run fields take the same defaults as the TOML parser.
#[derive(Default)]
pub struct SimulationBuilder {
    structure: Option<StructureSpec>,
    model: Option<ModelSpec>,
    integrator: Option<IntegratorSpec>,
    thermostat: Option<ThermostatSpec>,
    topology: Option<BondedTopologySpec>,
    constraints: Option<ConstraintsSpec>,
    barostat: Option<BarostatSpec>,
    electrostatics: Option<EwaldSpec>,
    run: RunSpec,
    unit_factors: Option<crate::units::UnitFactors>,
}

impl SimulationBuilder {
    pub fn new() -> Self {
        Self {
            run: RunSpec {
                steps: 100,
                initial_temperature: 1.0,
                seed: 42,
                neighbor_skin: 0.3,
                trajectory_file: None,
                trajectory_every: 100,
                log_file: None,
                log_every: 10,
                restore_file: None,
                restart_file: None,
                restart_every: 500,
            },
            ..Default::default()
        }
    }

    // ---- structure ----
    pub fn structure(mut self, s: StructureSpec) -> Self {
        self.structure = Some(s);
        self
    }
    /// FCC lattice with `element`, `replicas` unit cells and `lattice_constant`.
    pub fn fcc(self, lattice_constant: f64, replicas: [u32; 3], element: &str) -> Self {
        self.structure(StructureSpec::Fcc {
            lattice_constant,
            replicas,
            element: element.to_string(),
            mass: None,
            charge: None,
        })
    }

    // ---- model ----
    pub fn model(mut self, m: ModelSpec) -> Self {
        self.model = Some(m);
        self
    }
    /// Single-species shifted-force Lennard-Jones.
    pub fn lj(self, sigma: f64, epsilon: f64, cutoff: f64) -> Self {
        self.model(ModelSpec::LennardJones {
            sigma: Some(sigma),
            epsilon: Some(epsilon),
            cutoff,
            types: Vec::new(),
            overrides: Vec::new(),
            mixing: crate::io::input::MixingRuleSpec::LorentzBerthelot,
        })
    }

    // ---- integrator ----
    pub fn integrator(mut self, i: IntegratorSpec) -> Self {
        self.integrator = Some(i);
        self
    }
    pub fn velocity_verlet(self, dt: f64) -> Self {
        self.integrator(IntegratorSpec::VelocityVerlet { dt })
    }
    pub fn langevin(self, dt: f64, gamma: f64, target_temperature: f64, seed: u64) -> Self {
        self.integrator(IntegratorSpec::LangevinBaoab {
            dt,
            gamma,
            target_temperature,
            seed,
        })
    }

    // ---- optional sections ----
    pub fn thermostat(mut self, t: ThermostatSpec) -> Self {
        self.thermostat = Some(t);
        self
    }
    pub fn topology(mut self, t: BondedTopologySpec) -> Self {
        self.topology = Some(t);
        self
    }
    pub fn constraints(mut self, c: ConstraintsSpec) -> Self {
        self.constraints = Some(c);
        self
    }
    pub fn barostat(mut self, b: BarostatSpec) -> Self {
        self.barostat = Some(b);
        self
    }
    pub fn electrostatics(mut self, e: EwaldSpec) -> Self {
        self.electrostatics = Some(e);
        self
    }
    // ---- run ----
    pub fn run(mut self, r: RunSpec) -> Self {
        self.run = r;
        self
    }
    pub fn steps(mut self, s: usize) -> Self {
        self.run.steps = s;
        self
    }
    pub fn temperature(mut self, t: f64) -> Self {
        self.run.initial_temperature = t;
        self
    }
    pub fn seed(mut self, s: u64) -> Self {
        self.run.seed = s;
        self
    }
    pub fn skin(mut self, skin: f64) -> Self {
        self.run.neighbor_skin = skin;
        self
    }
    pub fn trajectory(mut self, file: impl Into<String>, every: usize) -> Self {
        self.run.trajectory_file = Some(file.into());
        self.run.trajectory_every = every;
        self
    }
    pub fn log(mut self, file: impl Into<String>, every: usize) -> Self {
        self.run.log_file = Some(file.into());
        self.run.log_every = every;
        self
    }
    /// Write a checkpoint to `file` every `every` steps (and at run end).
    pub fn restart(mut self, file: impl Into<String>, every: usize) -> Self {
        self.run.restart_file = Some(file.into());
        self.run.restart_every = every;
        self
    }
    /// Restore the state (positions/velocities/box/step) from `file`.
    pub fn restore(mut self, file: impl Into<String>) -> Self {
        self.run.restore_file = Some(file.into());
        self
    }
    /// Render the validated [`InputConfig`] (allows custom `base_dir`).
    pub fn to_config(self) -> Result<InputConfig> {
        let structure = self
            .structure
            .ok_or_else(|| ModenaError::Config("builder: missing structure".into()))?;
        let model = self
            .model
            .ok_or_else(|| ModenaError::Config("builder: missing model".into()))?;
        let integrator = self
            .integrator
            .ok_or_else(|| ModenaError::Config("builder: missing integrator".into()))?;
        let config = InputConfig {
            units: UnitSystem::Reduced,
            structure,
            model,
            integrator,
            thermostat: self.thermostat,
            topology: self.topology,
            constraints: self.constraints,
            barostat: self.barostat,
            electrostatics: self.electrostatics,
            run: self.run,
            unit_factors: self.unit_factors,
        };
        config.validate()?;
        Ok(config)
    }
    /// Validate + build a [`Simulation`] (base dir `"."`).
    pub fn build(self) -> Result<crate::runner::Simulation> {
        self.build_with_base(std::path::Path::new("."))
    }
    /// Validate + build a [`Simulation`] with `base_dir` for relative paths.
    pub fn build_with_base(self, base: &std::path::Path) -> Result<crate::runner::Simulation> {
        let config = self.to_config()?;
        crate::runner::Simulation::from_config(&config, base)
    }
}

//! Input configuration (TOML). Describes one simulation declaratively; the
//! runner builds the engine objects from it.
//!
//! Example:
//! ```toml
//! units = "reduced"
//!
//! [Structure]
//! type = "fcc"
//! lattice_constant = 1.5444
//! replicas = [4, 4, 4]
//! element = "Ar"
//!
//! [Model]
//! type = "lennard_jones"
//! sigma = 1.0
//! epsilon = 1.0
//! cutoff = 2.5
//!
//! [Integrator]
//! type = "velocity_verlet"
//! dt = 0.005
//!
//! [Run]
//! steps = 10000
//! initial_temperature = 1.0
//! seed = 42
//! trajectory_file = "traj.xyz"
//! trajectory_every = 100
//! log_file = "energy.csv"
//! log_every = 10
//! ```

use std::path::Path;

use serde::Deserialize;

use crate::error::{ModenaError, Result};
use crate::units::UnitSystem;

#[derive(Debug, Clone, Deserialize)]
pub struct InputConfig {
    #[serde(default)]
    pub units: UnitSystem,
    pub structure: StructureSpec,
    pub model: ModelSpec,
    pub integrator: IntegratorSpec,
    #[serde(default)]
    pub thermostat: Option<ThermostatSpec>,
    pub run: RunSpec,
}

#[derive(Debug, Clone, Deserialize)]
#[serde(tag = "type", rename_all = "snake_case")]
pub enum StructureSpec {
    /// FCC lattice generated in-place.
    Fcc {
        lattice_constant: f64,
        replicas: [u32; 3],
        #[serde(default = "default_element")]
        element: String,
    },
    /// Read positions from an XYZ file.
    Xyz {
        file: String,
        box_size: [f64; 3],
        #[serde(default = "default_element")]
        element: String,
    },
}

fn default_element() -> String {
    "X".to_string()
}

#[derive(Debug, Clone, Deserialize)]
#[serde(tag = "type", rename_all = "snake_case")]
pub enum ModelSpec {
    LennardJones {
        sigma: f64,
        epsilon: f64,
        cutoff: f64,
    },
}

#[derive(Debug, Clone, Deserialize)]
#[serde(tag = "type", rename_all = "snake_case")]
pub enum IntegratorSpec {
    VelocityVerlet {
        dt: f64,
    },
    LangevinBaoab {
        dt: f64,
        gamma: f64,
        target_temperature: f64,
        #[serde(default = "default_seed")]
        seed: u64,
    },
}

fn default_seed() -> u64 {
    42
}

#[derive(Debug, Clone, Deserialize)]
#[serde(tag = "type", rename_all = "snake_case")]
pub enum ThermostatSpec {
    Berendsen {
        target_temperature: f64,
        #[serde(default = "default_tau")]
        tau: f64,
    },
    VelocityRescaling {
        target_temperature: f64,
        #[serde(default = "default_tau")]
        tau: f64,
        #[serde(default = "default_seed")]
        seed: u64,
    },
    None,
}

fn default_tau() -> f64 {
    0.5
}

#[derive(Debug, Clone, Deserialize)]
pub struct RunSpec {
    pub steps: usize,
    #[serde(default = "default_initial_temp")]
    pub initial_temperature: f64,
    #[serde(default = "default_seed")]
    pub seed: u64,
    #[serde(default)]
    pub trajectory_file: Option<String>,
    #[serde(default = "default_traj_every")]
    pub trajectory_every: usize,
    #[serde(default)]
    pub log_file: Option<String>,
    #[serde(default = "default_log_every")]
    pub log_every: usize,
}

fn default_initial_temp() -> f64 {
    1.0
}
fn default_traj_every() -> usize {
    100
}
fn default_log_every() -> usize {
    10
}

impl InputConfig {
    pub fn from_path<P: AsRef<Path>>(path: P) -> Result<Self> {
        let text = std::fs::read_to_string(path)?;
        let config: InputConfig = toml::from_str(&text)?;
        if config.units == UnitSystem::Real {
            return Err(ModenaError::Config(
                "Real units are not wired in Phase 1; use `units = \"reduced\"`. See ROADMAP M2."
                    .into(),
            ));
        }
        Ok(config)
    }
}

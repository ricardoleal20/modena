//! MODENA — Molecular Dynamics Engine for Numerical Atomistics.
//!
//! A fast, extensible molecular-dynamics framework written in Rust. Phase 1
//! (this release) provides:
//!
//! - Lennard-Jones pair potential (energy-shifted, cutoff).
//! - Cell-list neighbor search for orthorhombic periodic boxes.
//! - Velocity-Verlet (NVE) and BAOAB-Langevin (NVT) integrators.
//! - Berendsen and stochastic velocity-rescaling thermostats.
//! - Rayon-parallel force evaluation.
//! - XYZ trajectory output + CSV energy logging.
//! - TOML-driven input config.
//! - A `ForceModel` trait seam for future ML surrogate potentials.
//!
//! Units are reduced (LJ) in Phase 1. See `ROADMAP.md` for the staged plan
//! (bonds/angles, Coulomb+Ewald, NPT, real units, Python bindings, ML, GPU).
//!
//! **Naming note:** the working name is MODENA. The alternate brand **RHEA**
//! (Rust High-performance Engine for Atomistics) is held in reserve; the crate
//! is structured so a rename is mechanical (crate name + this doc header).

#![allow(clippy::needless_doctest_main)]

pub mod error;
pub mod force_model;
pub mod integrator;
pub mod io;
pub mod lattice;
pub mod neighbor;
pub mod potential;
pub mod runner;
pub mod system;
pub mod thermo;
pub mod thermostat;
pub mod units;
pub mod vec3;

pub use error::{ModenaError, Result};
pub use force_model::{ForceModel, PairForceModel, Virial};
pub use integrator::{Integrator, LangevinBAOAB, VelocityVerlet};
pub use io::input::{
    InputConfig, IntegratorSpec, ModelSpec, RunSpec, StructureSpec, ThermostatSpec,
};
pub use io::log::LogWriter;
pub use io::xyz::{read_xyz, TrajectoryWriter};
pub use neighbor::NeighborList;
pub use potential::{LennardJones, PairPotential};
pub use runner::{run_from_file, RunSummary, Simulation};
pub use system::{SimBox, System};
pub use thermo::Energies;
pub use thermostat::{Berendsen, NoneThermostat, Thermostat, VelocityRescaling};
pub use units::UnitSystem;
pub use vec3::Vec3;

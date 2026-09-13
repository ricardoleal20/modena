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
// The potential validators use `!(x > 0.0)` rather than `x <= 0.0` deliberately:
// the negated form rejects NaN inputs (NaN > 0.0 is false → !(false) = true →
// error), which the direct comparison would silently accept. Clippy's
// `neg_cmp_op_on_partial_ord` lint dislikes the style; the NaN-safety is load-
// bearing, so we allow it crate-wide.
#![allow(clippy::neg_cmp_op_on_partial_ord)]

pub mod barostat;
pub mod bonded;
pub mod composite;
pub mod conservation;
pub mod constraints;
pub mod decomp;
pub mod error;
pub mod ewald;
pub mod force_model;
#[cfg(feature = "gpu")]
pub mod gpu_comm;
pub mod integrator;
pub mod io;
pub mod lattice;
#[cfg(feature = "mpi")]
pub mod mpi_comm;
pub mod neighbor;
pub mod nose_hoover;
pub mod potential;
pub mod replica_exchange;
pub mod runner;
pub mod system;
pub mod thermo;
pub mod thermostat;
pub mod units;
pub mod vec3;

pub use barostat::{
    Barostat, BerendsenBarostat, MonteCarloBarostat, NoneBarostat, ParrinelloRahman,
};
pub use bonded::{Angle, Bond, BondedForceModel, BondedTopology, Dihedral, DihedralTerm};
pub use composite::CompositeForceModel;
pub use conservation::{ConservationReport, ConservationVerdict};
pub use constraints::{Constraint, Constraints};
pub use decomp::{
    Communicator, DomainDecomposition, DomainForceModel, LocalCommunicator, RankGrid,
};
pub use error::{ModenaError, Result};
pub use ewald::EwaldForceModel;
pub use force_model::{ForceModel, PairForceModel, Virial};
#[cfg(feature = "gpu")]
pub use gpu_comm::GpuPairForceModel;
pub use integrator::{Integrator, LangevinBAOAB, VelocityVerlet};
pub use io::input::{
    BondedTopologySpec, ConstraintSpec, ConstraintsSpec, EwaldSpec, InputConfig, IntegratorSpec,
    ModelSpec, RunSpec, StructureSpec, ThermostatSpec,
};
pub use io::log::LogWriter;
pub use io::xyz::{read_xyz, TrajectoryWriter};
pub use lattice::{random_alloy, rocksalt, AlloySpecies};
pub use neighbor::NeighborList;
pub use nose_hoover::NoseHooverChain;
pub use potential::mixing::{MixingRule, PairParamTable};
pub use potential::{
    Buckingham, BuckinghamOverride, BuckinghamTypeParams, InterpKind, LennardJones, LjOverride,
    LjTypeParams, Morse, MorseOverride, MorseTypeParams, PairPotential, SoftSphere,
    SoftSphereOverride, SoftSphereTypeParams, Tabulated,
};
pub use replica_exchange::{ReplicaExchange, ReplicaExchangeStats};
pub use runner::{run_from_file, RunEvent, RunSummary, Simulation};
pub use system::{SimBox, System};
pub use thermo::Energies;
pub use thermostat::{Berendsen, NoneThermostat, Thermostat, VelocityRescaling};
pub use units::{UnitFactors, UnitSystem};
pub use vec3::Vec3;

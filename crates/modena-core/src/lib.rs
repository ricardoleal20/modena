//! MODENA — Molecular Dynamics Engine for Numerical Atomistics.
//!
//! An early development classical molecular dynamics engine in Rust. The library
//! includes pair and bonded force models, constraints, electrostatics, integrators,
//! ensemble controls, and a TOML-driven runner with XYZ and CSV output.
//!
//! Numerical regression coverage applies to selected models and workflows.
//! Constrained pressure remains diagnostic; constrained NPT and full anisotropic
//! coupling are unsupported. See `docs/VALIDATED_SCOPE.md` in the repository for
//! method coverage, units, continuation compatibility, and backend restrictions.

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
pub mod config;
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
pub mod spme;
pub mod system;
pub mod thermo;
pub mod thermostat;
pub mod units;
pub mod vec3;

pub use barostat::{
    Barostat, BerendsenBarostat, MonteCarloBarostat, NoneBarostat, ParrinelloRahman,
};
pub use bonded::{Angle, Bond, BondedForceModel, BondedTopology, Dihedral, DihedralTerm};
pub use composite::{CompositeForceModel, LongRangeForceModel};
pub use conservation::{ConservationReport, ConservationVerdict};
pub use constraints::{Constraint, Constraints};
pub use decomp::{
    Communicator, DomainDecomposition, DomainForceModel, LocalCommunicator, RankGrid,
};
pub use error::{ModenaError, Result};
pub use ewald::EwaldForceModel;
pub use force_model::{ForceModel, PairForceModel, Virial};
#[cfg(feature = "gpu")]
pub use gpu_comm::{GpuGather, GpuPairForceModel};
pub use integrator::{Integrator, LangevinBAOAB, VelocityVerlet};
pub use io::input::{
    BondedTopologySpec, ConstraintSpec, ConstraintsSpec, EwaldMethod, EwaldSpec, InputConfig,
    IntegratorSpec, ModelSpec, RestartMode, RunSpec, StructureSpec, ThermostatSpec,
};
pub use io::log::LogWriter;
pub use io::xyz::{read_xyz, TrajectoryWriter};
pub use lattice::{random_alloy, rocksalt, AlloySpecies};
pub use neighbor::{NeighborList, PairScale};
pub use nose_hoover::NoseHooverChain;
pub use potential::mixing::{MixingRule, PairParamTable};
pub use potential::{
    Buckingham, BuckinghamOverride, BuckinghamTypeParams, InterpKind, LennardJones, LjOverride,
    LjTypeParams, Morse, MorseOverride, MorseTypeParams, PairPotential, SoftSphere,
    SoftSphereOverride, SoftSphereTypeParams, Tabulated,
};
pub use replica_exchange::{ReplicaExchange, ReplicaExchangeStats};
pub use runner::{
    run_from_file, Ensemble, PressureStatus, RunEvent, RunMetadata, RunSummary, Simulation,
};
pub use spme::SpmeForceModel;
pub use system::{SimBox, System};
pub use thermo::Energies;
pub use thermostat::{Berendsen, NoneThermostat, Thermostat, VelocityRescaling};
pub use units::{UnitFactors, UnitSystem};
pub use vec3::Vec3;

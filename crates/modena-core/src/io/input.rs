//! Input configuration (TOML). Describes one simulation declaratively; the
//! runner builds the engine objects from it.
//!
//! Example (reduced, single-species LJ):
//! ```toml
//! units = "reduced"
//!
//! [structure]
//! type = "fcc"
//! lattice_constant = 1.5444
//! replicas = [4, 4, 4]
//! element = "Ar"
//!
//! [model]
//! type = "lennard_jones"
//! sigma = 1.0
//! epsilon = 1.0
//! cutoff = 2.5
//!
//! [integrator]
//! type = "velocity_verlet"
//! dt = 0.005
//!
//! [run]
//! steps = 10000
//! initial_temperature = 1.0
//! seed = 42
//! ```
//!
//! Real units + multi-species (M2):
//! ```toml
//! units = "real"
//!
//! [structure]
//! type = "random_alloy"
//! lattice_constant = 5.26
//! replicas = [4, 4, 4]
//! seed = 42
//! [[structure.species]]
//! element = "Ar"
//! mass = 39.948
//! concentration = 0.8
//! [[structure.species]]
//! element = "Ne"
//! mass = 20.180
//! concentration = 0.2
//!
//! [model]
//! type = "lennard_jones"
//! cutoff = 12.0
//! [[model.types]]
//! sigma = 3.40
//! epsilon = 0.238
//! [[model.types]]
//! sigma = 2.75
//! epsilon = 0.072
//! ```

use std::path::Path;

use serde::Deserialize;

use crate::error::{ModenaError, Result};
use crate::units::{UnitFactors, UnitSystem};

#[derive(Debug, Clone, Deserialize)]
pub struct InputConfig {
    #[serde(default)]
    pub units: UnitSystem,
    pub structure: StructureSpec,
    pub model: ModelSpec,
    pub integrator: IntegratorSpec,
    #[serde(default)]
    pub thermostat: Option<ThermostatSpec>,
    /// Bonded topology (bonds/angles/dihedrals). Optional — absent for a pure
    /// pair (non-bonded) system. M3.
    #[serde(default)]
    pub topology: Option<BondedTopologySpec>,
    /// Rigid distance constraints (SHAKE/RATTLE). Optional. M3.
    #[serde(default)]
    pub constraints: Option<ConstraintsSpec>,
    /// Pressure coupling (NPT). Optional. M5.
    #[serde(default)]
    pub barostat: Option<BarostatSpec>,
    /// Long-range Coulomb electrostatics (Ewald). Optional. M4.
    #[serde(default)]
    pub electrostatics: Option<EwaldSpec>,
    pub run: RunSpec,
    /// Conversion factors used to render output back to real units. Populated by
    /// [`from_str`](Self::from_str) when `units = "real"`; `None` in reduced
    /// mode (output is already in reduced units). Not deserialized from TOML.
    #[serde(skip)]
    pub unit_factors: Option<UnitFactors>,
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
        /// Per-element mass (real units: amu; reduced: 1.0). Optional — defaults
        /// to 1.0 (reduced) which is correct for single-species LJ.
        #[serde(default)]
        mass: Option<f64>,
        /// Charge on every atom (M4 electrostatics; reduced units). Defaults to 0.
        #[serde(default)]
        charge: Option<f64>,
    },
    /// Read positions from an XYZ file.
    Xyz {
        file: String,
        box_size: [f64; 3],
        #[serde(default = "default_element")]
        element: String,
        #[serde(default)]
        mass: Option<f64>,
        /// Charge on every atom (M4 electrostatics; reduced units). Defaults to 0.
        #[serde(default)]
        charge: Option<f64>,
    },
    /// FCC lattice with per-atom species assigned by concentration (random alloy).
    RandomAlloy {
        lattice_constant: f64,
        replicas: [u32; 3],
        /// Per-species definition: element, mass, charge, and concentration
        /// (concentrations must sum to 1).
        species: Vec<SpeciesSpec>,
        #[serde(default = "default_seed")]
        seed: u64,
    },
    /// Rocksalt (NaCl) lattice: ordered cation/anion sublattices with charges.
    /// M4 — the canonical ionic test structure. Nearest cation–anion distance
    /// is `lattice_constant / 2`.
    RockSalt {
        lattice_constant: f64,
        replicas: [u32; 3],
        cation: String,
        anion: String,
        #[serde(default)]
        cation_mass: Option<f64>,
        #[serde(default)]
        anion_mass: Option<f64>,
        cation_charge: f64,
        anion_charge: f64,
    },
}

#[derive(Debug, Clone, Deserialize)]
pub struct SpeciesSpec {
    pub element: String,
    /// Mass (real units: amu; reduced: relative to the reference mass).
    #[serde(default = "default_mass")]
    pub mass: f64,
    /// Fraction of atoms assigned this species (0..1; must sum to 1 across species).
    pub concentration: f64,
    /// Per-atom charge for this species (M4 electrostatics; reduced units).
    /// Defaults to 0 (uncharged).
    #[serde(default)]
    pub charge: f64,
}

// ---------------------------------------------------------------------------
// Bonded topology + constraints (M3).
// ---------------------------------------------------------------------------

/// One harmonic bond: `atoms = [i, j]`, equilibrium length `r0`, stiffness `k`.
/// `U = ½ k (r − r0)²`.
#[derive(Debug, Clone, Copy, Deserialize)]
pub struct BondSpec {
    pub atoms: [usize; 2],
    pub r0: f64,
    pub k: f64,
}

/// One harmonic angle: `atoms = [i, j, k]` (`j` is the vertex), equilibrium
/// angle `theta0` (radians), stiffness `k`. `U = ½ k (θ − θ0)²`.
#[derive(Debug, Clone, Copy, Deserialize)]
pub struct AngleSpec {
    pub atoms: [usize; 3],
    pub theta0: f64,
    pub k: f64,
}

/// One periodic dihedral term (OPLS/CHARMM form `k(1 + cos(nφ − δ))`).
#[derive(Debug, Clone, Copy, Deserialize)]
pub struct DihedralTermSpec {
    pub n: i32,
    pub k: f64,
    pub delta: f64,
}

/// One dihedral: `atoms = [i, j, k, l]` plus its Fourier terms.
#[derive(Debug, Clone, Deserialize)]
pub struct DihedralSpec {
    pub atoms: [usize; 4],
    pub terms: Vec<DihedralTermSpec>,
}

/// The full bonded topology.
#[derive(Debug, Clone, Default, Deserialize)]
pub struct BondedTopologySpec {
    #[serde(default)]
    pub bonds: Vec<BondSpec>,
    #[serde(default)]
    pub angles: Vec<AngleSpec>,
    /// Exclude bonded 1-2/1-3/1-4 pairs from the *pair* model (a protein force
    /// field's topology excludes its own covalent neighbors from the
    /// non-bonded sum). `false` by default — the validation decks (e.g. the
    /// SHAKE trimers, whose bonded pairs deliberately feel the full LJ)
    /// rely on the current behavior.
    #[serde(default)]
    pub exclude_bonded: bool,
    #[serde(default)]
    pub dihedrals: Vec<DihedralSpec>,
}

/// One rigid distance constraint: `atoms = [i, j]`, fixed length `length`.
#[derive(Debug, Clone, Copy, Deserialize)]
pub struct ConstraintSpec {
    pub atoms: [usize; 2],
    pub length: f64,
}

/// A set of SHAKE/RATTLE distance constraints.
#[derive(Debug, Clone, Default, Deserialize)]
pub struct ConstraintsSpec {
    #[serde(default)]
    pub bonds: Vec<ConstraintSpec>,
}

fn default_element() -> String {
    "X".to_string()
}

fn default_mass() -> f64 {
    1.0
}

/// Mixing rule for cross-species pair parameters.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Deserialize, Default)]
#[serde(rename_all = "snake_case")]
pub enum MixingRuleSpec {
    /// Lorentz–Berthelot (default): σ_ij = (σ_i+σ_j)/2, ε_ij = √(ε_i·ε_j).
    #[default]
    LorentzBerthelot,
}

/// A per-type parameter entry for a multi-species potential.
#[derive(Debug, Clone, Deserialize)]
pub struct TypeSpec {
    pub sigma: f64,
    pub epsilon: f64,
}

/// An explicit per-pair override.
#[derive(Debug, Clone, Deserialize)]
pub struct OverrideSpec {
    pub types: [u32; 2],
    pub sigma: f64,
    pub epsilon: f64,
}

#[derive(Debug, Clone, Deserialize)]
#[serde(tag = "type", rename_all = "snake_case")]
pub enum ModelSpec {
    /// Lennard-Jones 12-6. Single-species (`sigma`/`epsilon`) or multi-species
    /// (`types` + optional `overrides`). If `types` is present it takes
    /// precedence; the single-pair fields are the backward-compatible default.
    LennardJones {
        #[serde(default)]
        sigma: Option<f64>,
        #[serde(default)]
        epsilon: Option<f64>,
        cutoff: f64,
        #[serde(default)]
        types: Vec<TypeSpec>,
        #[serde(default)]
        overrides: Vec<OverrideSpec>,
        #[serde(default)]
        mixing: MixingRuleSpec,
    },
    /// Buckingham `A·exp(−r/ρ) − B/r⁶`.
    Buckingham {
        #[serde(default)]
        a: Option<f64>,
        #[serde(default)]
        b: Option<f64>,
        #[serde(default)]
        rho: Option<f64>,
        cutoff: f64,
        #[serde(default)]
        types: Vec<BuckinghamTypeSpec>,
        #[serde(default)]
        overrides: Vec<BuckinghamOverrideSpec>,
    },
    /// Morse `D·[1 − exp(−α(r−r₀))]²`.
    Morse {
        #[serde(default)]
        d: Option<f64>,
        #[serde(default)]
        alpha: Option<f64>,
        #[serde(default)]
        r0: Option<f64>,
        cutoff: f64,
        #[serde(default)]
        types: Vec<MorseTypeSpec>,
        #[serde(default)]
        overrides: Vec<MorseOverrideSpec>,
    },
    /// Soft-sphere `ε·(σ/r)ⁿ`.
    SoftSphere {
        #[serde(default)]
        sigma: Option<f64>,
        #[serde(default)]
        epsilon: Option<f64>,
        #[serde(default)]
        n: Option<i32>,
        cutoff: f64,
        #[serde(default)]
        types: Vec<SoftSphereTypeSpec>,
        #[serde(default)]
        overrides: Vec<SoftSphereOverrideSpec>,
    },
    /// Tabulated potential from a text file (`r U F` per line).
    Tabulated {
        file: String,
        #[serde(default = "default_interp")]
        interpolation: InterpKindSpec,
    },
}

#[derive(Debug, Clone, Deserialize)]
pub struct BuckinghamTypeSpec {
    pub a: f64,
    pub b: f64,
    pub rho: f64,
}
#[derive(Debug, Clone, Deserialize)]
pub struct BuckinghamOverrideSpec {
    pub types: [u32; 2],
    pub a: f64,
    pub b: f64,
    pub rho: f64,
}
#[derive(Debug, Clone, Deserialize)]
pub struct MorseTypeSpec {
    pub d: f64,
    pub alpha: f64,
    pub r0: f64,
}
#[derive(Debug, Clone, Deserialize)]
pub struct MorseOverrideSpec {
    pub types: [u32; 2],
    pub d: f64,
    pub alpha: f64,
    pub r0: f64,
}
#[derive(Debug, Clone, Deserialize)]
pub struct SoftSphereTypeSpec {
    pub sigma: f64,
    pub epsilon: f64,
    pub n: i32,
}
#[derive(Debug, Clone, Deserialize)]
pub struct SoftSphereOverrideSpec {
    pub types: [u32; 2],
    pub sigma: f64,
    pub epsilon: f64,
    pub n: i32,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum InterpKindSpec {
    Linear,
}

fn default_interp() -> InterpKindSpec {
    InterpKindSpec::Linear
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
    /// Nosé–Hoover chain (true canonical NVT). M5.
    NoseHooverChain {
        dt: f64,
        target_temperature: f64,
        #[serde(default = "default_tau")]
        tau: f64,
        #[serde(default = "default_chain_len")]
        chain_len: usize,
    },
}

fn default_chain_len() -> usize {
    4
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

/// Pressure coupling (M5). `compressibility` is the isothermal compressibility
/// β_T (1/pressure); in reduced LJ units a value ~1 is typical.
#[derive(Debug, Clone, Deserialize)]
#[serde(tag = "type", rename_all = "snake_case")]
pub enum BarostatSpec {
    /// Berendsen weak-coupling barostat (good for equilibration).
    Berendsen {
        target_pressure: f64,
        #[serde(default = "default_tau")]
        tau: f64,
        #[serde(default = "default_compressibility")]
        compressibility: f64,
        #[serde(default)]
        anisotropic: bool,
    },
    /// Parrinello–Rahman barostat (samples the NPT ensemble).
    ParrinelloRahman {
        target_pressure: f64,
        #[serde(default = "default_tau")]
        tau: f64,
        /// Thermostat temperature for the piston mass coupling.
        target_temperature: f64,
        #[serde(default)]
        anisotropic: bool,
    },
}

fn default_compressibility() -> f64 {
    1.0
}

/// Ewald electrostatics parameters (M4). All lengths/wavenumbers are in the
/// user's input units (Å / 1·Å for `units = "real"`; reduced for
/// `units = "reduced"`); they are converted to reduced units at the I/O
/// boundary. Charges are supplied per-atom via the structure spec (Fcc `charge`,
/// `SpeciesSpec.charge`, or `RockSalt` cation/anion charges), in units of the
/// elementary charge `e`.
#[derive(Debug, Clone, Deserialize)]
pub struct EwaldSpec {
    /// Ewald screening parameter κ (1/length).
    pub kappa: f64,
    /// Real-space cutoff `r_c` (length).
    pub real_cutoff: f64,
    /// Reciprocal-space cutoff `k_c` (1/length).
    pub kspace_cutoff: f64,
}

/// Coulomb's constant in real (kcal·Å/mol) units: `e²/(4πε₀) =
/// 332.06371 kcal·Å/mol`. With charges in units of `e` and distances in Å, the
/// Coulomb energy is `332.06371 · q_i q_j / r` (kcal/mol). Used to derive the
/// reduced-charge factor `q* = q · √(COULOMB_KCAL_ANGSTROM · σ / ε)`.
pub const COULOMB_KCAL_ANGSTROM: f64 = 332.06371;

#[derive(Debug, Clone, Deserialize)]
pub struct RunSpec {
    pub steps: usize,
    #[serde(default = "default_initial_temp")]
    pub initial_temperature: f64,
    #[serde(default = "default_seed")]
    pub seed: u64,
    /// Verlet neighbor-list skin (same units as the cutoffs). The list is
    /// rebuilt only when any atom's minimum-image displacement since the last
    /// build exceeds `skin / 2`. `0` disables the skin (rebuild every step).
    /// Default `0.3` (LAMMPS's default skin).
    #[serde(default = "default_neighbor_skin")]
    pub neighbor_skin: f64,
    #[serde(default)]
    pub trajectory_file: Option<String>,
    #[serde(default = "default_traj_every")]
    pub trajectory_every: usize,
    #[serde(default)]
    pub log_file: Option<String>,
    #[serde(default = "default_log_every")]
    pub log_every: usize,
}

impl Default for RunSpec {
    /// The same defaults the TOML parser applies (`trajectory_every` 100,
    /// `log_every` 10, `neighbor_skin` 0.3, seed 42, T 1.0).
    fn default() -> Self {
        Self {
            steps: 100,
            initial_temperature: 1.0,
            seed: 42,
            neighbor_skin: 0.3,
            trajectory_file: None,
            trajectory_every: 100,
            log_file: None,
            log_every: 10,
        }
    }
}

fn default_initial_temp() -> f64 {
    1.0
}
fn default_neighbor_skin() -> f64 {
    0.3
}
fn default_traj_every() -> usize {
    100
}
fn default_log_every() -> usize {
    10
}

impl InputConfig {
    /// Parse a TOML config file. If `units = "real"`, every numeric field is
    /// converted to reduced units in-place using factors derived from the
    /// model's reference (σ, ε) and the structure's reference mass. The core
    /// therefore always runs in reduced units.
    pub fn from_path<P: AsRef<Path>>(path: P) -> Result<Self> {
        let text = std::fs::read_to_string(path)?;
        Self::from_str(&text)
    }

    /// Parse a config from a TOML string (I/O-boundary unit conversion applied).
    #[allow(clippy::should_implement_trait)]
    pub fn from_str(text: &str) -> Result<Self> {
        let mut config: InputConfig = toml::from_str(text)?;
        if config.units == UnitSystem::Real {
            let factors = config.real_unit_factors()?;
            config.convert_to_reduced(&factors);
            config.unit_factors = Some(factors);
        }
        // `step % log_every` panics on a zero divisor mid-run (and, on the MCP
        // path, strands the run in `Running` forever).
        config.validate()?;
        Ok(config)
    }

    /// Semantic validation of a parsed config (independent of unit conversion).
    pub(crate) fn validate(&self) -> Result<()> {
        if self.run.log_every == 0 {
            return Err(ModenaError::Config(
                "run.log_every must be >= 1 (0 would divide by zero)".into(),
            ));
        }
        if self.run.trajectory_every == 0 {
            return Err(ModenaError::Config(
                "run.trajectory_every must be >= 1 (0 would divide by zero)".into(),
            ));
        }
        if self.run.initial_temperature < 0.0 {
            return Err(ModenaError::Config(
                "run.initial_temperature must be >= 0".into(),
            ));
        }
        Ok(())
    }

    /// Derive the real→reduced conversion factors from the model's reference
    /// (σ, ε) and the structure's reference mass. The reference σ/ε are taken
    /// from the first model type (or the single-pair fields); the reference
    /// mass from the first species (or 1.0 amu).
    fn real_unit_factors(&self) -> Result<UnitFactors> {
        let (sigma, epsilon) = self.model.reference_sigma_epsilon()?;
        let mass_ref = self.structure.reference_mass();
        Ok(UnitFactors::new(sigma, epsilon, mass_ref))
    }

    /// Convert every numeric config field from real to reduced units in-place.
    fn convert_to_reduced(&mut self, f: &UnitFactors) {
        // Structure: lattice_constant (length), masses (mass), charges (charge).
        match &mut self.structure {
            StructureSpec::Fcc {
                lattice_constant,
                mass,
                charge,
                ..
            } => {
                *lattice_constant = f.to_reduced_length(*lattice_constant);
                if let Some(m) = mass {
                    *m = f.to_reduced_mass(*m);
                }
                if let Some(q) = charge {
                    *q = f.to_reduced_charge(*q);
                }
            }
            StructureSpec::Xyz {
                box_size,
                mass,
                charge,
                ..
            } => {
                for b in box_size.iter_mut() {
                    *b = f.to_reduced_length(*b);
                }
                if let Some(m) = mass {
                    *m = f.to_reduced_mass(*m);
                }
                if let Some(q) = charge {
                    *q = f.to_reduced_charge(*q);
                }
            }
            StructureSpec::RandomAlloy {
                lattice_constant,
                species,
                ..
            } => {
                *lattice_constant = f.to_reduced_length(*lattice_constant);
                for s in species.iter_mut() {
                    s.mass = f.to_reduced_mass(s.mass);
                    s.charge = f.to_reduced_charge(s.charge);
                }
            }
            StructureSpec::RockSalt {
                lattice_constant,
                cation_mass,
                anion_mass,
                cation_charge,
                anion_charge,
                ..
            } => {
                *lattice_constant = f.to_reduced_length(*lattice_constant);
                if let Some(m) = cation_mass {
                    *m = f.to_reduced_mass(*m);
                }
                if let Some(m) = anion_mass {
                    *m = f.to_reduced_mass(*m);
                }
                *cation_charge = f.to_reduced_charge(*cation_charge);
                *anion_charge = f.to_reduced_charge(*anion_charge);
            }
        }

        // Electrostatics: κ and k_c are 1/length (×σ), real_cutoff is length.
        if let Some(ew) = &mut self.electrostatics {
            let inv_len = 1.0 / f.to_reduced_length(1.0); // 1/Å → 1/reduced-length
            ew.kappa *= inv_len;
            ew.kspace_cutoff *= inv_len;
            ew.real_cutoff = f.to_reduced_length(ew.real_cutoff);
        }

        // Model: σ (length), ε (energy), cutoff (length), and per-type analogues.
        self.model.convert_to_reduced(f);

        // Integrator: dt (time), target_temperature (temp).
        match &mut self.integrator {
            IntegratorSpec::VelocityVerlet { dt } => {
                *dt = f.to_reduced_time(*dt);
            }
            IntegratorSpec::LangevinBaoab {
                dt,
                target_temperature,
                gamma,
                ..
            } => {
                *dt = f.to_reduced_time(*dt);
                *target_temperature = f.to_reduced_temperature(*target_temperature);
                // gamma has units 1/time.
                *gamma = f.to_reduced_time(1.0 / *gamma).recip();
            }
            IntegratorSpec::NoseHooverChain {
                dt,
                target_temperature,
                tau,
                ..
            } => {
                *dt = f.to_reduced_time(*dt);
                *target_temperature = f.to_reduced_temperature(*target_temperature);
                *tau = f.to_reduced_time(*tau);
            }
        }

        // Barostat: target_pressure (pressure), compressibility (1/pressure),
        // tau (time), target_temperature (temp).
        if let Some(b) = &mut self.barostat {
            match b {
                BarostatSpec::Berendsen {
                    target_pressure,
                    tau,
                    compressibility,
                    ..
                } => {
                    *target_pressure = f.to_reduced_pressure(*target_pressure);
                    *tau = f.to_reduced_time(*tau);
                    // κ has units 1/pressure → κ_reduced = κ_real · P_factor.
                    *compressibility *= f.to_real_pressure(1.0);
                }
                BarostatSpec::ParrinelloRahman {
                    target_pressure,
                    tau,
                    target_temperature,
                    ..
                } => {
                    *target_pressure = f.to_reduced_pressure(*target_pressure);
                    *tau = f.to_reduced_time(*tau);
                    *target_temperature = f.to_reduced_temperature(*target_temperature);
                }
            }
        }

        // Thermostat: target_temperature (temp), tau (time).
        if let Some(t) = &mut self.thermostat {
            match t {
                ThermostatSpec::Berendsen {
                    target_temperature,
                    tau,
                } => {
                    *target_temperature = f.to_reduced_temperature(*target_temperature);
                    *tau = f.to_reduced_time(*tau);
                }
                ThermostatSpec::VelocityRescaling {
                    target_temperature,
                    tau,
                    ..
                } => {
                    *target_temperature = f.to_reduced_temperature(*target_temperature);
                    *tau = f.to_reduced_time(*tau);
                }
                ThermostatSpec::None => {}
            }
        }

        // Run: initial_temperature (temp).
        self.run.initial_temperature = f.to_reduced_temperature(self.run.initial_temperature);

        // Bonded topology: r0/length (length), bond k (energy/length²), angle
        // k and dihedral k (energy). Angles/dihedral phases are in radians
        // (dimensionless).
        if let Some(topo) = &mut self.topology {
            for b in &mut topo.bonds {
                b.r0 = f.to_reduced_length(b.r0);
                b.k = f.to_reduced_energy(b.k) * f.to_reduced_length(1.0).powi(-2);
            }
            for a in &mut topo.angles {
                a.k = f.to_reduced_energy(a.k);
            }
            for d in &mut topo.dihedrals {
                for t in &mut d.terms {
                    t.k = f.to_reduced_energy(t.k);
                }
            }
        }
        if let Some(cons) = &mut self.constraints {
            for c in &mut cons.bonds {
                c.length = f.to_reduced_length(c.length);
            }
        }
    }
}

impl ModelSpec {
    /// Reference (σ, ε) for unit conversion, in real units. For multi-type
    /// models, the first type's parameters are the reference.
    fn reference_sigma_epsilon(&self) -> Result<(f64, f64)> {
        match self {
            ModelSpec::LennardJones {
                sigma,
                epsilon,
                types,
                ..
            } => {
                if let Some(t) = types.first() {
                    Ok((t.sigma, t.epsilon))
                } else if let (Some(s), Some(e)) = (sigma, epsilon) {
                    Ok((*s, *e))
                } else {
                    Err(ModenaError::Config(
                        "LennardJones needs either `types` or `sigma`/`epsilon`".into(),
                    ))
                }
            }
            ModelSpec::SoftSphere {
                sigma,
                epsilon,
                types,
                ..
            } => {
                if let Some(t) = types.first() {
                    Ok((t.sigma, t.epsilon))
                } else if let (Some(s), Some(e)) = (sigma, epsilon) {
                    Ok((*s, *e))
                } else {
                    Err(ModenaError::Config(
                        "SoftSphere needs either `types` or `sigma`/`epsilon`".into(),
                    ))
                }
            }
            // Buckingham/Morse have no single σ; use a 1 Å / 1 kcal/mol reference
            // (the absolute scale doesn't matter for a non-LJ potential as long
            // as the user supplies consistent units — the conversion is a no-op
            // when σ=ε=1, which is what reduced-mode already guarantees).
            ModelSpec::Buckingham { .. } | ModelSpec::Morse { .. } => Ok((1.0, 1.0)),
            ModelSpec::Tabulated { .. } => Ok((1.0, 1.0)),
        }
    }

    /// Convert model parameters from real to reduced units.
    fn convert_to_reduced(&mut self, f: &UnitFactors) {
        match self {
            ModelSpec::LennardJones {
                sigma,
                epsilon,
                cutoff,
                types,
                overrides,
                ..
            } => {
                *cutoff = f.to_reduced_length(*cutoff);
                if let Some(s) = sigma {
                    *s = f.to_reduced_length(*s);
                }
                if let Some(e) = epsilon {
                    *e = f.to_reduced_energy(*e);
                }
                for t in types.iter_mut() {
                    t.sigma = f.to_reduced_length(t.sigma);
                    t.epsilon = f.to_reduced_energy(t.epsilon);
                }
                for o in overrides.iter_mut() {
                    o.sigma = f.to_reduced_length(o.sigma);
                    o.epsilon = f.to_reduced_energy(o.epsilon);
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
                *cutoff = f.to_reduced_length(*cutoff);
                if let Some(v) = a {
                    *v = f.to_reduced_energy(*v);
                }
                if let Some(v) = rho {
                    *v = f.to_reduced_length(*v);
                }
                // B has units energy·length⁶.
                if let Some(v) = b {
                    *v = f.to_reduced_energy(*v) * f.to_reduced_length(1.0).powi(6);
                }
                for t in types.iter_mut() {
                    t.a = f.to_reduced_energy(t.a);
                    t.rho = f.to_reduced_length(t.rho);
                    t.b = f.to_reduced_energy(t.b) * f.to_reduced_length(1.0).powi(6);
                }
                for o in overrides.iter_mut() {
                    o.a = f.to_reduced_energy(o.a);
                    o.rho = f.to_reduced_length(o.rho);
                    o.b = f.to_reduced_energy(o.b) * f.to_reduced_length(1.0).powi(6);
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
                *cutoff = f.to_reduced_length(*cutoff);
                if let Some(v) = d {
                    *v = f.to_reduced_energy(*v);
                }
                if let Some(v) = alpha {
                    *v = f.to_reduced_length(1.0 / *v).recip();
                }
                if let Some(v) = r0 {
                    *v = f.to_reduced_length(*v);
                }
                for t in types.iter_mut() {
                    t.d = f.to_reduced_energy(t.d);
                    t.alpha = f.to_reduced_length(1.0 / t.alpha).recip();
                    t.r0 = f.to_reduced_length(t.r0);
                }
                for o in overrides.iter_mut() {
                    o.d = f.to_reduced_energy(o.d);
                    o.alpha = f.to_reduced_length(1.0 / o.alpha).recip();
                    o.r0 = f.to_reduced_length(o.r0);
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
                *cutoff = f.to_reduced_length(*cutoff);
                if let Some(v) = sigma {
                    *v = f.to_reduced_length(*v);
                }
                if let Some(v) = epsilon {
                    *v = f.to_reduced_energy(*v);
                }
                let _ = n;
                for t in types.iter_mut() {
                    t.sigma = f.to_reduced_length(t.sigma);
                    t.epsilon = f.to_reduced_energy(t.epsilon);
                }
                for o in overrides.iter_mut() {
                    o.sigma = f.to_reduced_length(o.sigma);
                    o.epsilon = f.to_reduced_energy(o.epsilon);
                }
            }
            ModelSpec::Tabulated { .. } => {
                // The table file is read in reduced units (caller's responsibility
                // to supply a reduced-unit table when using real units, since the
                // tabulated form has no single reference σ/ε).
            }
        }
    }
}

impl StructureSpec {
    /// Reference mass for unit conversion (amu). The first species' mass, or
    /// the explicit `mass` field, or 1.0 amu as a fallback.
    fn reference_mass(&self) -> f64 {
        match self {
            StructureSpec::Fcc { mass, .. } | StructureSpec::Xyz { mass, .. } => {
                mass.unwrap_or(1.0)
            }
            StructureSpec::RandomAlloy { species, .. } => {
                species.first().map(|s| s.mass).unwrap_or(1.0)
            }
            StructureSpec::RockSalt { cation_mass, .. } => cation_mass.unwrap_or(1.0),
        }
    }
}

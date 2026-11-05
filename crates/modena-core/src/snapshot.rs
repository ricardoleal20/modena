//! Strict native shifted-force LJ snapshots, with Angstrom/eV numerical inputs.
//! No simulation, velocities, thermodynamic unit conversion, or output files.

use serde::{Deserialize, Serialize};

use crate::{
    ForceModel, LennardJones, ModenaError, NeighborList, PairForceModel, Result, SimBox, System,
    Vec3, Virial,
};

pub const SNAPSHOT_PROTOCOL: &str = "modena-native-lj-v1";
pub const MAX_SNAPSHOT_BYTES: usize = 8 * 1024 * 1024;
pub const MAX_SNAPSHOT_ATOMS: usize = 100_000;

#[derive(Debug, Clone, Deserialize, Serialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct SnapshotUnits {
    pub length: String,
    pub energy: String,
    pub force: String,
}

impl Default for SnapshotUnits {
    fn default() -> Self {
        Self {
            length: "angstrom".into(),
            energy: "eV".into(),
            force: "eV/angstrom".into(),
        }
    }
}

#[derive(Debug, Clone, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct LjSnapshot {
    pub protocol: String,
    pub units: SnapshotUnits,
    pub atomic_numbers: Vec<u8>,
    pub positions_angstrom: Vec<[f64; 3]>,
    pub cell_angstrom: [[f64; 3]; 3],
    pub pbc: [bool; 3],
    pub sigma_angstrom: f64,
    pub epsilon_ev: f64,
    pub cutoff_angstrom: f64,
    pub skin_angstrom: f64,
}

#[derive(Debug, Serialize)]
pub struct SnapshotEvaluation {
    pub protocol: &'static str,
    pub units: SnapshotUnits,
    pub energy_ev: f64,
    pub forces_ev_per_angstrom: Vec<[f64; 3]>,
}

fn invalid(message: &str) -> ModenaError {
    ModenaError::Input(message.into())
}

/// Evaluate through the existing native pair path. All lengths are Angstrom,
/// all energies eV; dimensional consistency gives forces in eV/Angstrom.
pub fn evaluate_lj_snapshot(input: &LjSnapshot) -> Result<SnapshotEvaluation> {
    if input.protocol != SNAPSHOT_PROTOCOL || input.units != SnapshotUnits::default() {
        return Err(invalid("unsupported snapshot protocol or units"));
    }
    let n = input.atomic_numbers.len();
    if n == 0 || n > MAX_SNAPSHOT_ATOMS || input.positions_angstrom.len() != n {
        return Err(invalid(
            "snapshot needs 1..=100000 atoms and matching positions",
        ));
    }
    let z = input.atomic_numbers[0];
    if !(1..=118).contains(&z) || input.atomic_numbers.iter().any(|&a| a != z) {
        return Err(invalid(
            "snapshot supports homogeneous atomic numbers 1..=118",
        ));
    }
    if input.pbc != [true; 3] {
        return Err(invalid("snapshot requires full periodicity"));
    }
    for (i, row) in input.cell_angstrom.iter().enumerate() {
        for (j, &value) in row.iter().enumerate() {
            if !value.is_finite() || (i == j && value <= 0.) || (i != j && value != 0.) {
                return Err(invalid(
                    "snapshot requires a positive finite diagonal orthorhombic cell",
                ));
            }
        }
    }
    if input
        .positions_angstrom
        .iter()
        .flatten()
        .any(|v| !v.is_finite())
    {
        return Err(invalid("snapshot positions must be finite"));
    }
    if !input.sigma_angstrom.is_finite()
        || input.sigma_angstrom <= 0.
        || !input.epsilon_ev.is_finite()
        || input.epsilon_ev < 0.
        || !input.cutoff_angstrom.is_finite()
        || input.cutoff_angstrom <= 0.
        || !input.skin_angstrom.is_finite()
        || input.skin_angstrom < 0.
    {
        return Err(invalid("invalid finite LJ parameters or skin"));
    }
    let lengths = std::array::from_fn(|i| input.cell_angstrom[i][i]);
    let radius = input.cutoff_angstrom + input.skin_angstrom;
    if !radius.is_finite()
        || lengths.iter().any(|&l| radius > 0.5 * l)
        || !radius.powi(2).is_finite()
        || radius.powi(2) == 0.
    {
        return Err(invalid(
            "cutoff plus skin must be finite and at most half each box length",
        ));
    }
    // The shared native grid is dense. Bound its allocation before building,
    // including extreme finite cells and tiny cutoffs supplied over JSON.
    let cells = lengths.iter().try_fold(1usize, |count, &l| {
        let dim = (l / radius).floor();
        if dim > 1_000_000. {
            return None;
        }
        count
            .checked_mul((dim as usize).max(1))
            .filter(|&c| c <= 1_000_000)
    });
    if cells.is_none() {
        return Err(invalid("snapshot neighbor grid exceeds 1000000 cells"));
    }
    let mut system = System::new(SimBox::from_array(lengths));
    system.resize(n);
    for (i, p) in input.positions_angstrom.iter().enumerate() {
        system.set_pos(i, Vec3(*p));
    }
    system.wrap_positions();
    if system.positions.iter().any(|v| !v.is_finite()) {
        return Err(invalid("snapshot position wrapping overflow"));
    }
    let mut neighbors = NeighborList::with_skin(input.cutoff_angstrom, input.skin_angstrom);
    neighbors.build(&system);
    for &(i, j) in neighbors.pairs() {
        if system
            .box_size
            .min_image(system.pos(i) - system.pos(j))
            .norm2()
            < 1e-14
        {
            return Err(invalid("coincident or singular minimum-image atoms"));
        }
    }
    let model = PairForceModel::new(LennardJones::new(
        input.sigma_angstrom,
        input.epsilon_ev,
        input.cutoff_angstrom,
    )?);
    let energy = model.compute(&mut system, &mut Virial::default(), &mut neighbors)?;
    if !energy.is_finite() || system.forces.iter().any(|v| !v.is_finite()) {
        return Err(ModenaError::Physics(
            "nonfinite native LJ snapshot output".into(),
        ));
    }
    Ok(SnapshotEvaluation {
        protocol: SNAPSHOT_PROTOCOL,
        units: SnapshotUnits::default(),
        energy_ev: energy,
        forces_ev_per_angstrom: system
            .forces
            .chunks_exact(3)
            .map(|f| [f[0], f[1], f[2]])
            .collect(),
    })
}

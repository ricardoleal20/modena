//! Frozen example-only records. All measured quantities are checked before reduction.
use modena_core::stress::{ConstraintSolverOptions, TensorAvailability};
use modena_core::{Constraints, SimBox, System, Vec3};
use serde::{Deserialize, Serialize};
use std::{
    collections::BTreeMap,
    fs::{File, OpenOptions},
    io::Write,
    path::Path,
};
pub type R<T> = Result<T, String>;

/// Completed numerical work only. A killed process can leave an unknown suffix.
#[derive(Default)]
pub struct Work {
    system_steps: u64,
    molecule_steps: u64,
    receipt: Option<std::path::PathBuf>,
}
impl Work {
    pub fn begin() -> Self {
        Self {
            receipt: std::env::var_os("B07_WORK_PATH").map(Into::into),
            ..Self::default()
        }
    }
    pub fn step(&mut self, molecules: usize) -> R<()> {
        self.system_steps += 1;
        self.molecule_steps += molecules as u64;
        if self.system_steps.is_multiple_of(256) {
            self.persist(false)?;
        }
        Ok(())
    }
    pub fn persist(&self, complete: bool) -> R<()> {
        if let Some(p) = &self.receipt {
            let temp = p.with_extension("tmp");
            let row = serde_json::json!({"system_steps":self.system_steps,"molecule_steps":self.molecule_steps,"complete":complete});
            std::fs::write(&temp, serde_json::to_vec(&row).map_err(err)?).map_err(err)?;
            std::fs::rename(temp, p).map_err(err)?;
        }
        Ok(())
    }
}
pub fn err(e: impl std::fmt::Display) -> String {
    e.to_string()
}
pub fn finite(xs: &[f64]) -> R<()> {
    if xs.iter().all(|x| x.is_finite()) {
        Ok(())
    } else {
        Err("nonfinite validation data".into())
    }
}
pub fn tensor(t: &TensorAvailability) -> R<[f64; 6]> {
    let TensorAvailability::Complete(t) = t else {
        return Err("incomplete mechanical tensor".into());
    };
    finite(&t.0)?;
    Ok(t.0)
}
pub fn write_new(path: &Path, value: &impl Serialize) -> R<()> {
    if path.try_exists().map_err(err)? {
        return Err(format!("output already exists: {}", path.display()));
    }
    // A failed serialization leaves its bytes beside the intended path. Linking
    // the finished sibling publishes atomically without replacing an earlier result.
    let mut partial = path.as_os_str().to_os_string();
    partial.push(".partial");
    let staged = Path::new(&partial);
    let mut file = OpenOptions::new()
        .create_new(true)
        .write(true)
        .open(staged)
        .map_err(err)?;
    serde_json::to_writer(&mut file, value).map_err(err)?;
    file.write_all(b"\n").map_err(err)?;
    file.sync_all().map_err(err)?;
    std::fs::hard_link(staged, path).map_err(err)?;
    std::fs::remove_file(staged).map_err(err)
}
pub fn read<T: serde::de::DeserializeOwned>(path: &Path) -> R<T> {
    serde_json::from_reader(File::open(path).map_err(err)?).map_err(err)
}
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum CaseStatus {
    Prepared,
    Passed,
    Failed,
    Inconclusive,
    Unsupported,
    Pending,
    BudgetExceeded,
}
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct GateResult {
    pub name: String,
    pub observed: Vec<f64>,
    pub expected: Vec<f64>,
    pub absolute_error: Vec<f64>,
    pub tolerance: Vec<f64>,
    pub status: CaseStatus,
}
pub fn gate(name: &str, obs: &[f64], expected: &[f64], tolerance: &[f64]) -> R<GateResult> {
    if obs.len() != expected.len()
        || obs.len() != tolerance.len()
        || obs.is_empty()
        || tolerance.iter().any(|x| *x < 0.)
    {
        return Err(format!("gate shape or tolerance: {name}"));
    }
    for v in [obs, expected, tolerance] {
        finite(v)?;
    }
    let diff: Vec<_> = obs
        .iter()
        .zip(expected)
        .map(|(a, b)| (a - b).abs())
        .collect();
    finite(&diff)?;
    let status = if diff.iter().zip(tolerance).all(|(e, t)| e <= t) {
        CaseStatus::Passed
    } else {
        CaseStatus::Failed
    };
    Ok(GateResult {
        name: name.into(),
        observed: obs.into(),
        expected: expected.into(),
        absolute_error: diff,
        tolerance: tolerance.into(),
        status,
    })
}
pub fn scaled(name: &str, observed: &[f64], expected: &[f64], tol: f64) -> R<GateResult> {
    gate(
        name,
        observed,
        expected,
        &expected
            .iter()
            .map(|x| tol * (1. + x.abs()))
            .collect::<Vec<_>>(),
    )
}
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Sample {
    pub seed: u64,
    pub step: u64,
    pub time: f64,
    pub dt: f64,
    pub endpoint: [f64; 6],
    pub interval: [f64; 6],
    pub com_energy: f64,
    pub rotation_energy: f64,
    pub orientation_second: [f64; 6],
}
impl Sample {
    pub fn validate(&self) -> R<()> {
        finite(&self.endpoint)?;
        finite(&self.interval)?;
        finite(&self.orientation_second)?;
        finite(&[self.time, self.dt, self.com_energy, self.rotation_energy])?;
        if self.time < 0. || self.dt <= 0. || self.com_energy < 0. || self.rotation_energy < 0. {
            return Err("invalid sample".into());
        }
        Ok(())
    }
}
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Estimate {
    pub mean: f64,
    pub se: f64,
    pub ci95: [f64; 2],
    pub n_eff: f64,
    pub block_samples: usize,
    pub retained_blocks: usize,
}
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ValidationFixture {
    pub id: String,
    pub positions: Vec<[f64; 3]>,
    pub velocities: Vec<[f64; 3]>,
    pub masses: Vec<f64>,
    pub charges: Vec<f64>,
    pub cell: [[f64; 3]; 3],
    pub constraints: Constraints,
    pub kappa: Option<f64>,
    pub real_cutoff: Option<f64>,
    pub reciprocal_cutoff: Option<f64>,
    pub pair_scales: Vec<(usize, usize, f64)>,
}
impl ValidationFixture {
    pub fn from_system(id: String, s: &System, constraints: Constraints) -> Self {
        Self {
            id,
            positions: (0..s.n()).map(|i| s.pos(i).0).collect(),
            velocities: (0..s.n()).map(|i| s.vel(i).0).collect(),
            masses: s.masses.clone(),
            charges: s.charges.clone(),
            cell: s.box_size.basis().map(|v| v.0),
            constraints,
            kappa: None,
            real_cutoff: None,
            reciprocal_cutoff: None,
            pair_scales: vec![],
        }
    }
    pub fn system(&self) -> R<System> {
        let n = self.masses.len();
        if n == 0
            || self.velocities.len() != n
            || self.positions.len() != n
            || self.charges.len() != n
            || self.masses.iter().any(|m| *m <= 0.)
        {
            return Err("fixture shape/masses".into());
        }
        for a in self
            .positions
            .iter()
            .chain(&self.velocities)
            .chain(&self.cell)
        {
            finite(a)?
        }
        finite(&self.masses)?;
        finite(&self.charges)?;
        for v in [self.kappa, self.real_cutoff, self.reciprocal_cutoff]
            .into_iter()
            .flatten()
        {
            finite(&[v])?;
            if v <= 0. {
                return Err("nonpositive electrostatic parameter".into());
            }
        }
        let mut seen = std::collections::BTreeSet::new();
        for &(i, j, v) in &self.pair_scales {
            if i >= n || j >= n || i >= j || !seen.insert((i, j)) || !(0. ..=1.).contains(&v) {
                return Err("pair scale".into());
            }
        }
        let [a, b, c] = self.cell.map(Vec3);
        if !a.dot(b.cross(c)).is_finite() || a.dot(b.cross(c)) <= 0. {
            return Err("invalid oriented cell".into());
        }
        if [a.0[1], a.0[2], b.0[0], b.0[2], c.0[0], c.0[1]] != [0.; 6]
            || [a.0[0], b.0[1], c.0[2]].iter().any(|v| *v <= 0.)
        {
            return Err("frozen fixture requires orthorhombic cell".into());
        }
        let mut s = System::new(SimBox::from_array([a.0[0], b.0[1], c.0[2]]));
        s.resize(n);
        s.masses.clone_from(&self.masses);
        s.charges.clone_from(&self.charges);
        for i in 0..n {
            s.set_pos(i, Vec3(self.positions[i]));
            s.set_vel(i, Vec3(self.velocities[i]));
            s.species[i] = u32::from(self.kappa.is_some() && self.id != "pair-oo" && i % 3 != 0);
        }
        self.constraints.validate(&s).map_err(err)?;
        Ok(s)
    }
}
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct WorkLimits {
    pub wall_seconds: u64,
    pub max_workers: usize,
    pub rust_rss_bytes: u64,
    pub reference_rss_bytes: u64,
    pub artifact_bytes: u64,
    pub molecule_steps: u64,
    pub projected_draws: u64,
    pub reference_jobs: usize,
}
impl WorkLimits {
    pub const FROZEN: Self = Self {
        wall_seconds: 7200,
        max_workers: 4,
        rust_rss_bytes: 512 << 20,
        reference_rss_bytes: 1 << 30,
        artifact_bytes: 5 << 30,
        molecule_steps: 400_000_000,
        projected_draws: 1_000_000,
        reference_jobs: 24,
    };
}
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ValidationProtocol {
    pub version: u32,
    pub solver: ConstraintSolverOptions,
    pub strain_eps: Vec<f64>,
    pub gas_seeds: Vec<u64>,
    pub gas_dt: Vec<f64>,
    pub burn_time: f64,
    pub measure_times: Vec<f64>,
    pub sample_interval: f64,
    pub fixtures: Vec<ValidationFixture>,
    pub limits: WorkLimits,
    pub specification_sha256: String,
    pub source_sha256: BTreeMap<String, String>,
    pub snapshot_sha256: BTreeMap<String, String>,
}
impl ValidationProtocol {
    pub fn validate(&self) -> R<()> {
        if self.version != 1
            || self.solver != ConstraintSolverOptions::TIGHT
            || self.strain_eps != [1e-4, 3e-5, 1e-5, 3e-6]
            || self.gas_seeds != (1001..=1008).collect::<Vec<_>>()
            || self.gas_dt != [0.04, 0.02, 0.01, 0.005]
            || self.measure_times != [400., 800., 1600.]
            || self.burn_time != 50.
            || self.sample_interval != 0.1
            || self.limits != WorkLimits::FROZEN
        {
            return Err("protocol constants changed".into());
        }
        let mut ids = std::collections::BTreeSet::new();
        for f in &self.fixtures {
            if !ids.insert(&f.id) {
                return Err("duplicate fixture".into());
            }
            f.system()?;
        }
        let valid_digest = |s: &str| s.len() == 64 && s.bytes().all(|b| b.is_ascii_hexdigit());
        if !valid_digest(&self.specification_sha256)
            || !self.source_sha256.iter().any(|(k, v)| {
                k.ends_with("/docs/superpowers/specs/2026-11-06-b07-tensor-pressure-design.md")
                    && v == &self.specification_sha256
            })
            || self
                .source_sha256
                .values()
                .chain(self.snapshot_sha256.values())
                .any(|v| !valid_digest(v))
        {
            return Err("invalid protocol digest".into());
        }
        if ids.len() != 19 || self.source_sha256.is_empty() || self.snapshot_sha256.is_empty() {
            return Err("empty freeze inventory".into());
        }
        Ok(())
    }
}
#[derive(Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ExecutionProvenance {
    pub source_head: String,
    pub compiler: String,
    pub build_profile: String,
    pub features: Vec<String>,
    pub command: Vec<String>,
    pub protocol_sha256: String,
    pub binary_sha256: String,
    pub external_source_revision: Option<String>,
    pub external_packages: Vec<String>,
    pub threads: usize,
    pub started_utc: String,
    pub ended_utc: String,
    pub elapsed_seconds: f64,
    pub peak_rss_bytes: u64,
    pub molecule_steps: u64,
    pub projected_draws: u64,
    pub reference_jobs: usize,
}
#[derive(Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ValidationReport {
    pub protocol: ValidationProtocol,
    pub provenance: ExecutionProvenance,
    pub gates: Vec<GateResult>,
    pub raw_files: BTreeMap<String, String>,
    pub status: CaseStatus,
}

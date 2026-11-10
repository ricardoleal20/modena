//! Offline calibration records, pure preparation, and independent reference evaluation.
#![allow(dead_code)] // Wire records are consumed by later offline phases.
use std::{
    collections::BTreeMap,
    io::Write,
    path::{Path, PathBuf},
    process::Command,
};
pub mod allocator;
mod candidates;
mod cost;
mod evaluator;
mod fixtures;
mod metrics;
pub mod physical;
pub mod profile;
mod selection;
#[allow(unused_imports)]
pub use cost::{
    begin_selection, begin_selection_with_reference_observer, cost_smoke,
    freeze_prepared_selection, measure_feasible_candidates,
    measure_feasible_candidates_with_observer, measure_ordinary, measure_ordinary_with_observer,
    ConfirmationOutcome, CostEvent, CostStage, FrozenSelection, OrdinaryCostOutcome,
    PreparedReferences, TimingBudget,
};
#[allow(unused_imports)]
pub use selection::ordinary_timing_summary;
#[allow(unused_imports)]
pub use selection::{finish_confirmation, select, selected_fragment};
mod reference;
#[allow(unused_imports)]
pub use candidates::{
    candidate_grid, evaluate_prepared_candidate, evaluate_prepared_candidate_with_observer,
    feasible_targets, preflight_selection_sweep, CandidateBudget, CandidateEvent, CandidateOutcome,
    TranslationCertificate,
};
#[allow(unused_imports)]
pub use evaluator::{evaluate_spme, SpmeEvaluator};
#[allow(unused_imports)]
// Public pure APIs are also used by the integration test and later phases.
pub use fixtures::{frozen_protocol, materialize, validate_frozen, validate_snapshot};
#[allow(unused_imports)]
pub use metrics::{metrics, passes_target};
#[cfg(test)]
#[allow(unused_imports)]
pub use reference::converged_reference_test;
#[allow(unused_imports)]
pub use reference::{
    converged_prepared_reference_with_observer, converged_reference,
    converged_reference_with_budget, reference, reference_with_budget, ReferenceAudit,
    ReferenceBudget, ReferenceEvent, ReferenceOrigin, ReferenceOutcome, ReferenceUsage,
};
pub const SCHEMA: &str = "modena-spme-calibration-v1";
#[derive(Clone, Debug, PartialEq, serde::Serialize, serde::Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Snapshot {
    pub id: String,
    pub role: FixtureRole,
    pub family: String,
    pub units: String,
    pub pbc: [bool; 3],
    pub seed: Option<u64>,
    #[serde(
        serialize_with = "serialize_finite",
        deserialize_with = "deserialize_finite"
    )]
    pub amplitude: Option<f64>,
    #[serde(
        serialize_with = "serialize_finite",
        deserialize_with = "deserialize_finite"
    )]
    pub cell: [f64; 3],
    #[serde(
        serialize_with = "serialize_finite",
        deserialize_with = "deserialize_finite"
    )]
    pub positions: Vec<[f64; 3]>,
    #[serde(
        serialize_with = "serialize_finite",
        deserialize_with = "deserialize_finite"
    )]
    pub charges: Vec<f64>,
    pub scales: Vec<Scale>,
}
#[derive(Clone, Debug, PartialEq, serde::Serialize, serde::Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Scale {
    pub i: usize,
    pub j: usize,
    #[serde(
        serialize_with = "serialize_finite",
        deserialize_with = "deserialize_finite"
    )]
    pub coulomb: f64,
}
#[derive(Clone, Debug, PartialEq, serde::Serialize, serde::Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Parameters {
    #[serde(
        serialize_with = "serialize_finite",
        deserialize_with = "deserialize_finite"
    )]
    pub kappa: f64,
    #[serde(
        serialize_with = "serialize_finite",
        deserialize_with = "deserialize_finite"
    )]
    pub real_cutoff: f64,
    #[serde(
        serialize_with = "serialize_finite",
        deserialize_with = "deserialize_finite"
    )]
    pub kspace_cutoff: f64,
}
#[derive(Clone, Debug, PartialEq, serde::Serialize, serde::Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Evaluation {
    #[serde(
        serialize_with = "serialize_finite",
        deserialize_with = "deserialize_finite"
    )]
    pub energy: f64,
    #[serde(
        serialize_with = "serialize_finite",
        deserialize_with = "deserialize_finite"
    )]
    pub strain_w: f64,
    #[serde(
        serialize_with = "serialize_finite",
        deserialize_with = "deserialize_finite"
    )]
    pub forces: Vec<[f64; 3]>,
}
#[derive(Clone, Debug, PartialEq, serde::Serialize, serde::Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ReferenceSettings {
    #[serde(
        serialize_with = "serialize_finite",
        deserialize_with = "deserialize_finite"
    )]
    pub kappa: f64,
    #[serde(
        serialize_with = "serialize_finite",
        deserialize_with = "deserialize_finite"
    )]
    pub real_cutoff: f64,
    #[serde(
        serialize_with = "serialize_finite",
        deserialize_with = "deserialize_finite"
    )]
    pub kspace_cutoff: f64,
}

#[derive(Clone, Debug, PartialEq, serde::Serialize, serde::Deserialize)]
pub enum FixtureRole {
    Selection,
    Confirmation,
    RegressionOnly,
}
#[derive(Clone, Debug, PartialEq, serde::Serialize, serde::Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Limits {
    pub max_atoms: usize,                 // 256
    pub max_mesh_points: usize,           // 1_048_576
    pub max_candidates: usize,            // 126 per family
    pub max_accuracy_calls: u64,          // 25_000 per family, selection+confirmation
    pub max_timing_calls: u64,            // 20_000 per family, includes cold/rebuild/warmup
    pub max_profile_calls: u64,           // 1_000 per family
    pub max_real_terms_per_call: u64,     // 100_000_000
    pub max_atom_k_terms_per_call: u64,   // 100_000_000
    pub max_real_terms_per_family: u64,   // 5_000_000_000
    pub max_atom_k_terms_per_family: u64, // 500_000_000
    pub max_owned_heap_bytes: u64,        // 268_435_456 (256 MiB)
    pub max_peak_rss_bytes: u64,          // 536_870_912 (512 MiB), acceptance cap
    pub max_family_wall_seconds: u64,     // 1_800, cumulative verify/sweep/confirm/profile
    pub max_threads: usize,               // 4
}
#[derive(Clone, Debug, PartialEq, serde::Serialize, serde::Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Tolerances {
    #[serde(
        serialize_with = "serialize_finite",
        deserialize_with = "deserialize_finite"
    )]
    pub reference_force_rms_abs_units: f64, // 1e-8 F_unit
    #[serde(
        serialize_with = "serialize_finite",
        deserialize_with = "deserialize_finite"
    )]
    pub reference_force_target_fraction: f64, // 0.01 times tightest target times D
    #[serde(
        serialize_with = "serialize_finite",
        deserialize_with = "deserialize_finite"
    )]
    pub reference_force_max_units: f64, // 1e-7 F_unit
    #[serde(
        serialize_with = "serialize_finite",
        deserialize_with = "deserialize_finite"
    )]
    pub reference_energy_per_atom_units: f64, // 1e-9 E_unit
    #[serde(
        serialize_with = "serialize_finite",
        deserialize_with = "deserialize_finite"
    )]
    pub reference_w_per_atom_units: f64, // 1e-9 E_unit
    #[serde(
        serialize_with = "serialize_finite",
        deserialize_with = "deserialize_finite"
    )]
    pub relative_denominator_floor: f64, // 1e-8 F_unit
    #[serde(
        serialize_with = "serialize_finite",
        deserialize_with = "deserialize_finite"
    )]
    pub candidate_force_rms_abs_factor: f64, // 1 times target F_unit
    #[serde(
        serialize_with = "serialize_finite",
        deserialize_with = "deserialize_finite"
    )]
    pub candidate_force_max_factor: f64, // 10 times target F_unit
    #[serde(
        serialize_with = "serialize_finite",
        deserialize_with = "deserialize_finite"
    )]
    pub candidate_energy_factor: f64, // 1 times target E_unit per atom
    #[serde(
        serialize_with = "serialize_finite",
        deserialize_with = "deserialize_finite"
    )]
    pub candidate_w_factor: f64, // 10 times target E_unit per atom
    #[serde(
        serialize_with = "serialize_finite",
        deserialize_with = "deserialize_finite"
    )]
    pub cost_tie_fraction: f64, // 0.05
}
#[derive(Clone, Debug, PartialEq, serde::Serialize, serde::Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Protocol {
    pub schema: String,      // modena-spme-calibration-v1
    pub fixture_set: String, // c2-spme-v1
    pub family: String,      // cubic-64, ortho-64, special-3, cubic-216, ortho-216
    pub replicas: u32,       // 2 or 3; special-3 prepared alongside replicas 2
    #[serde(
        serialize_with = "serialize_finite",
        deserialize_with = "deserialize_finite"
    )]
    pub a: f64,
    #[serde(
        serialize_with = "serialize_finite",
        deserialize_with = "deserialize_finite"
    )]
    pub q_unit: f64, // 2.2, 1.0
    #[serde(
        serialize_with = "serialize_finite",
        deserialize_with = "deserialize_finite"
    )]
    pub targets: Vec<f64>, // [1e-3,1e-4]
    #[serde(
        serialize_with = "serialize_finite",
        deserialize_with = "deserialize_finite"
    )]
    pub kappa_a: Vec<f64>, // [1.2,1.6,2.0,2.4,3.0,3.6]
    #[serde(
        serialize_with = "serialize_finite",
        deserialize_with = "deserialize_finite"
    )]
    pub cutoff_lmin: Vec<f64>, // [.30,.40,.49]
    #[serde(
        serialize_with = "serialize_finite",
        deserialize_with = "deserialize_finite"
    )]
    pub kcut_a: Vec<f64>, // [8,12,16,24,32,48,64]
    #[serde(
        serialize_with = "serialize_finite",
        deserialize_with = "deserialize_finite"
    )]
    pub reference_ladder: Vec<[f64; 2]>, // [R/a,K*a], four or prechosen six levels
    #[serde(
        serialize_with = "serialize_finite",
        deserialize_with = "deserialize_finite"
    )]
    pub reference_kappa_a: f64, // 1.1
    #[serde(
        serialize_with = "serialize_finite",
        deserialize_with = "deserialize_finite"
    )]
    pub crosscheck_kappa_a: [f64; 2], // [.9,1.3]
    #[serde(
        serialize_with = "serialize_finite",
        deserialize_with = "deserialize_finite"
    )]
    pub translations: Vec<[f64; 3]>, // exact seven rows in section 2
    pub limits: Limits,
    pub tolerances: Tolerances,
    pub threads: usize, // 1 or 4
    pub warmups: usize,
    pub batches: usize,
    pub batch_calls: usize, // 5,5,20
}
#[derive(Clone, Debug, PartialEq, serde::Serialize, serde::Deserialize)]
#[serde(deny_unknown_fields)]
pub struct WorkCounts {
    pub real_terms: u64,
    pub atom_k_terms: u64,
}
#[derive(Clone, Debug, PartialEq, serde::Serialize, serde::Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ErrorMetrics {
    #[serde(
        serialize_with = "serialize_finite",
        deserialize_with = "deserialize_finite"
    )]
    pub force_rms_abs: f64,
    #[serde(
        serialize_with = "serialize_finite",
        deserialize_with = "deserialize_finite"
    )]
    pub force_rms_rel: Option<f64>,
    #[serde(
        serialize_with = "serialize_finite",
        deserialize_with = "deserialize_finite"
    )]
    pub reference_force_rms: f64,
    #[serde(
        serialize_with = "serialize_finite",
        deserialize_with = "deserialize_finite"
    )]
    pub max_atom_force_error: f64,
    #[serde(
        serialize_with = "serialize_finite",
        deserialize_with = "deserialize_finite"
    )]
    pub max_component_force_error: f64,
    #[serde(
        serialize_with = "serialize_finite",
        deserialize_with = "deserialize_finite"
    )]
    pub energy_per_atom_signed: f64,
    #[serde(
        serialize_with = "serialize_finite",
        deserialize_with = "deserialize_finite"
    )]
    pub energy_per_atom_abs: f64,
    #[serde(
        serialize_with = "serialize_finite",
        deserialize_with = "deserialize_finite"
    )]
    pub w_per_atom_signed: f64,
    #[serde(
        serialize_with = "serialize_finite",
        deserialize_with = "deserialize_finite"
    )]
    pub w_per_atom_abs: f64,
    #[serde(
        serialize_with = "serialize_finite",
        deserialize_with = "deserialize_finite"
    )]
    pub net_force: [f64; 3],
    #[serde(
        serialize_with = "serialize_finite",
        deserialize_with = "deserialize_finite"
    )]
    pub net_force_norm: f64,
}
#[derive(Clone, Debug, PartialEq, serde::Serialize, serde::Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ReferenceLevel {
    pub settings: ReferenceSettings,
    pub evaluation: Evaluation,
    pub delta_from_previous: Option<ErrorMetrics>,
    pub work: WorkCounts,
}
#[derive(Clone, Debug, PartialEq, serde::Serialize, serde::Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ReferenceRecord {
    pub snapshot_id: String,
    pub snapshot_sha256: String,
    pub levels: Vec<ReferenceLevel>,
    pub crosschecks: Vec<ReferenceLevel>,
    pub chosen_level: Option<usize>,
    pub status: Status,
}
#[derive(Clone, Debug, PartialEq, serde::Serialize, serde::Deserialize)]
#[serde(deny_unknown_fields)]
pub struct MetricRow {
    pub snapshot_id: String,
    pub shift_index: usize,
    pub metrics: ErrorMetrics,
    #[serde(
        serialize_with = "serialize_finite",
        deserialize_with = "deserialize_finite"
    )]
    pub passing_targets: Vec<f64>,
    pub status: Status,
}
#[derive(Clone, Debug, PartialEq, serde::Serialize, serde::Deserialize)]
#[serde(deny_unknown_fields)]
pub struct TimingRecord {
    pub reset_included: bool, // true for ordinary primary batches
    pub cold_ns: u64,
    pub rebuild_ns: u64,
    pub batch_total_ns: Vec<u64>,
    pub calls_per_batch: usize,
    #[serde(
        serialize_with = "serialize_finite",
        deserialize_with = "deserialize_finite"
    )]
    pub median_ns_per_call: f64,
    #[serde(
        serialize_with = "serialize_finite",
        deserialize_with = "deserialize_finite"
    )]
    pub q1_ns_per_call: f64,
    #[serde(
        serialize_with = "serialize_finite",
        deserialize_with = "deserialize_finite"
    )]
    pub q3_ns_per_call: f64,
}
#[derive(Clone, Debug, PartialEq, serde::Serialize, serde::Deserialize)]
#[serde(deny_unknown_fields)]
pub struct CandidateReport {
    pub id: usize,
    pub family: String,
    pub parameters: Parameters,
    pub actual_mesh: [usize; 3],
    #[serde(
        serialize_with = "serialize_finite",
        deserialize_with = "deserialize_finite"
    )]
    pub aliases: Vec<f64>,
    pub rows: Vec<MetricRow>,
    pub timing: Option<TimingRecord>,
    pub status: Status,
}
#[derive(Clone, Debug, PartialEq, serde::Serialize, serde::Deserialize)]
#[serde(deny_unknown_fields)]
pub struct TargetDecision {
    #[serde(
        serialize_with = "serialize_finite",
        deserialize_with = "deserialize_finite"
    )]
    pub target: f64,
    pub selected_id: Option<usize>,
    pub selection_status: Status,
    pub confirmation_rows: Vec<MetricRow>,
    pub confirmation_status: Status,
}
#[derive(Clone, Debug, PartialEq, serde::Serialize, serde::Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PhaseRecord {
    pub candidate_id: usize,
    pub snapshot_id: String,
    pub threads: usize,
    pub cold: bool,
    pub calls: usize,
    pub phases_ns: Vec<[u64; 12]>,
}
#[derive(Clone, Debug, PartialEq, serde::Serialize, serde::Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Provenance {
    pub git_head: String,
    pub source_sha256: std::collections::BTreeMap<String, String>,
    pub protocol_sha256: String,
    pub snapshot_sha256: String,
    pub reference_algorithm: String, // all-image-full-pm-k-v1
    pub rustc_verbose: String,
    pub cargo_version: String,
    pub os: String,
    pub cpu: String,
    pub threads: usize,
    pub profile: String,
    pub features: Vec<String>,
    pub command: Vec<String>,
    pub utc_started: String,
    pub family_elapsed_ns: u64,
    pub measured_peak_rss_bytes: Option<u64>,
}
#[derive(Clone, Debug, PartialEq, serde::Serialize, serde::Deserialize)]
#[serde(deny_unknown_fields)]
pub struct CampaignReport {
    pub protocol: Protocol,
    pub provenance: Provenance,
    pub snapshot_ids: Vec<String>,
    pub references: Vec<ReferenceRecord>,
    pub candidates: Vec<CandidateReport>,
    pub decisions: Vec<TargetDecision>,
    pub phase_records: Vec<PhaseRecord>,
    pub status: Status,
}
#[derive(Clone, Debug, PartialEq, serde::Serialize, serde::Deserialize)]
pub enum Status {
    Prepared,
    Verified,
    Passed,
    Failed,
    Pending,
    CalibratedSnapshots,
    NotCalibrated,
    RelativeUndefined,
    ReferenceUnconverged,
    ReferenceStagnated,
    ReferencePrecisionLimited,
    BudgetExceeded,
    UnsupportedDomain,
    TargetNotReached,
    ConfirmationFailed,
    MeshBoundary,
    OfflineHelperCompletePhaseProfilingPending,
    FullScopeMeasured,
}
#[derive(Clone, Debug, PartialEq, serde::Serialize, serde::Deserialize)]
pub enum CalibrationFailure {
    InvalidInput(String),
    UnsupportedDomain(String),
    BudgetExceeded(String),
    ReferenceUnconverged(String),
    ReferenceStagnated(String),
    ReferencePrecisionLimited(String),
    RelativeUndefined(String),
    TargetNotReached(
        #[serde(
            serialize_with = "serialize_finite",
            deserialize_with = "deserialize_finite"
        )]
        f64,
    ),
    ConfirmationFailed(
        #[serde(
            serialize_with = "serialize_finite",
            deserialize_with = "deserialize_finite"
        )]
        f64,
    ),
    Io(String),
}
impl std::fmt::Display for CalibrationFailure {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{self:?}")
    }
}
impl std::error::Error for CalibrationFailure {}
impl CalibrationFailure {
    pub fn exit_code(&self) -> i32 {
        match self {
            Self::InvalidInput(_) | Self::UnsupportedDomain(_) | Self::BudgetExceeded(_) => 2,
            Self::ReferenceUnconverged(_)
            | Self::ReferenceStagnated(_)
            | Self::ReferencePrecisionLimited(_)
            | Self::RelativeUndefined(_) => 3,
            Self::TargetNotReached(_) | Self::ConfirmationFailed(_) => 4,
            Self::Io(_) => 1,
        }
    }
}
pub(crate) fn invalid(s: impl Into<String>) -> CalibrationFailure {
    CalibrationFailure::InvalidInput(s.into())
}
pub(crate) fn io(e: impl std::fmt::Display) -> CalibrationFailure {
    CalibrationFailure::Io(e.to_string())
}
/// Check products and caps before any mesh/term allocation or evaluation.
pub fn check_work(
    mesh: [usize; 3],
    work: &WorkCounts,
    heap: u64,
    limits: &Limits,
) -> Result<(), CalibrationFailure> {
    let points = mesh
        .into_iter()
        .try_fold(1usize, |n, x| if x == 0 { None } else { n.checked_mul(x) })
        .ok_or_else(|| CalibrationFailure::BudgetExceeded("mesh product overflow/zero".into()))?;
    if points > limits.max_mesh_points
        || work.real_terms > limits.max_real_terms_per_call
        || work.atom_k_terms > limits.max_atom_k_terms_per_call
        || heap > limits.max_owned_heap_bytes
    {
        return Err(CalibrationFailure::BudgetExceeded(
            "mesh, term or owned heap cap".into(),
        ));
    }
    Ok(())
}
/// Candidate cutoff restriction is deliberately separate from all-image reference R.
pub fn validate_parameters(p: &Parameters, cell: [f64; 3]) -> Result<(), CalibrationFailure> {
    if cell.iter().any(|x| !x.is_finite() || *x <= 0.)
        || [p.kappa, p.real_cutoff, p.kspace_cutoff]
            .iter()
            .any(|x| !x.is_finite() || *x <= 0.)
    {
        return Err(invalid(
            "finite positive candidate parameters and cell required",
        ));
    }
    if p.real_cutoff >= 0.5 * cell.into_iter().fold(f64::INFINITY, f64::min) {
        return Err(invalid("candidate rc must be below half shortest edge"));
    }
    Ok(())
}
pub fn validate_reference_settings(p: &ReferenceSettings) -> Result<(), CalibrationFailure> {
    if [p.kappa, p.real_cutoff, p.kspace_cutoff]
        .iter()
        .any(|x| !x.is_finite() || *x <= 0.)
    {
        Err(invalid("finite positive reference settings required"))
    } else {
        Ok(())
    }
}
#[derive(Clone, Debug, PartialEq, serde::Serialize, serde::Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PreparedManifest {
    pub schema: String,
    pub status: Status,
    pub git_head: String,
    pub protocol_sha256: String,
    pub snapshots_sha256: String,
    pub specification_sha256: String,
    pub fixture_sha256: BTreeMap<String, String>,
    pub source_sha256: BTreeMap<String, String>,
}
pub fn source_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../..")
        .canonicalize()
        .expect("repository root")
}
/// Host prerequisite: direct shasum argv, never a shell and never RTK.
pub fn sha256(path: &Path) -> Result<String, CalibrationFailure> {
    let path = path.canonicalize().map_err(io)?;
    let out = Command::new("shasum")
        .args(["-a", "256"])
        .arg(&path)
        .output()
        .map_err(io)?;
    if !out.status.success() {
        return Err(io("shasum failed"));
    }
    let value = String::from_utf8(out.stdout).map_err(io)?;
    let expected_path = path.to_str().ok_or_else(|| io("non-UTF8 hash path"))?;
    let (hash, name) = value
        .trim_end_matches('\n')
        .split_once("  ")
        .ok_or_else(|| io("invalid single-file shasum output"))?;
    if hash.len() != 64
        || !hash.bytes().all(|x| x.is_ascii_hexdigit())
        || name != expected_path
        || value.lines().count() != 1
    {
        return Err(io("invalid shasum digest/path"));
    }
    Ok(hash.to_ascii_lowercase())
}
fn json_bytes<T: serde::Serialize>(v: &T) -> Result<Vec<u8>, CalibrationFailure> {
    let mut bytes = serde_json::to_vec_pretty(v).map_err(io)?;
    bytes.push(b'\n');
    Ok(bytes)
}
fn write_new(path: &Path, bytes: &[u8]) -> Result<(), CalibrationFailure> {
    let mut f = std::fs::OpenOptions::new()
        .create_new(true)
        .write(true)
        .open(path)
        .map_err(io)?;
    f.write_all(bytes).map_err(io)
}
fn source_hashes() -> Result<BTreeMap<String, String>, CalibrationFailure> {
    fn gather(root: &Path, dir: &Path, paths: &mut Vec<String>) -> Result<(), CalibrationFailure> {
        for entry in std::fs::read_dir(dir).map_err(io)? {
            let path = entry.map_err(io)?.path();
            if path.is_dir() {
                gather(root, &path, paths)?;
            } else if path.extension().is_some_and(|x| x == "rs") {
                paths.push(
                    path.strip_prefix(root)
                        .map_err(io)?
                        .to_str()
                        .ok_or_else(|| io("non-UTF8 source path"))?
                        .into(),
                );
            }
        }
        Ok(())
    }
    let root = source_root();
    let mut paths = vec![
        "Cargo.toml".into(),
        "Cargo.lock".into(),
        "crates/modena-core/Cargo.toml".into(),
        "crates/modena-core/examples/spme_convergence.rs".into(),
    ];
    gather(&root, &root.join("crates/modena-core/src"), &mut paths)?;
    gather(
        &root,
        &root.join("crates/modena-core/examples/spme_calibration"),
        &mut paths,
    )?;
    paths.sort();
    paths
        .into_iter()
        .map(|p| Ok((p.clone(), sha256(&root.join(p))?)))
        .collect()
}

pub fn prepare(out: &Path, protocol: Protocol) -> Result<PreparedManifest, CalibrationFailure> {
    let snapshots = materialize(&protocol)?;
    validate_frozen(&protocol, &snapshots)?;
    // Resolve all external prerequisites before creating the output directory.
    let sources = source_hashes()?;
    let spec = sha256(&source_root().join("docs/scientific-hardening/B05-SPME-PROTOCOL.md"))?;
    let git = Command::new("git")
        .args(["rev-parse", "HEAD"])
        .current_dir(source_root())
        .output()
        .map_err(io)?;
    if !git.status.success() {
        return Err(io("git HEAD unavailable"));
    }
    let head = String::from_utf8(git.stdout)
        .map_err(io)?
        .trim()
        .to_string();
    std::fs::create_dir(out).map_err(io)?;
    write_new(&out.join("protocol.json"), &json_bytes(&protocol)?)?;
    write_new(&out.join("snapshots.json"), &json_bytes(&snapshots)?)?;
    let fixture_dir = out.join("fixtures");
    std::fs::create_dir(&fixture_dir).map_err(io)?;
    let mut hashes = BTreeMap::new();
    for s in &snapshots {
        let path = fixture_dir.join(format!("{}.json", s.id));
        write_new(&path, &json_bytes(s)?)?;
        hashes.insert(s.id.clone(), sha256(&path)?);
    }
    let manifest = PreparedManifest {
        schema: SCHEMA.into(),
        status: Status::Prepared,
        git_head: head,
        protocol_sha256: sha256(&out.join("protocol.json"))?,
        snapshots_sha256: sha256(&out.join("snapshots.json"))?,
        specification_sha256: spec,
        fixture_sha256: hashes,
        source_sha256: sources,
    };
    write_new(&out.join("manifest.json"), &json_bytes(&manifest)?)?;
    Ok(manifest)
}
pub fn validate_prepared(
    dir: &Path,
) -> Result<(Protocol, Vec<Snapshot>, PreparedManifest), CalibrationFailure> {
    let read = |p: &str| std::fs::read(dir.join(p)).map_err(io);
    let m: PreparedManifest = serde_json::from_slice(&read("manifest.json")?).map_err(io)?;
    if m.schema != SCHEMA
        || m.status != Status::Prepared
        || m.protocol_sha256 != sha256(&dir.join("protocol.json"))?
        || m.snapshots_sha256 != sha256(&dir.join("snapshots.json"))?
        || m.source_sha256 != source_hashes()?
        || m.specification_sha256
            != sha256(&source_root().join("docs/scientific-hardening/B05-SPME-PROTOCOL.md"))?
    {
        return Err(invalid(
            "prepared protocol/snapshot/source/specification hash mismatch",
        ));
    }
    let p: Protocol = serde_json::from_slice(&read("protocol.json")?).map_err(io)?;
    let s: Vec<Snapshot> = serde_json::from_slice(&read("snapshots.json")?).map_err(io)?;
    validate_frozen(&p, &s)?;
    if m.fixture_sha256.len() != s.len() {
        return Err(invalid("fixture hash count"));
    }
    for snap in &s {
        let path = dir.join("fixtures").join(format!("{}.json", snap.id));
        if m.fixture_sha256.get(&snap.id) != Some(&sha256(&path)?)
            || std::fs::read(path).map_err(io)? != json_bytes(snap)?
        {
            return Err(invalid(format!("{} fixture hash mismatch", snap.id)));
        }
    }
    Ok((p, s, m))
}

trait Finite {
    fn finite(&self) -> bool;
}
impl Finite for f64 {
    fn finite(&self) -> bool {
        self.is_finite()
    }
}
impl<T: Finite> Finite for Option<T> {
    fn finite(&self) -> bool {
        self.as_ref().is_none_or(Finite::finite)
    }
}
impl<T: Finite> Finite for Vec<T> {
    fn finite(&self) -> bool {
        self.iter().all(Finite::finite)
    }
}
impl<T: Finite, const N: usize> Finite for [T; N] {
    fn finite(&self) -> bool {
        self.iter().all(Finite::finite)
    }
}
fn serialize_finite<T: Finite + serde::Serialize, S: serde::Serializer>(
    value: &T,
    serializer: S,
) -> Result<S::Ok, S::Error> {
    if !value.finite() {
        return Err(serde::ser::Error::custom("nonfinite floating point record"));
    }
    value.serialize(serializer)
}
fn deserialize_finite<'de, T: Finite + serde::Deserialize<'de>, D: serde::Deserializer<'de>>(
    deserializer: D,
) -> Result<T, D::Error> {
    let value = T::deserialize(deserializer)?;
    if !value.finite() {
        return Err(serde::de::Error::custom("nonfinite floating point record"));
    }
    Ok(value)
}

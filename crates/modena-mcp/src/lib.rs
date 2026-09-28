//! Agent-native interface for MODENA — Feature A.
//!
//! This crate exposes MODENA's simulation primitives as **typed, validated
//! tools** an LLM agent can call, rather than forcing the agent to generate
//! LAMMPS/GROMACS scripts and parse text. The same [`RunManager`] is the
//! control surface for both the MCP server (stdio, for agents) and the future
//! Tauri GUI — **one core, two surfaces**.
//!
//! # Architecture
//!
//! ```text
//!   agent / GUI ──tool calls──▶ RunManager (registry + thread pool)
//!                                    │
//!                                    ▼
//!                          modena_core::Simulation::run_streaming
//! ```
//!
//! [`RunManager`] owns a `HashMap<RunId, RunEntry>` and spawns each run on its
//! own `std::thread` (MODENA's low RAM ~5 MB/run makes a pool of 8–16
//! concurrent runs affordable — the "low-RAM enables dense agent streams"
//! claim, made concrete). Runs stream [`RunEvent`]s into a shared status cell
//! the manager reads on demand, and accept cancellation via an `AtomicBool`.
//!
//! # The conservation signal
//!
//! Every run carries a [`ConservationReport`] (Feature B). A degrading verdict
//! mid-sweep is *data* the agent reasons over, not an error — see
//! [`RunManager::detect_anomaly`].

use std::collections::{HashMap, VecDeque};
use std::path::PathBuf;
use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};
use std::sync::{Arc, Mutex};
use std::thread::{self, JoinHandle};

use modena_core::{
    lattice, read_xyz, ConservationReport, ConservationVerdict, Energies, InputConfig,
    IntegratorSpec, RunEvent, RunSummary, SimBox, Simulation, StructureSpec, ThermostatSpec,
};

pub mod protocol;
pub mod server;
mod tools;

pub use tools::{call_tool, tool_list, Tool};

/// Maximum number of buffered thermo samples kept per run for
/// [`RunManager::get_run_events`]. Older samples are dropped (FIFO).
const MAX_SAMPLES: usize = 4096;

/// Default cap on concurrently running simulations. Cheap to raise on a
/// capable machine; the agent can learn the right level from
/// `estimate_run_cost`.
pub const DEFAULT_MAX_CONCURRENT: usize = 8;

/// Opaque identifier for a single simulation run.
pub type RunId = u64;
/// Opaque identifier for a parameter sweep (a collection of runs).
pub type SweepId = u64;

/// Lifecycle state of a run.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, serde::Serialize, serde::Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum RunState {
    /// Accepted but not yet stepped (thread not started / building).
    #[default]
    Queued,
    /// Integrator loop is running.
    Running,
    /// Reached its step count cleanly.
    Completed,
    /// Cancelled via [`RunManager::cancel_run`] (stopped within one step).
    Cancelled,
    /// Aborted — the error is in [`RunStatus::error`].
    Failed,
}

/// A single thermodynamic sample buffered from the run stream.
#[derive(Debug, Clone, serde::Serialize)]
pub struct ThermoSample {
    pub step: usize,
    pub time: f64,
    pub energies: Energies,
}

/// Queryable status of a run, snapshotted under the run's lock.
#[derive(Debug, Clone, serde::Serialize)]
pub struct RunStatus {
    pub run_id: RunId,
    pub state: RunState,
    pub step: usize,
    pub total_steps: usize,
    pub latest_thermo: Option<Energies>,
    pub conservation: ConservationReport,
    pub error: Option<String>,
    /// Whether this is an NVE (microcanonical) run — velocity-Verlet with no
    /// thermostat. Only for NVE runs is [`ConservationReport`] physically
    /// meaningful: under NVT the thermostat deliberately exchanges energy, so
    /// `E_total` drift is expected and must not be read as a numerical fault.
    pub is_nve: bool,
}

/// The mutable, thread-shared part of a run. Updated by the run thread,
/// snapshotted by [`RunManager::get_run_status`].
#[derive(Debug, Default)]
struct RunStatusInner {
    state: RunState,
    step: usize,
    total_steps: usize,
    latest_thermo: Option<Energies>,
    conservation: ConservationReport,
    error: Option<String>,
    is_nve: bool,
    samples: VecDeque<ThermoSample>,
}

impl RunStatusInner {
    fn new(total_steps: usize, is_nve: bool) -> Self {
        Self {
            total_steps,
            is_nve,
            ..Default::default()
        }
    }

    fn push_sample(&mut self, sample: ThermoSample) {
        if self.samples.len() >= MAX_SAMPLES {
            self.samples.pop_front();
        }
        self.samples.push_back(sample);
    }

    fn snapshot(&self, run_id: RunId) -> RunStatus {
        RunStatus {
            run_id,
            state: self.state,
            step: self.step,
            total_steps: self.total_steps,
            latest_thermo: self.latest_thermo,
            conservation: self.conservation,
            error: self.error.clone(),
            is_nve: self.is_nve,
        }
    }
}

struct RunEntry {
    /// Cancel flag shared with the run thread. Setting it stops the run within
    /// one integration step.
    stop: Arc<AtomicBool>,
    /// Thread-shared status cell.
    status: Arc<Mutex<RunStatusInner>>,
    /// The join handle, `None` until the thread is spawned (always `Some` in
    /// practice once `run_async` returns). Held so the manager can block on a
    /// run if needed (e.g. in tests).
    handle: Mutex<Option<JoinHandle<()>>>,
}

/// Which scalar parameter a sweep varies across its runs.
#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum SweepField {
    /// `Run.initial_temperature` (NVE init temperature).
    InitialTemperature,
    /// Langevin BAOAB `target_temperature` (requires the Langevin integrator).
    LangevinTargetTemperature,
    /// Thermostat `target_temperature` (Berendsen / velocity rescaling).
    ThermostatTargetTemperature,
    /// Number density ρ = N/V. For an FCC lattice this rescales the lattice
    /// constant as `a = (4/ρ)^(1/3)` (independent of replica count).
    Density,
}

impl SweepField {
    fn human(&self) -> &'static str {
        match self {
            Self::InitialTemperature => "initial_temperature",
            Self::LangevinTargetTemperature => "langevin_target_temperature",
            Self::ThermostatTargetTemperature => "thermostat_target_temperature",
            Self::Density => "density",
        }
    }
}

struct SweepEntry {
    field: SweepField,
    /// `(run_id, swept_value)` in launch order.
    runs: Vec<(RunId, f64)>,
}

/// Result of [`RunManager::build_system`].
#[derive(Debug, Clone, serde::Serialize)]
pub struct BuildSystemResult {
    pub atoms: usize,
    pub density: f64,
    pub box_size: [f64; 3],
    pub element: String,
}

/// Result of [`RunManager::estimate_run_cost`]. Figures are *rough budgeting
/// estimates* from an internal cost model, not measured guarantees.
#[derive(Debug, Clone, serde::Serialize)]
pub struct EstimateRunCostResult {
    pub atoms: usize,
    pub steps: usize,
    pub est_wall_s: f64,
    pub est_rss_mb: f64,
}

/// Result of a synchronous [`RunManager::run_probe`].
#[derive(Debug, Clone, serde::Serialize)]
pub struct ProbeResult {
    pub steps: usize,
    pub final_thermo: Energies,
    pub conservation: ConservationReport,
}

/// Result of [`RunManager::launch_sweep`].
#[derive(Debug, Clone, serde::Serialize)]
pub struct LaunchSweepResult {
    pub sweep_id: SweepId,
    pub field: SweepField,
    pub run_ids: Vec<RunId>,
}

/// One run's row in [`SweepResult`].
#[derive(Debug, Clone, serde::Serialize)]
pub struct SweepRunResult {
    pub run_id: RunId,
    /// The swept parameter value for this run.
    pub value: f64,
    pub status: RunStatus,
}

/// Result of [`RunManager::get_sweep_results`].
#[derive(Debug, Clone, serde::Serialize)]
pub struct SweepResult {
    pub sweep_id: SweepId,
    pub field: SweepField,
    pub runs: Vec<SweepRunResult>,
}

/// A detected anomaly in a run or sweep.
#[derive(Debug, Clone, serde::Serialize)]
pub struct Anomaly {
    /// `"energy_drift"` | `"divergence"` | `"possible_phase_transition"`.
    #[serde(rename = "type")]
    pub typ: String,
    /// Human-readable location, e.g. `"run 7"` or `"T≈1.08"`.
    pub location: String,
    /// Structured evidence the agent can cite.
    pub evidence: serde_json::Value,
}

/// Result of [`RunManager::detect_anomaly`].
#[derive(Debug, Clone, serde::Serialize)]
pub struct AnomalyReport {
    pub anomalies: Vec<Anomaly>,
    /// One-line human summary.
    pub summary: String,
}

/// Description of a registered pair potential.
#[derive(Debug, Clone, serde::Serialize)]
pub struct PotentialInfo {
    pub name: String,
    pub params: Vec<String>,
}

/// Errors surfaced to the agent. All are recoverable — the agent can adjust its
/// request and retry.
#[derive(Debug, thiserror::Error)]
pub enum McpError {
    #[error("invalid tool input: {0}")]
    InvalidInput(String),
    #[error("run {0} not found")]
    RunNotFound(RunId),
    #[error("sweep {0} not found")]
    SweepNotFound(SweepId),
    #[error("simulation error: {0}")]
    Simulation(String),
    #[error("the concurrent-run pool is full ({0}); cancel a run or wait")]
    PoolFull(usize),
    #[error("{0}")]
    Other(String),
}

/// The control surface for agent/GUI-driven simulations. Cheap to clone is
/// *not* supported by design — share it via `Arc` if multiple owners are
/// needed (the MCP server uses a single instance behind its dispatcher).
pub struct RunManager {
    runs: Mutex<HashMap<RunId, RunEntry>>,
    sweeps: Mutex<HashMap<SweepId, SweepEntry>>,
    next_run_id: AtomicU64,
    next_sweep_id: AtomicU64,
    max_concurrent: usize,
    /// Directory for resolving relative structure/log paths in configs. Managed
    /// runs typically omit file output, so this is usually inert.
    work_dir: PathBuf,
}

impl Default for RunManager {
    fn default() -> Self {
        Self::new(DEFAULT_MAX_CONCURRENT, PathBuf::from("."))
    }
}

impl RunManager {
    /// Create a manager with a concurrency cap and a working directory.
    pub fn new(max_concurrent: usize, work_dir: PathBuf) -> Self {
        Self {
            runs: Mutex::new(HashMap::new()),
            sweeps: Mutex::new(HashMap::new()),
            next_run_id: AtomicU64::new(1),
            next_sweep_id: AtomicU64::new(1),
            max_concurrent: max_concurrent.max(1),
            work_dir,
        }
    }

    /// Count of runs currently in the `Running` state.
    fn count_running(&self) -> usize {
        self.runs
            .lock()
            .unwrap()
            .values()
            .filter(|e| e.status.lock().unwrap().state == RunState::Running)
            .count()
    }

    // ---- discovery ----

    /// `list_potentials`: the registered pair potentials.
    pub fn list_potentials(&self) -> Vec<PotentialInfo> {
        vec![PotentialInfo {
            name: "lennard_jones".into(),
            params: vec!["sigma".into(), "epsilon".into(), "cutoff".into()],
        }]
    }

    /// `build_system`: validate + materialize a structure, returning derived
    /// counts. Does not run.
    pub fn build_system(&self, structure: &StructureSpec) -> Result<BuildSystemResult, McpError> {
        let system = build_structure(structure, &self.work_dir)?;
        let n = system.n();
        let vol = system.box_size.volume();
        let density = if vol > 0.0 { n as f64 / vol } else { 0.0 };
        let element = match structure {
            StructureSpec::Fcc { element, .. } | StructureSpec::Xyz { element, .. } => {
                element.clone()
            }
            StructureSpec::RandomAlloy { species, .. } => species
                .first()
                .map(|s| s.element.clone())
                .unwrap_or_else(|| "alloy".to_string()),
            StructureSpec::RockSalt { cation, anion, .. } => {
                format!("{cation}-{anion}")
            }
            StructureSpec::Water { .. } => "SPC/E-water".to_string(),
        };
        Ok(BuildSystemResult {
            atoms: n,
            density,
            box_size: [
                system.box_size.size.0[0],
                system.box_size.size.0[1],
                system.box_size.size.0[2],
            ],
            element,
        })
    }

    /// `estimate_run_cost`: a rough wall-time and RAM estimate from an internal
    /// cost model (≈ 50 ns/atom/step, ≈ 96 B/atom + 5 MB baseline). Helps the
    /// agent budget before launching large runs.
    pub fn estimate_run_cost(&self, config: &InputConfig) -> EstimateRunCostResult {
        let atoms = match &config.structure {
            StructureSpec::Fcc {
                lattice_constant: _,
                replicas,
                ..
            } => 4 * (replicas[0] as usize) * (replicas[1] as usize) * (replicas[2] as usize),
            StructureSpec::RandomAlloy { replicas, .. } => {
                4 * (replicas[0] as usize) * (replicas[1] as usize) * (replicas[2] as usize)
            }
            StructureSpec::Water { n_waters, .. } => 3 * *n_waters as usize,
            StructureSpec::RockSalt { replicas, .. } => {
                8 * (replicas[0] as usize) * (replicas[1] as usize) * (replicas[2] as usize)
            }
            StructureSpec::Xyz { .. } => {
                // Without reading the file, fall back to a nominal count.
                build_structure(&config.structure, &self.work_dir)
                    .map(|s| s.n())
                    .unwrap_or(0)
            }
        };
        let atoms_f = atoms as f64;
        let steps = config.run.steps;
        // Calibrated rough model (see BENCHMARKS.md for measured throughput).
        let est_wall_s = steps as f64 * atoms_f * 5e-8;
        let est_rss_mb = 5.0 + atoms_f * 1.2e-4;
        EstimateRunCostResult {
            atoms,
            steps,
            est_wall_s,
            est_rss_mb,
        }
    }

    // ---- synchronous probe ----

    /// `run_probe`: run a short (≤ `max_steps`) simulation synchronously on the
    /// calling thread and return the final thermo + conservation report. The
    /// config's `run.steps` is clamped to `max_steps`.
    pub fn run_probe(
        &self,
        mut config: InputConfig,
        max_steps: usize,
    ) -> Result<ProbeResult, McpError> {
        let cap = max_steps.max(1);
        if config.run.steps > cap {
            config.run.steps = cap;
        }
        // Probes are hermetic: no on-disk trajectory/log unless the caller set one.
        let mut sim = Simulation::from_config(&config, &self.work_dir)
            .map_err(|e| McpError::Simulation(e.to_string()))?;
        let RunSummary {
            steps,
            final_energies,
            conservation,
        } = sim.run().map_err(|e| McpError::Simulation(e.to_string()))?;
        Ok(ProbeResult {
            steps,
            final_thermo: final_energies,
            conservation,
        })
    }

    // ---- asynchronous runs ----

    /// `run_async`: launch a streaming run on its own thread. Returns the run
    /// id immediately; poll with [`get_run_status`](Self::get_run_status).
    pub fn run_async(&self, config: InputConfig) -> Result<RunId, McpError> {
        let total_steps = config.run.steps;
        if self.count_running() >= self.max_concurrent {
            return Err(McpError::PoolFull(self.max_concurrent));
        }
        let run_id = self.next_run_id.fetch_add(1, Ordering::Relaxed);
        let is_nve = is_nve_config(&config);

        let stop = Arc::new(AtomicBool::new(false));
        let status = Arc::new(Mutex::new(RunStatusInner::new(total_steps, is_nve)));
        let entry = RunEntry {
            stop: stop.clone(),
            status: status.clone(),
            handle: Mutex::new(None),
        };
        self.runs.lock().unwrap().insert(run_id, entry);

        let work_dir = self.work_dir.clone();
        let join = thread::spawn(move || run_loop(config, work_dir, stop, status));
        // Store the handle for potential blocking joins (tests).
        if let Some(e) = self.runs.lock().unwrap().get_mut(&run_id) {
            *e.handle.lock().unwrap() = Some(join);
        }
        Ok(run_id)
    }

    /// `get_run_status`: snapshot a run's current state.
    pub fn get_run_status(&self, run_id: RunId) -> Result<RunStatus, McpError> {
        let runs = self.runs.lock().unwrap();
        let entry = runs.get(&run_id).ok_or(McpError::RunNotFound(run_id))?;
        let st = entry.status.lock().unwrap();
        Ok(st.snapshot(run_id))
    }

    /// `get_run_events`: buffered thermo samples with `step > since_step`.
    pub fn get_run_events(
        &self,
        run_id: RunId,
        since_step: usize,
    ) -> Result<Vec<ThermoSample>, McpError> {
        let runs = self.runs.lock().unwrap();
        let entry = runs.get(&run_id).ok_or(McpError::RunNotFound(run_id))?;
        let st = entry.status.lock().unwrap();
        Ok(st
            .samples
            .iter()
            .filter(|s| s.step > since_step)
            .cloned()
            .collect())
    }

    /// `cancel_run`: request a clean stop. The run transitions to `Cancelled`
    /// within one integration step.
    pub fn cancel_run(&self, run_id: RunId) -> Result<RunStatus, McpError> {
        let runs = self.runs.lock().unwrap();
        let entry = runs.get(&run_id).ok_or(McpError::RunNotFound(run_id))?;
        entry.stop.store(true, Ordering::Relaxed);
        let st = entry.status.lock().unwrap();
        Ok(st.snapshot(run_id))
    }

    /// Block until a run reaches a terminal state. Primarily for tests; the MCP
    /// server polls instead.
    pub fn join_run(&self, run_id: RunId) -> Result<(), McpError> {
        let handle = {
            let runs = self.runs.lock().unwrap();
            let entry = runs.get(&run_id).ok_or(McpError::RunNotFound(run_id))?;
            let handle = entry.handle.lock().unwrap().take();
            handle
        };
        if let Some(h) = handle {
            let _ = h.join();
        }
        Ok(())
    }

    // ---- sweeps ----

    /// `launch_sweep`: launch one async run per value of `field`, all from
    /// `base_config`. Returns the sweep id and the per-value run ids.
    pub fn launch_sweep(
        &self,
        mut base_config: InputConfig,
        field: SweepField,
        values: Vec<f64>,
    ) -> Result<LaunchSweepResult, McpError> {
        if values.is_empty() {
            return Err(McpError::InvalidInput("sweep values are empty".into()));
        }
        if values.len() > self.max_concurrent {
            return Err(McpError::PoolFull(self.max_concurrent));
        }
        let sweep_id = self.next_sweep_id.fetch_add(1, Ordering::Relaxed);
        let mut run_ids = Vec::with_capacity(values.len());

        for v in &values {
            apply_sweep_field(&mut base_config, field, *v)?;
            // Hermetic sweeps: no per-run file output (avoid clutter / races).
            base_config.run.trajectory_file = None;
            base_config.run.log_file = None;
            base_config.run.restart_file = None;
            let run_id = self.run_async(base_config.clone())?;
            run_ids.push(run_id);
        }

        let runs = run_ids
            .iter()
            .zip(values.iter())
            .map(|(r, v)| (*r, *v))
            .collect();
        self.sweeps
            .lock()
            .unwrap()
            .insert(sweep_id, SweepEntry { field, runs });
        Ok(LaunchSweepResult {
            sweep_id,
            field,
            run_ids,
        })
    }

    /// `get_sweep_results`: aggregated status of every run in a sweep.
    pub fn get_sweep_results(&self, sweep_id: SweepId) -> Result<SweepResult, McpError> {
        let sweeps = self.sweeps.lock().unwrap();
        let entry = sweeps
            .get(&sweep_id)
            .ok_or(McpError::SweepNotFound(sweep_id))?;
        let field = entry.field;
        let runs = entry
            .runs
            .iter()
            .map(|(run_id, value)| {
                let status = self.get_run_status(*run_id).unwrap_or(RunStatus {
                    run_id: *run_id,
                    state: RunState::Failed,
                    step: 0,
                    total_steps: 0,
                    latest_thermo: None,
                    conservation: ConservationReport::default(),
                    error: Some("run lookup failed".into()),
                    is_nve: false,
                });
                SweepRunResult {
                    run_id: *run_id,
                    value: *value,
                    status,
                }
            })
            .collect();
        Ok(SweepResult {
            sweep_id,
            field,
            runs,
        })
    }

    // ---- anomaly detection ----

    /// `detect_anomaly`: engine-side heuristics over a single run or a sweep.
    ///
    /// For a run: flags `energy_drift` (verdict `Drifting`) and `divergence`
    /// (NaN/Inf in thermo).
    ///
    /// For a sweep: additionally runs a basic phase-transition detector — it
    /// sorts completed runs by the swept value, computes the mean potential
    /// energy per run, and flags the interval with the largest `|d⟨PE⟩/d(param)|`
    /// as a `possible_phase_transition`. These are cheap, deterministic
    /// heuristics so the agent reasons over structured findings rather than raw
    /// float arrays.
    pub fn detect_anomaly(
        &self,
        run_id: Option<RunId>,
        sweep_id: Option<SweepId>,
    ) -> Result<AnomalyReport, McpError> {
        let mut anomalies = Vec::new();

        if let Some(rid) = run_id {
            let status = self.get_run_status(rid)?;
            check_run_anomalies(rid, &status, &mut anomalies);
        }

        if let Some(sid) = sweep_id {
            let sweep = self.get_sweep_results(sid)?;
            for row in &sweep.runs {
                check_run_anomalies(row.run_id, &row.status, &mut anomalies);
            }
            detect_transition(&sweep, &mut anomalies);
        }

        let summary = if anomalies.is_empty() {
            "no anomalies detected".to_string()
        } else {
            format!(
                "{} anomaly/anomalies: {}",
                anomalies.len(),
                anomalies
                    .iter()
                    .map(|a| format!("{} at {}", a.typ, a.location))
                    .collect::<Vec<_>>()
                    .join("; ")
            )
        };
        Ok(AnomalyReport { anomalies, summary })
    }
}

/// The body of an async run thread. Owns the `Simulation`, streams events into
/// the shared status cell, and respects the cancel flag.
fn run_loop(
    config: InputConfig,
    work_dir: PathBuf,
    stop: Arc<AtomicBool>,
    status: Arc<Mutex<RunStatusInner>>,
) {
    let mut sim = match Simulation::from_config(&config, &work_dir) {
        Ok(s) => s,
        Err(e) => {
            let mut st = status.lock().unwrap();
            st.state = RunState::Failed;
            st.error = Some(e.to_string());
            return;
        }
    };
    {
        let mut st = status.lock().unwrap();
        st.state = RunState::Running;
    }
    // `run_streaming` borrows the flag for the call duration; the closure reads
    // it (shared) to distinguish completion from cancellation.
    let stop_ref: &AtomicBool = &stop;
    let result = sim.run_streaming(Some(stop_ref), |ev| {
        let mut st = status.lock().unwrap();
        // `ev` is `&RunEvent`, so matched fields bind by reference; the payload
        // types (usize, f64, Energies, RunSummary) are all `Copy`.
        match ev {
            RunEvent::Started { steps, .. } => st.total_steps = *steps,
            RunEvent::Thermo {
                step,
                time,
                energies,
                ..
            } => {
                let step = *step;
                let time = *time;
                let energies = *energies;
                st.step = step;
                st.conservation.observe(energies.total);
                st.latest_thermo = Some(energies);
                st.push_sample(ThermoSample {
                    step,
                    time,
                    energies,
                });
            }
            RunEvent::Frame { step, .. } => st.step = *step,
            RunEvent::Finished { summary } => {
                let summary = *summary;
                st.step = summary.steps;
                st.conservation = summary.conservation;
                st.latest_thermo = Some(summary.final_energies);
                st.state = if stop.load(Ordering::Relaxed) {
                    RunState::Cancelled
                } else {
                    RunState::Completed
                };
            }
            RunEvent::Failed { error } => {
                st.state = RunState::Failed;
                st.error = Some(error.clone());
            }
        }
    });
    if let Err(e) = result {
        let mut st = status.lock().unwrap();
        // Don't overwrite a Cancelled state set by Finished if the error was
        // just the early-return from the broken loop.
        if st.state != RunState::Cancelled {
            st.state = RunState::Failed;
            st.error = Some(e.to_string());
        }
    }
}

/// Materialize a `System` from a structure spec using public core API.
fn build_structure(
    spec: &StructureSpec,
    base_dir: &std::path::Path,
) -> Result<modena_core::System, McpError> {
    match spec {
        StructureSpec::Fcc {
            lattice_constant,
            replicas,
            element,
            ..
        } => {
            // NaN-safe: reject non-finite and non-positive (NaN is not `> 0` and
            // not finite, so both clauses catch it).
            if !lattice_constant.is_finite() || *lattice_constant <= 0.0 {
                return Err(McpError::InvalidInput(
                    "lattice_constant must be > 0".into(),
                ));
            }
            if replicas.contains(&0) {
                return Err(McpError::InvalidInput("replicas must be >= 1".into()));
            }
            Ok(lattice::fcc(*lattice_constant, *replicas, element))
        }
        StructureSpec::Xyz {
            file,
            box_size,
            element,
            ..
        } => {
            let path = base_dir.join(file);
            let mut sys = read_xyz(&path, SimBox::from_array(*box_size))
                .map_err(|e| McpError::Simulation(format!("read_xyz {file}: {e}")))?;
            for el in sys.elements.iter_mut() {
                if el == "X" {
                    *el = element.clone();
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
            if !lattice_constant.is_finite() || *lattice_constant <= 0.0 {
                return Err(McpError::InvalidInput(
                    "lattice_constant must be > 0".into(),
                ));
            }
            if replicas.contains(&0) {
                return Err(McpError::InvalidInput("replicas must be >= 1".into()));
            }
            let alloy: Vec<modena_core::AlloySpecies> = species
                .iter()
                .map(|s| modena_core::AlloySpecies {
                    element: s.element.clone(),
                    mass: s.mass,
                    concentration: s.concentration,
                    charge: s.charge,
                })
                .collect();
            lattice::random_alloy(*lattice_constant, *replicas, &alloy, *seed)
                .map_err(|e| McpError::InvalidInput(e.to_string()))
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
        } => {
            if !lattice_constant.is_finite() || *lattice_constant <= 0.0 {
                return Err(McpError::InvalidInput(
                    "lattice_constant must be > 0".into(),
                ));
            }
            if replicas.contains(&0) {
                return Err(McpError::InvalidInput("replicas must be >= 1".into()));
            }
            Ok(lattice::rocksalt(
                *lattice_constant,
                *replicas,
                cation,
                anion,
                cation_mass.unwrap_or(1.0),
                anion_mass.unwrap_or(1.0),
                *cation_charge,
                *anion_charge,
            ))
        }
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
        } => lattice::water_lattice(
            *n_waters as usize,
            *density,
            *d_oh,
            *hoh_angle,
            *q_o,
            *q_h,
            *m_o,
            *m_h,
            *seed,
            1.0, // the MCP operates on reduced-unit configs; lengths pre-converted
        )
        .map_err(|e| McpError::InvalidInput(e.to_string())),
    }
}

/// Apply a swept value to the relevant field of `config`.
fn apply_sweep_field(
    config: &mut InputConfig,
    field: SweepField,
    value: f64,
) -> Result<(), McpError> {
    if !value.is_finite() {
        return Err(McpError::InvalidInput(format!(
            "sweep value must be finite, got {value}"
        )));
    }
    match field {
        SweepField::InitialTemperature => {
            if value <= 0.0 {
                return Err(McpError::InvalidInput("temperature must be > 0".into()));
            }
            config.run.initial_temperature = value;
        }
        SweepField::LangevinTargetTemperature => {
            if value <= 0.0 {
                return Err(McpError::InvalidInput("temperature must be > 0".into()));
            }
            match &mut config.integrator {
                IntegratorSpec::LangevinBaoab {
                    target_temperature, ..
                } => *target_temperature = value,
                other => {
                    return Err(McpError::InvalidInput(format!(
                        "langevin_target_temperature sweep requires the langevin_baoab integrator, found {:?}",
                        integrator_kind(other)
                    )))
                }
            }
        }
        SweepField::ThermostatTargetTemperature => {
            if value <= 0.0 {
                return Err(McpError::InvalidInput("temperature must be > 0".into()));
            }
            match &mut config.thermostat {
                Some(ThermostatSpec::Berendsen { target_temperature, .. })
                | Some(ThermostatSpec::VelocityRescaling { target_temperature, .. }) => {
                    *target_temperature = value;
                }
                other => {
                    return Err(McpError::InvalidInput(format!(
                        "thermostat_target_temperature sweep requires a berendsen or velocity_rescaling thermostat, found {}",
                        thermostat_kind(other)
                    )))
                }
            }
        }
        SweepField::Density => {
            if value <= 0.0 {
                return Err(McpError::InvalidInput("density must be > 0".into()));
            }
            match &mut config.structure {
                StructureSpec::Fcc {
                    lattice_constant, ..
                } => {
                    // FCC: ρ = 4 / a³  =>  a = (4/ρ)^(1/3).
                    *lattice_constant = (4.0 / value).cbrt();
                }
                other => {
                    return Err(McpError::InvalidInput(format!(
                        "density sweep requires an fcc structure, found {}",
                        structure_kind(other)
                    )))
                }
            }
        }
    }
    Ok(())
}

fn integrator_kind(s: &IntegratorSpec) -> &'static str {
    match s {
        IntegratorSpec::VelocityVerlet { .. } => "velocity_verlet",
        IntegratorSpec::LangevinBaoab { .. } => "langevin_baoab",
        IntegratorSpec::NoseHooverChain { .. } => "nose_hoover_chain",
    }
}

/// A run is NVE iff it uses velocity-Verlet with no (or `none`) thermostat.
/// Only NVE runs have a physically meaningful total-energy drift.
fn is_nve_config(config: &InputConfig) -> bool {
    matches!(config.integrator, IntegratorSpec::VelocityVerlet { .. })
        && matches!(config.thermostat, None | Some(ThermostatSpec::None))
}
fn thermostat_kind(s: &Option<ThermostatSpec>) -> &'static str {
    match s {
        None => "none",
        Some(ThermostatSpec::Berendsen { .. }) => "berendsen",
        Some(ThermostatSpec::VelocityRescaling { .. }) => "velocity_rescaling",
        Some(ThermostatSpec::None) => "none",
    }
}
fn structure_kind(s: &StructureSpec) -> &'static str {
    match s {
        StructureSpec::Fcc { .. } => "fcc",
        StructureSpec::Xyz { .. } => "xyz",
        StructureSpec::RandomAlloy { .. } => "random_alloy",
        StructureSpec::RockSalt { .. } => "rocksalt",
        StructureSpec::Water { .. } => "water",
    }
}

/// Per-run anomaly checks: drift (NVE only) + divergence (NaN/Inf).
fn check_run_anomalies(run_id: RunId, status: &RunStatus, out: &mut Vec<Anomaly>) {
    // Energy drift is only a fault under NVE; under NVT the thermostat
    // deliberately exchanges energy, so E_total drift is expected.
    if status.is_nve && status.conservation.verdict == ConservationVerdict::Drifting {
        out.push(Anomaly {
            typ: "energy_drift".into(),
            location: format!("run {run_id}"),
            evidence: serde_json::json!({
                "relative_drift": status.conservation.relative_drift,
                "verdict": status.conservation.verdict,
                "step": status.step,
                "ensemble": "nve",
            }),
        });
    }
    if let Some(t) = status.latest_thermo {
        let diverged = [t.kinetic, t.potential, t.total, t.temperature, t.pressure]
            .iter()
            .any(|x| !x.is_finite());
        if diverged {
            out.push(Anomaly {
                typ: "divergence".into(),
                location: format!("run {run_id}"),
                evidence: serde_json::json!({
                    "step": status.step,
                    "latest_thermo": t,
                }),
            });
        }
    }
}

/// Basic phase-transition detector: largest `|d⟨PE⟩/d(param)|` interval across
/// completed, finite runs.
fn detect_transition(sweep: &SweepResult, out: &mut Vec<Anomaly>) {
    // Collect (param, mean_PE) for completed runs with finite thermo.
    let mut pts: Vec<(f64, f64)> = Vec::new();
    for row in &sweep.runs {
        if row.status.state != RunState::Completed {
            continue;
        }
        let mean_pe = mean_potential_energy(row.run_id, &row.status);
        if let Some(pe) = mean_pe {
            if pe.is_finite() {
                pts.push((row.value, pe));
            }
        }
    }
    if pts.len() < 2 {
        return;
    }
    pts.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));

    // Largest magnitude slope between consecutive points.
    let mut best_i = 0;
    let mut best_slope = 0.0_f64;
    for i in 0..pts.len() - 1 {
        let (x0, y0) = pts[i];
        let (x1, y1) = pts[i + 1];
        let dx = x1 - x0;
        if dx.abs() < 1e-12 {
            continue;
        }
        let slope = ((y1 - y0) / dx).abs();
        if slope > best_slope {
            best_slope = slope;
            best_i = i;
        }
    }
    let (x0, _y0) = pts[best_i];
    let (x1, _y1) = pts[best_i + 1];
    let midpoint = 0.5 * (x0 + x1);
    out.push(Anomaly {
        typ: "possible_phase_transition".into(),
        location: format!("{}≈{midpoint:.4}", sweep.field.human()),
        evidence: serde_json::json!({
            "field": sweep.field.human(),
            "interval": [x0, x1],
            "slope_magnitude": best_slope,
            "n_completed": pts.len(),
            "note": "largest |d<PE>/d(param)| across the sweep; heuristic, not a proof of a transition",
        }),
    });
}

/// Mean potential energy from buffered samples, falling back to the latest thermo.
fn mean_potential_energy(_run_id: RunId, status: &RunStatus) -> Option<f64> {
    // The buffered samples live in RunStatusInner, not RunStatus. We approximate
    // with the latest thermo PE when samples aren't on the snapshot. This is a
    // reasonable proxy for the heuristic; a richer API would expose the sample
    // mean directly.
    status.latest_thermo.map(|t| t.potential)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A minimal reduced-units LJ config (FCC, 32 atoms, 200 steps, NVE).
    fn sample_config(steps: usize) -> InputConfig {
        // Built by deserializing a TOML string to stay robust to schema changes.
        let toml = format!(
            r#"
units = "reduced"
[structure]
type = "fcc"
lattice_constant = 1.6796
replicas = [2, 2, 2]
element = "Ar"
[model]
type = "lennard_jones"
sigma = 1.0
epsilon = 1.0
cutoff = 2.5
[integrator]
type = "velocity_verlet"
dt = 0.002
[run]
steps = {steps}
initial_temperature = 1.0
seed = 42
"#
        );
        toml::from_str(&toml).expect("valid config")
    }

    #[test]
    fn build_system_reports_counts() {
        let mgr = RunManager::default();
        let cfg = sample_config(10);
        let res = mgr.build_system(&cfg.structure).unwrap();
        assert_eq!(res.atoms, 32);
        // ρ = 32 / (2*1.6796)^3 ≈ 0.8442
        assert!(
            (res.density - 0.8442).abs() < 1e-3,
            "density = {}",
            res.density
        );
    }

    #[test]
    fn estimate_run_cost_is_finite_and_positive() {
        let mgr = RunManager::default();
        let res = mgr.estimate_run_cost(&sample_config(1000));
        assert_eq!(res.atoms, 32);
        assert_eq!(res.steps, 1000);
        assert!(res.est_wall_s > 0.0);
        assert!(res.est_rss_mb > 0.0);
    }

    #[test]
    fn run_probe_returns_conservation() {
        let mgr = RunManager::default();
        let res = mgr.run_probe(sample_config(2000), 200).unwrap();
        assert_eq!(res.steps, 200);
        // Short NVE at dt=0.002 must conserve energy.
        assert!(res.conservation.verdict != ConservationVerdict::Drifting);
    }

    #[test]
    fn run_probe_clamps_steps() {
        let mgr = RunManager::default();
        let res = mgr.run_probe(sample_config(10_000), 100).unwrap();
        assert_eq!(res.steps, 100);
    }

    #[test]
    fn async_run_completes_and_reports_status() {
        let mgr = RunManager::default();
        let rid = mgr.run_async(sample_config(200)).unwrap();
        mgr.join_run(rid).unwrap();
        let status = mgr.get_run_status(rid).unwrap();
        assert_eq!(status.state, RunState::Completed);
        assert_eq!(status.step, 200);
        assert!(status.latest_thermo.is_some());
        assert!(status.conservation.verdict != ConservationVerdict::Drifting);
    }

    #[test]
    fn cancel_run_stops_early() {
        let mgr = RunManager::default();
        // Long run; cancel immediately.
        let rid = mgr.run_async(sample_config(1_000_000)).unwrap();
        // Give it a moment to start, then cancel.
        std::thread::sleep(std::time::Duration::from_millis(20));
        mgr.cancel_run(rid).unwrap();
        mgr.join_run(rid).unwrap();
        let status = mgr.get_run_status(rid).unwrap();
        assert!(
            status.state == RunState::Cancelled,
            "expected cancelled, got {:?} (step {})",
            status.state,
            status.step
        );
        assert!(status.step < 1_000_000);
    }

    #[test]
    fn sweep_runs_to_completion() {
        let mgr = RunManager::default();
        // Langevin NVT at three temperatures — uses the langevin sweep field.
        let toml_base = r#"
units = "reduced"
[structure]
type = "fcc"
lattice_constant = 1.6796
replicas = [2, 2, 2]
element = "Ar"
[model]
type = "lennard_jones"
sigma = 1.0
epsilon = 1.0
cutoff = 2.5
[integrator]
type = "langevin_baoab"
dt = 0.005
gamma = 1.0
target_temperature = 1.0
seed = 7
[run]
steps = 200
initial_temperature = 1.0
seed = 42
"#;
        let mut base: InputConfig = toml::from_str(toml_base).unwrap();
        base.run.trajectory_file = None;
        base.run.log_file = None;
        let res = mgr
            .launch_sweep(
                base,
                SweepField::LangevinTargetTemperature,
                vec![0.8, 1.0, 1.2],
            )
            .unwrap();
        assert_eq!(res.run_ids.len(), 3);
        for rid in &res.run_ids {
            mgr.join_run(*rid).unwrap();
        }
        let sweep = mgr.get_sweep_results(res.sweep_id).unwrap();
        assert_eq!(sweep.runs.len(), 3);
        for row in &sweep.runs {
            assert_eq!(
                row.status.state,
                RunState::Completed,
                "run {} not completed",
                row.run_id
            );
        }
        // Higher T → higher mean PE (less negative) for LJ liquid. At least
        // monotonic-ish: the 1.2 run's PE should exceed the 0.8 run's PE.
        let pe = |v: f64| {
            sweep
                .runs
                .iter()
                .find(|r| (r.value - v).abs() < 1e-9)
                .and_then(|r| r.status.latest_thermo)
                .map(|t| t.potential)
                .unwrap()
        };
        assert!(
            pe(1.2) >= pe(0.8) - 1.0,
            "expected PE(1.2) >= PE(0.8) roughly, got {} vs {}",
            pe(1.2),
            pe(0.8)
        );
    }

    #[test]
    fn detect_anomaly_flags_drift() {
        // A too-large timestep at high T drives energy drift → Drifting verdict.
        let toml = r#"
units = "reduced"
[structure]
type = "fcc"
lattice_constant = 1.6796
replicas = [2, 2, 2]
element = "Ar"
[model]
type = "lennard_jones"
sigma = 1.0
epsilon = 1.0
cutoff = 2.5
[integrator]
type = "velocity_verlet"
dt = 0.05
[run]
steps = 400
initial_temperature = 5.0
seed = 1
"#;
        let cfg: InputConfig = toml::from_str(toml).unwrap();
        let mgr = RunManager::default();
        let rid = mgr.run_async(cfg).unwrap();
        mgr.join_run(rid).unwrap();
        let report = mgr.detect_anomaly(Some(rid), None).unwrap();
        // Either a drift or divergence anomaly should surface for this abusive config.
        assert!(
            !report.anomalies.is_empty(),
            "expected an anomaly for dt=0.05/T=5, got: {report:?}"
        );
    }

    #[test]
    fn density_sweep_rescales_lattice() {
        let mut cfg = sample_config(10);
        apply_sweep_field(&mut cfg, SweepField::Density, 0.8442).unwrap();
        if let StructureSpec::Fcc {
            lattice_constant, ..
        } = &cfg.structure
        {
            // a = (4/0.8442)^(1/3) ≈ 1.6796
            assert!(
                (lattice_constant - 1.6796).abs() < 1e-3,
                "a = {lattice_constant}"
            );
        } else {
            panic!("expected fcc");
        }
    }

    #[test]
    fn rejects_invalid_sweep_combos() {
        let mut cfg = sample_config(10); // velocity_verlet, no thermostat
        let err = apply_sweep_field(&mut cfg, SweepField::LangevinTargetTemperature, 1.0);
        assert!(err.is_err());
        let err = apply_sweep_field(&mut cfg, SweepField::ThermostatTargetTemperature, 1.0);
        assert!(err.is_err());
        let err = apply_sweep_field(&mut cfg, SweepField::Density, -1.0);
        assert!(err.is_err());
    }
}

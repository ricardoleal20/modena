//! Pinned ASE EMT provider. No Python dependency is needed to build.
use std::io::{Read, Seek, SeekFrom, Write};
use std::process::{Command, Stdio};
use std::sync::Mutex;
use std::time::{Duration, Instant};

use crate::{ForceModel, ModenaError, NeighborList, Result, System, UnitFactors, Virial};
use serde::{Deserialize, Serialize};

const PROTOCOL: &str = "modena-ase-emt-v1";
const WORKER: &str = include_str!("external/ase_emt_worker.py");
const MAX_STDOUT: u64 = 16 * 1024 * 1024;
const MAX_STDERR: u64 = 64 * 1024;
mod session;
/// Exact SI conversion from eV per particle to kcal per mole.
pub const EV_TO_KCAL_MOL: f64 = 1.602176634e-19 * 6.02214076e23 / 4184.0;

fn default_python() -> String {
    "python3".into()
}
fn full_pbc() -> [bool; 3] {
    [true; 3]
}
fn default_timeout() -> u64 {
    30_000
}

/// Process transport does not contribute to the physical provider identity.
#[derive(Debug, Clone, Copy, Default, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum AseEmtTransport {
    #[default]
    OneShot,
    Persistent,
}

/// Fixed EMT parameters and required physical/model pins.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct AseEmtSpec {
    #[serde(default = "default_python")]
    pub python: String,
    pub ase_version: String,
    pub model_sha256: String,
    pub atomic_numbers: Vec<u8>,
    pub composition: String,
    #[serde(default = "full_pbc")]
    pub pbc: [bool; 3],
    #[serde(default)]
    pub asap_cutoff: bool,
    #[serde(default = "default_timeout")]
    pub timeout_ms: u64,
    #[serde(default)]
    pub transport: AseEmtTransport,
}

/// Explicit wire unit convention.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct ProviderUnits {
    pub positions: String,
    pub energy: String,
    pub forces: String,
    pub stress: String,
}

/// Physical calculator identity; the executable path is transport only.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct ProviderIdentity {
    pub protocol: String,
    pub provider: String,
    pub model: String,
    pub task: String,
    pub composition: String,
    pub units: ProviderUnits,
    pub capabilities: Vec<String>,
    pub ase_version: String,
    pub numpy_version: String,
    pub model_sha256: String,
    pub asap_cutoff: bool,
}

/// Snapshot output in physical wire units (before reduced-unit conversion).
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ProviderEvaluation {
    pub energy_ev: f64,
    pub forces_ev_per_angstrom: Vec<[f64; 3]>,
    pub stress_ev_per_angstrom3: [[f64; 3]; 3],
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct InfoResponse {
    protocol: String,
    action: String,
    identity: ProviderIdentity,
}
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct EvaluateResponse {
    protocol: String,
    action: String,
    identity: ProviderIdentity,
    evaluation: ProviderEvaluation,
}

#[derive(Serialize)]
struct EvaluateRequest<'a> {
    protocol: &'static str,
    action: &'static str,
    atomic_numbers: &'a [u8],
    positions_angstrom: ScaledPositions<'a>,
    cell_angstrom: [[f64; 3]; 3],
    pbc: [bool; 3],
    asap_cutoff: bool,
}
// Borrow the caller's snapshot; triples are scaled on the stack during
// serialization instead of materializing another N-atom buffer or Value tree.
struct ScaledPositions<'a> {
    positions: &'a [f64],
    sigma: f64,
}
impl Serialize for ScaledPositions<'_> {
    fn serialize<S: serde::Serializer>(
        &self,
        serializer: S,
    ) -> std::result::Result<S::Ok, S::Error> {
        use serde::ser::SerializeSeq;
        let mut seq = serializer.serialize_seq(Some(self.positions.len() / 3))?;
        for p in self.positions.chunks_exact(3) {
            seq.serialize_element(&[p[0] * self.sigma, p[1] * self.sigma, p[2] * self.sigma])?;
        }
        seq.end()
    }
}

fn failure(message: impl Into<String>) -> ModenaError {
    // Parser diagnostics may quote arbitrary response keys. Bound every
    // provider diagnostic, not only the captured subprocess stderr.
    let mut message = message.into();
    if message.len() > 4096 {
        let mut end = 4096;
        while !message.is_char_boundary(end) {
            end -= 1;
        }
        message.truncate(end);
    }
    ModenaError::Physics(format!("ase_emt: {message}"))
}
fn hash_valid(hash: &str) -> bool {
    hash.len() == 64
        && hash
            .bytes()
            .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
}
fn validate_identity(id: &ProviderIdentity) -> Result<()> {
    if id.protocol != PROTOCOL
        || id.provider != "ase_emt"
        || id.model != "EMT"
        || id.task != "energy_forces_stress"
        || id.composition != "total"
        || id.asap_cutoff
        || id.units
            != (ProviderUnits {
                positions: "angstrom".into(),
                energy: "eV".into(),
                forces: "eV/angstrom".into(),
                stress: "eV/angstrom^3".into(),
            })
        || id.capabilities
            != [
                "energy",
                "forces",
                "stress",
                "full_pbc",
                "provider_neighbors",
                "conservative",
                "total_momentum",
            ]
        || id.ase_version.trim().is_empty()
        || id.numpy_version.trim().is_empty()
        || !hash_valid(&id.model_sha256)
    {
        return Err(failure("unsupported provider identity/units/capabilities"));
    }
    Ok(())
}

/// Discover the installed fixed calculator without accepting any model pins.
pub fn describe_ase_emt(python: &str, timeout_ms: u64) -> Result<ProviderIdentity> {
    let data = invoke(
        python,
        timeout_ms,
        &serde_json::json!({"protocol":PROTOCOL,"action":"info"}),
    )?;
    let response: InfoResponse = serde_json::from_slice(&data)
        .map_err(|e| failure(format!("invalid info response: {e}")))?;
    if response.protocol != PROTOCOL || response.action != "info" {
        return Err(failure("invalid info envelope"));
    }
    validate_identity(&response.identity)?;
    Ok(response.identity)
}

/// Pinned provider with physical length and energy scales.
pub struct AseEmtForceModel {
    spec: AseEmtSpec,
    factors: UnitFactors,
    identity: ProviderIdentity,
    session: Option<Mutex<session::Session>>,
}
impl AseEmtForceModel {
    pub fn new(spec: AseEmtSpec, factors: UnitFactors) -> Result<Self> {
        if spec.ase_version.trim().is_empty() {
            return Err(failure("ase_version pin is required"));
        }
        if !hash_valid(&spec.model_sha256) {
            return Err(failure("model_sha256 must be 64 lowercase hex digits"));
        }
        if spec.composition != "total" {
            return Err(failure("composition must be total"));
        }
        if spec.pbc != [true; 3] {
            return Err(failure("pbc must be fully periodic"));
        }
        if spec.asap_cutoff {
            return Err(failure("asap_cutoff must be false"));
        }
        if spec.atomic_numbers.is_empty()
            || spec
                .atomic_numbers
                .iter()
                .any(|z| ![13, 28, 29, 46, 47, 78, 79].contains(z))
        {
            return Err(failure(
                "atomic_numbers must contain supported standard metals",
            ));
        }
        if [
            factors.sigma,
            factors.epsilon,
            factors.mass_ref,
            factors.time,
            factors.temperature,
            factors.pressure,
            factors.charge,
        ]
        .iter()
        .any(|x| !x.is_finite() || *x <= 0.0)
        {
            return Err(failure("positive finite physical unit scales required"));
        }
        let mut session = match spec.transport {
            AseEmtTransport::OneShot => None,
            AseEmtTransport::Persistent => {
                Some(session::Session::new(&spec.python, spec.timeout_ms)?)
            }
        };
        let identity = if let Some(session) = session.as_mut() {
            let response: InfoResponse =
                session.request(&serde_json::json!({"protocol":PROTOCOL,"action":"info"}))?;
            if response.protocol != PROTOCOL || response.action != "info" {
                return Err(failure("invalid info envelope"));
            }
            validate_identity(&response.identity)?;
            response.identity
        } else {
            describe_ase_emt(&spec.python, spec.timeout_ms)?
        };
        if identity.ase_version != spec.ase_version || identity.model_sha256 != spec.model_sha256 {
            return Err(failure("ase_version or model_sha256 pin mismatch"));
        }
        Ok(Self {
            spec,
            factors,
            identity,
            session: session.map(Mutex::new),
        })
    }
    pub fn identity(&self) -> &ProviderIdentity {
        &self.identity
    }
    pub fn spec(&self) -> &AseEmtSpec {
        &self.spec
    }

    /// Evaluate a fresh snapshot without changing any system or accumulator.
    pub fn evaluate(&self, system: &System) -> Result<ProviderEvaluation> {
        let n = system.n();
        // Even the smallest JSON numbers/triples need at least ten bytes per
        // atom. Reject impossible frames before snapshot validation/encoding.
        if self.session.is_some() && n > MAX_STDOUT as usize / 10 {
            return Err(failure("snapshot exceeds session input frame limit"));
        }
        if system.atomic_numbers != self.spec.atomic_numbers || n != self.spec.atomic_numbers.len()
        {
            return Err(failure(
                "atomic_numbers do not match the pinned composition",
            ));
        }
        if system.positions.len() != 3 * n
            || system.forces.len() != 3 * n
            || system.charges.len() != n
            || system.charges.iter().any(|q| *q != 0.0)
            || system.constraint_count() != 0
        {
            return Err(failure("invalid snapshot buffers, charges or constraints"));
        }
        let cell = system
            .box_size
            .basis()
            .map(|v| v.0.map(|x| x * self.factors.sigma));
        let [a, b, c] = cell.map(crate::Vec3);
        let volume = a.dot(b.cross(c));
        if !volume.is_finite()
            || volume <= 0.0
            || system
                .positions
                .iter()
                .any(|x| !(x * self.factors.sigma).is_finite())
            || cell.iter().flatten().any(|x| !x.is_finite())
        {
            return Err(failure(
                "snapshot requires finite positions and positive cell volume",
            ));
        }
        let request = EvaluateRequest {
            protocol: PROTOCOL,
            action: "evaluate",
            atomic_numbers: &system.atomic_numbers,
            positions_angstrom: ScaledPositions {
                positions: &system.positions,
                sigma: self.factors.sigma,
            },
            cell_angstrom: cell,
            pbc: self.spec.pbc,
            asap_cutoff: self.spec.asap_cutoff,
        };
        if let Some(session) = &self.session {
            let mut session = session
                .lock()
                .map_err(|_| failure("session lock poisoned"))?;
            let response = session.request(&request)?;
            let result = self.validate_evaluation(response, n, volume);
            if result.is_err() {
                session.terminate();
            }
            result
        } else {
            let data = invoke(&self.spec.python, self.spec.timeout_ms, &request)?;
            let response = serde_json::from_slice(&data)
                .map_err(|e| failure(format!("invalid evaluation response: {e}")))?;
            self.validate_evaluation(response, n, volume)
        }
    }
    fn validate_evaluation(
        &self,
        response: EvaluateResponse,
        n: usize,
        volume: f64,
    ) -> Result<ProviderEvaluation> {
        if response.protocol != PROTOCOL
            || response.action != "evaluate"
            || response.identity != self.identity
        {
            return Err(failure("evaluation identity or envelope mismatch"));
        }
        validate_identity(&response.identity)?;
        let evaluation = response.evaluation;
        if evaluation.forces_ev_per_angstrom.len() != n
            || !evaluation.energy_ev.is_finite()
            || evaluation
                .forces_ev_per_angstrom
                .iter()
                .flatten()
                .chain(evaluation.stress_ev_per_angstrom3.iter().flatten())
                .any(|x| !x.is_finite())
        {
            return Err(failure(
                "invalid evaluation property shape or nonfinite value",
            ));
        }
        for i in 0..3 {
            for j in 0..3 {
                if (evaluation.stress_ev_per_angstrom3[i][j]
                    - evaluation.stress_ev_per_angstrom3[j][i])
                    .abs()
                    > 1e-12
                {
                    return Err(failure("stress must be symmetric"));
                }
            }
        }
        self.converted(&evaluation, volume)?;
        Ok(evaluation)
    }
    fn terminate_session(&self) {
        if let Some(session) = &self.session {
            session
                .lock()
                .unwrap_or_else(|e| e.into_inner())
                .terminate();
        }
    }
    fn converted(
        &self,
        e: &ProviderEvaluation,
        volume_angstrom3: f64,
    ) -> Result<(f64, Vec<f64>, f64)> {
        let energy_scale = EV_TO_KCAL_MOL / self.factors.epsilon;
        let energy = e.energy_ev * energy_scale;
        let forces: Vec<f64> = e
            .forces_ev_per_angstrom
            .iter()
            .flatten()
            .map(|f| f * self.factors.sigma * energy_scale)
            .collect();
        let trace = (0..3).map(|i| e.stress_ev_per_angstrom3[i][i]).sum::<f64>();
        let w = -volume_angstrom3 * trace * energy_scale;
        if !energy.is_finite() || !w.is_finite() || forces.iter().any(|x| !x.is_finite()) {
            return Err(failure("nonfinite converted properties"));
        }
        Ok((energy, forces, w))
    }
}
impl ForceModel for AseEmtForceModel {
    fn uses_provider_neighbors(&self) -> bool {
        true
    }
    fn conserves_total_momentum(&self, _system: &System) -> bool {
        true
    }
    fn compute(
        &self,
        system: &mut System,
        virial: &mut Virial,
        _neighbors: &mut NeighborList,
    ) -> Result<f64> {
        let evaluation = self.evaluate(system)?;
        let volume = system.box_size.volume() * self.factors.sigma.powi(3);
        let (energy, forces, w) = self
            .converted(&evaluation, volume)
            .inspect_err(|_| self.terminate_session())?;
        // Stage even the additions: overflow or malformed existing accumulators
        // must not leave a partially applied provider contribution.
        let accumulated: Vec<f64> = system
            .forces
            .iter()
            .zip(&forces)
            .map(|(a, b)| a + b)
            .collect();
        let total_w = virial.w + w;
        if !total_w.is_finite() || accumulated.iter().any(|x| !x.is_finite()) {
            self.terminate_session();
            return Err(failure("nonfinite force/virial accumulation"));
        }
        system.forces.copy_from_slice(&accumulated);
        virial.w = total_w;
        Ok(energy)
    }
}

fn diagnostic(file: &mut std::fs::File) -> String {
    let _ = file.seek(SeekFrom::Start(0));
    let mut bytes = Vec::new();
    let _ = file.take(4096).read_to_end(&mut bytes);
    String::from_utf8_lossy(&bytes).into_owned()
}
fn invoke(python: &str, timeout_ms: u64, request: &impl Serialize) -> Result<Vec<u8>> {
    if python.trim().is_empty() || timeout_ms == 0 {
        return Err(failure("python and positive timeout_ms required"));
    }
    let mut input = tempfile::tempfile()?;
    let mut output = tempfile::tempfile()?;
    let mut errors = tempfile::tempfile()?;
    serde_json::to_writer(&mut input, request).map_err(|e| failure(e.to_string()))?;
    input.flush()?;
    input.seek(SeekFrom::Start(0))?;
    let mut child = Command::new(python)
        .arg("-c")
        .arg(WORKER)
        .stdin(Stdio::from(input))
        .stdout(Stdio::from(output.try_clone()?))
        .stderr(Stdio::from(errors.try_clone()?))
        .spawn()
        .map_err(|e| failure(format!("cannot launch Python: {e}")))?;
    let started = Instant::now();
    let outcome = loop {
        let poll = (|| -> Result<Option<std::process::ExitStatus>> {
            if output.metadata()?.len() > MAX_STDOUT || errors.metadata()?.len() > MAX_STDERR {
                return Err(failure("subprocess output limit exceeded"));
            }
            if started.elapsed() >= Duration::from_millis(timeout_ms) {
                return Err(failure(format!("subprocess timeout after {timeout_ms} ms")));
            }
            Ok(child.try_wait()?)
        })();
        match poll {
            Ok(Some(status)) => break Ok(status),
            Err(e) => break Err(e),
            Ok(None) => std::thread::sleep(Duration::from_millis(5)),
        }
    };
    let status = match outcome {
        Ok(s) => s,
        Err(e) => {
            let _ = child.kill();
            let _ = child.wait();
            return Err(failure(format!("{e}; stderr: {}", diagnostic(&mut errors))));
        }
    };
    if !status.success() {
        return Err(failure(format!(
            "subprocess failed ({status}); stderr: {}",
            diagnostic(&mut errors)
        )));
    }
    if output.metadata()?.len() > MAX_STDOUT || errors.metadata()?.len() > MAX_STDERR {
        return Err(failure("subprocess output limit exceeded"));
    }
    output.seek(SeekFrom::Start(0))?;
    let mut data = Vec::new();
    output.take(MAX_STDOUT + 1).read_to_end(&mut data)?;
    if data.len() as u64 > MAX_STDOUT {
        return Err(failure("subprocess output limit exceeded"));
    }
    Ok(data)
}

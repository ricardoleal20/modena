//! Immutable, source-bound phase receipts and durable event streams.
use super::*;
use std::io::{BufRead, BufReader};

#[derive(serde::Serialize, serde::Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Event {
    pub sequence: u64,
    pub manifest_sha256: String,
    pub kind: String,
    pub data: serde_json::Value,
}
pub struct Journal {
    file: std::fs::File,
    epoch: String,
    sequence: u64,
}
impl Journal {
    pub fn create(path: &Path, epoch: String) -> Result<Self, CalibrationFailure> {
        Ok(Self {
            file: std::fs::OpenOptions::new()
                .write(true)
                .create_new(true)
                .open(path)
                .map_err(io)?,
            epoch,
            sequence: 0,
        })
    }
    /// Serialization streams directly into the file: no growing retained buffer.
    /// A partial write fails this phase; no later compute is admitted.
    pub fn record<T: serde::Serialize>(
        &mut self,
        kind: &str,
        data: &T,
    ) -> Result<(), CalibrationFailure> {
        #[derive(serde::Serialize)]
        struct Wire<'a, T> {
            sequence: u64,
            manifest_sha256: &'a str,
            kind: &'a str,
            data: &'a T,
        }
        serde_json::to_writer(
            &mut self.file,
            &Wire {
                sequence: self.sequence,
                manifest_sha256: &self.epoch,
                kind,
                data,
            },
        )
        .map_err(io)?;
        self.file.write_all(b"\n").map_err(io)?;
        self.file.sync_all().map_err(io)?;
        self.sequence = self
            .sequence
            .checked_add(1)
            .ok_or_else(|| invalid("journal sequence overflow"))?;
        Ok(())
    }
    pub fn owned_heap_bytes(&self) -> u64 {
        self.epoch.capacity() as u64
    }
}
pub fn read_json<T: serde::de::DeserializeOwned>(path: &Path) -> Result<T, CalibrationFailure> {
    let file = std::fs::File::open(path).map_err(io)?;
    serde_json::from_reader(file).map_err(|e| invalid(format!("{}: {e}", path.display())))
}
pub fn check_hash(path: &Path, expected: &str) -> Result<(), CalibrationFailure> {
    if expected.len() != 64
        || !expected.bytes().all(|c| c.is_ascii_hexdigit())
        || sha256(path)? != expected
    {
        return Err(invalid("artifact differs from authoritative phase digest"));
    }
    Ok(())
}
pub fn read_events(
    path: &Path,
    expected: &str,
    epoch: &str,
) -> Result<Vec<Event>, CalibrationFailure> {
    let mut events = Vec::new();
    visit_events(path, expected, epoch, |event| {
        events.push(event);
        Ok(())
    })?;
    Ok(events)
}
pub fn visit_events<F>(
    path: &Path,
    expected: &str,
    epoch: &str,
    mut visit: F,
) -> Result<(), CalibrationFailure>
where
    F: FnMut(Event) -> Result<(), CalibrationFailure>,
{
    check_hash(path, expected)?;
    for (sequence, line) in BufReader::new(std::fs::File::open(path).map_err(io)?)
        .lines()
        .enumerate()
    {
        let event: Event =
            serde_json::from_str(&line.map_err(io)?).map_err(|e| invalid(e.to_string()))?;
        if event.sequence != sequence as u64 || event.manifest_sha256 != epoch {
            return Err(invalid("journal sequence/source mismatch"));
        }
        visit(event)?;
    }
    Ok(())
}

#[derive(Clone, Debug, serde::Serialize, serde::Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Artifact {
    pub path: PathBuf,
    pub sha256: String,
}
impl Artifact {
    pub fn capture(path: &Path) -> Result<Self, CalibrationFailure> {
        Ok(Self {
            path: path.canonicalize().map_err(io)?,
            sha256: sha256(path)?,
        })
    }
    pub fn validate(&self) -> Result<(), CalibrationFailure> {
        check_hash(&self.path, &self.sha256)
    }
}
#[derive(Clone, Debug, serde::Serialize, serde::Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Usage {
    pub candidate: candidates::CandidateUsage,
    pub timing_calls: u64,
    pub profile_calls: u64,
}
impl Usage {
    pub fn capture(b: &CandidateBudget, timing: u64, profile: u64) -> Self {
        Self {
            candidate: b.usage(),
            timing_calls: timing,
            profile_calls: profile,
        }
    }
    pub fn resume(
        &self,
        p: &Protocol,
        parent_active_ns: u64,
    ) -> Result<CandidateBudget, CalibrationFailure> {
        let c = &self.candidate;
        if c.accuracy_calls > p.limits.max_accuracy_calls
            || self.timing_calls > p.limits.max_timing_calls
            || self.profile_calls > p.limits.max_profile_calls
            || c.candidate_real_pair_upper_bound
                .checked_add(c.reference.work.real_terms)
                .is_none_or(|x| x > p.limits.max_real_terms_per_family)
        {
            return Err(CalibrationFailure::BudgetExceeded(
                "continuing family counters exceed caps".into(),
            ));
        }
        let reference = ReferenceBudget::resume(
            &p.limits,
            c.reference.work.clone(),
            c.reference.calls,
            std::time::Duration::from_nanos(parent_active_ns.max(c.reference.elapsed_ns)),
            0,
        )?;
        let mut b = CandidateBudget::new(&p.limits, reference);
        b.accuracy_calls = c.accuracy_calls;
        b.candidate_real_pair_upper_bound = c.candidate_real_pair_upper_bound;
        Ok(b)
    }
}
/// Supplemental record. Canonical raw scientific schemas remain unchanged.
/// Only the parent's expected digest confers import authority; embedded hashes
/// authenticate descendants after that external trust anchor is checked.
#[derive(Clone, Debug, serde::Serialize, serde::Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Receipt {
    pub schema: String,
    pub phase: String,
    pub prepared: PathBuf,
    pub manifest_sha256: String,
    pub prior: Option<Artifact>,
    pub inputs: BTreeMap<String, Artifact>,
    pub artifacts: BTreeMap<String, Artifact>,
    pub usage: Usage,
    pub failure: Option<CalibrationFailure>,
    pub finalization_failure: Option<CalibrationFailure>,
}
impl Receipt {
    pub(super) fn owned_heap_bytes(&self) -> u64 {
        let artifact = |a: &Artifact| a.path.capacity() as u64 + a.sha256.capacity() as u64;
        self.schema.capacity() as u64
            + self.phase.capacity() as u64
            + self.prepared.capacity() as u64
            + self.manifest_sha256.capacity() as u64
            + self.prior.as_ref().map_or(0, artifact)
            + self
                .inputs
                .iter()
                .chain(&self.artifacts)
                .map(|(key, a)| 256 + key.capacity() as u64 + artifact(a))
                .sum::<u64>()
            + self
                .failure
                .iter()
                .chain(&self.finalization_failure)
                .map(candidates::failure_bytes)
                .sum::<u64>()
    }
    pub fn load(path: &Path, expected: &str) -> Result<Self, CalibrationFailure> {
        Self::load_depth(path, expected, 0)
    }
    fn load_depth(path: &Path, expected: &str, depth: usize) -> Result<Self, CalibrationFailure> {
        if depth > 16 {
            return Err(invalid("phase chain too deep"));
        }
        check_hash(path, expected)?;
        let r: Self = read_json(path)?;
        validate_prepared(&r.prepared)?;
        if r.schema != "modena-spme-phase-v1"
            || r.manifest_sha256 != sha256(&r.prepared.join("manifest.json"))?
        {
            return Err(invalid("phase prepared/source identity mismatch"));
        }
        for artifact in r.artifacts.values() {
            artifact.validate()?;
        }
        for input in r.inputs.values() {
            Self::load_depth(&input.path, &input.sha256, depth + 1)?;
        }
        if let Some(prior) = &r.prior {
            let previous = Self::load_depth(&prior.path, &prior.sha256, depth + 1)?;
            if previous.manifest_sha256 != r.manifest_sha256 {
                return Err(invalid("phase chain changed prepared epoch"));
            }
            let old = &previous.usage;
            let new = &r.usage;
            if new.candidate.accuracy_calls < old.candidate.accuracy_calls
                || new.candidate.reference.calls < old.candidate.reference.calls
                || new.candidate.reference.work.real_terms < old.candidate.reference.work.real_terms
                || new.candidate.reference.work.atom_k_terms
                    < old.candidate.reference.work.atom_k_terms
                || new.candidate.candidate_real_pair_upper_bound
                    < old.candidate.candidate_real_pair_upper_bound
                || new.candidate.reference.elapsed_ns < old.candidate.reference.elapsed_ns
                || new.timing_calls < old.timing_calls
                || new.profile_calls < old.profile_calls
            {
                return Err(invalid("phase chain reset continuing usage"));
            }
        }
        Ok(r)
    }
    pub fn artifact(&self, key: &str) -> Result<&Artifact, CalibrationFailure> {
        self.artifacts
            .get(key)
            .ok_or_else(|| invalid(format!("missing {key} phase artifact")))
    }
    pub fn require_success(&self) -> Result<(), CalibrationFailure> {
        if let Some(f) = self.failure.as_ref().or(self.finalization_failure.as_ref()) {
            return Err(f.clone());
        }
        Ok(())
    }
    pub fn write_new(&self, path: &Path) -> Result<Artifact, CalibrationFailure> {
        super::write_new(path, &json_bytes(self)?)?;
        Artifact::capture(path)
    }
}

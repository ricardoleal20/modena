//! Opt-in representative phase diagnostics, independent of primary cost/ranking.
//! Task7 supplies the authoritative campaign journal and phase-chain adapter.
use super::*;
use modena_core::{
    force_model::{ForceModel, Virial},
    spme::{mesh_for_box, SpmeProfile},
    system::SimBox,
};

/// Diagnostic-only ID, outside the frozen selectable IDs 0..125.
pub const BASELINE_PROFILE_ID: usize = 126;
pub const PHASE_NAMES: [&str; 12] = [
    "neighbor",
    "real_space",
    "pair_correction",
    "workspace_setup",
    "influence_setup",
    "spread",
    "fft_forward",
    "influence",
    "fft_inverse",
    "gather",
    "self_background",
    "total",
];
const CALLS_PER_MODEL: u64 = 106;
const CONTROL_CALLS: u64 = 42;

/// No unchecked constructor/deserialization: only baseline controls or frozen
/// confirmation can authorize an offline profile model.
#[derive(serde::Serialize)]
pub struct ValidatedProfileModel {
    candidate_id: usize,
    parameters: Parameters,
    snapshot: Snapshot,
    prepared_manifest_sha256: String,
}
impl ValidatedProfileModel {
    pub(super) fn owned_heap_bytes(&self) -> u64 {
        self.prepared_manifest_sha256.capacity() as u64 + candidates::snapshot_bytes(&self.snapshot)
    }
    pub fn candidate_id(&self) -> usize {
        self.candidate_id
    }
    pub fn parameters(&self) -> &Parameters {
        &self.parameters
    }
    pub fn snapshot_id(&self) -> &str {
        &self.snapshot.id
    }
}
fn representative(s: &[Snapshot]) -> Result<&Snapshot, CalibrationFailure> {
    s.iter()
        .find(|s| s.role == FixtureRole::Selection)
        .ok_or_else(|| invalid("profile representative missing"))
}
pub fn confirmed_profile_models(
    prepared: &Path,
    frozen: &FrozenSelection,
    confirmation: &ConfirmationOutcome,
) -> Result<Vec<ValidatedProfileModel>, CalibrationFailure> {
    let (_, snapshots, _) = validate_prepared(prepared)?;
    let snapshot = representative(&snapshots)?;
    let digest = sha256(&prepared.join("manifest.json"))?;
    Ok(frozen
        .profile_candidates(prepared, confirmation)?
        .into_iter()
        .map(|(id, parameters)| ValidatedProfileModel {
            candidate_id: id,
            parameters,
            snapshot: snapshot.clone(),
            prepared_manifest_sha256: digest.clone(),
        })
        .collect())
}
fn cap(message: &str) -> CalibrationFailure {
    CalibrationFailure::BudgetExceeded(message.into())
}
fn plus(a: u64, b: u64) -> Result<u64, CalibrationFailure> {
    a.checked_add(b)
        .ok_or_else(|| cap("profile accounting overflow"))
}
fn product(a: u64, b: u64) -> Result<u64, CalibrationFailure> {
    a.checked_mul(b)
        .ok_or_else(|| cap("profile accounting overflow"))
}
fn project(
    calls: u64,
    current: u64,
    limit: u64,
    atoms: usize,
    p: &Protocol,
    b: &mut CandidateBudget,
) -> Result<(u64, u64), CalibrationFailure> {
    b.reference.checkpoint()?;
    let total = plus(current, calls)?;
    let pairs = product(atoms as u64, atoms.saturating_sub(1) as u64)? / 2;
    if pairs > p.limits.max_real_terms_per_call {
        return Err(cap("profile per-call pair work"));
    }
    let pairs = plus(b.candidate_real_pair_upper_bound, product(pairs, calls)?)?;
    if total > limit
        || plus(pairs, b.reference.work.real_terms)? > p.limits.max_real_terms_per_family
    {
        return Err(cap("profile call/family pair cap"));
    }
    Ok((total, pairs))
}
fn reserve_accuracy(
    atoms: usize,
    p: &Protocol,
    b: &mut CandidateBudget,
) -> Result<(), CalibrationFailure> {
    let (calls, pairs) = project(
        1,
        b.accuracy_calls,
        p.limits.max_accuracy_calls,
        atoms,
        p,
        b,
    )?;
    b.accuracy_calls = calls;
    b.candidate_real_pair_upper_bound = pairs;
    Ok(())
}
#[derive(Clone, Debug, Default, serde::Serialize, serde::Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ProfileBudget {
    pub calls: u64,
}
impl ProfileBudget {
    fn reserve(
        &mut self,
        atoms: usize,
        p: &Protocol,
        b: &mut CandidateBudget,
    ) -> Result<(), CalibrationFailure> {
        let (calls, pairs) = project(1, self.calls, p.limits.max_profile_calls, atoms, p, b)?;
        self.calls = calls;
        b.candidate_real_pair_upper_bound = pairs;
        Ok(())
    }
}
/// Full raw functional control; this is not an independent accuracy gate.
#[derive(serde::Serialize)]
pub struct BaselineControl {
    pub label: String,
    pub step: f64,
    pub evaluation: Evaluation,
    pub real_pairs: Vec<(usize, usize)>,
    pub mesh: [usize; 3],
    pub analytic: Option<f64>,
    pub finite_difference: Option<f64>,
}
#[derive(serde::Serialize)]
pub struct BaselineValidation {
    pub model: Option<ValidatedProfileModel>,
    pub controls: Vec<BaselineControl>,
    pub failure: Option<CalibrationFailure>,
}
fn membership(s: &Snapshot, rc: f64) -> Vec<(usize, usize)> {
    let mut pairs = Vec::new();
    for i in 0..s.charges.len() {
        for j in 0..i {
            let r2 = (0..3)
                .map(|d| {
                    let x = s.positions[i][d] - s.positions[j][d];
                    let x = x - s.cell[d] * (x / s.cell[d]).round();
                    x * x
                })
                .sum::<f64>();
            if r2 < rc * rc {
                pairs.push((j, i));
            }
        }
    }
    pairs
}
fn near(a: f64, b: f64) -> bool {
    a.is_finite() && b.is_finite() && (a - b).abs() <= 1e-12 * b.abs().max(1.)
}
fn functional_failure(message: &str) -> CalibrationFailure {
    CalibrationFailure::ReferencePrecisionLimited(message.into())
}
/// Baseline analytic functional validation on the timing representative. Caller
/// must not invoke after a failed family reference/physics/budget phase. Every
/// completed raw control reaches the observer before the next risky compute.
#[derive(serde::Serialize)]
pub enum BaselineEvent<'a> {
    Admitted {
        label: &'a str,
        usage: super::candidates::CandidateUsage,
    },
    Completed(&'a BaselineControl),
}
pub fn validate_baseline_with_observer<F>(
    prepared: &Path,
    b: &mut CandidateBudget,
    mut observer: F,
) -> Result<BaselineValidation, CalibrationFailure>
where
    F: FnMut(BaselineEvent<'_>) -> Result<(), CalibrationFailure>,
{
    let (p, snapshots, _) = validate_prepared(prepared)?;
    let s = representative(&snapshots)?;
    let parameters = Parameters {
        kappa: 1. / p.a,
        real_cutoff: 0.4 * s.cell.into_iter().fold(f64::INFINITY, f64::min),
        kspace_cutoff: 12. / p.a,
    };
    let digest = sha256(&prepared.join("manifest.json"))?;
    project(
        CONTROL_CALLS,
        b.accuracy_calls,
        p.limits.max_accuracy_calls,
        s.charges.len(),
        &p,
        b,
    )?;
    // Conservative retained controls/pair identities/full forces/snapshot copies;
    // opaque model/FFT storage is covered by the process allocator.
    let atoms = s.charges.len() as u64;
    let per_control = plus(product(atoms, 128)?, product(product(atoms, atoms)?, 16)?)?;
    if plus(
        plus(b.external_owned_heap_bytes, b.reference.owned_heap_baseline)?,
        product(per_control, CONTROL_CALLS + 4)?,
    )? > p.limits.max_owned_heap_bytes
    {
        return Err(cap("baseline control retained heap"));
    }
    candidates::checked_mesh(s.cell, parameters.kspace_cutoff, &p.limits)?;
    let mut out = BaselineValidation {
        model: None,
        controls: Vec::with_capacity(CONTROL_CALLS as usize),
        failure: None,
    };
    let result = (|| {
        let mut evaluator = SpmeEvaluator::with_protocol(s, parameters.clone(), &p)?;
        let central_pairs = membership(s, parameters.real_cutoff);
        let central_mesh = evaluator.actual_mesh();
        let mut evaluate = |snapshot: &Snapshot,
                            label: String,
                            step: f64,
                            preload: bool,
                            analytic: Option<f64>,
                            other: Option<f64>|
         -> Result<Evaluation, CalibrationFailure> {
            let pairs = membership(snapshot, parameters.real_cutoff);
            let mesh = mesh_for_box(&SimBox::from_array(snapshot.cell), parameters.kspace_cutoff);
            if pairs != central_pairs || mesh != central_mesh {
                return Err(invalid(
                    "MeshBoundary or real cutoff membership changed during baseline control",
                ));
            }
            reserve_accuracy(snapshot.charges.len(), &p, b)?;
            observer(BaselineEvent::Admitted {
                label: &label,
                usage: b.usage(),
            })?;
            evaluator.set_snapshot(snapshot)?;
            let evaluation = if preload {
                for (i, f) in evaluator.system.forces.iter_mut().enumerate() {
                    *f = 0.125 + (i as f64 * 0.3).sin();
                }
                let mut virial = Virial { w: 1.75 };
                let energy = evaluator
                    .model
                    .compute(&mut evaluator.system, &mut virial, &mut evaluator.neighbors)
                    .map_err(|e| invalid(e.to_string()))?;
                Evaluation {
                    energy,
                    strain_w: virial.w,
                    forces: evaluator
                        .system
                        .forces
                        .chunks_exact(3)
                        .map(|f| [f[0], f[1], f[2]])
                        .collect(),
                }
            } else {
                evaluator.evaluate()?
            };
            let fd = other.map(|ep| -(ep - evaluation.energy) / (2. * step));
            out.controls.push(BaselineControl {
                label,
                step,
                evaluation: evaluation.clone(),
                real_pairs: pairs,
                mesh,
                analytic,
                finite_difference: fd,
            });
            observer(BaselineEvent::Completed(
                out.controls
                    .last()
                    .ok_or_else(|| invalid("baseline control missing"))?,
            ))?;
            if let (Some(a), Some(fd)) = (analytic, fd) {
                if !fd.is_finite() || (a - fd).abs() > 2e-8 {
                    return Err(functional_failure(
                        "baseline analytic/FD residual exceeds 2e-8",
                    ));
                }
            }
            Ok(evaluation)
        };
        let central = evaluate(s, "central".into(), 0., false, None, None)?;
        let additive = evaluate(s, "preloaded-additive".into(), 0., true, None, None)?;
        if !near(additive.energy, central.energy)
            || !near(additive.strain_w, central.strain_w + 1.75)
            || additive
                .forces
                .iter()
                .flatten()
                .zip(central.forces.iter().flatten())
                .enumerate()
                .any(|(i, (a, z))| !near(*a, *z + (0.125 + (i as f64 * 0.3).sin())))
        {
            return Err(functional_failure("baseline additive mismatch"));
        }
        for h in [1e-5, 5e-6] {
            for atom in 0..3 {
                for axis in 0..3 {
                    let mut plus = s.clone();
                    let mut minus = s.clone();
                    plus.positions[atom][axis] += h;
                    minus.positions[atom][axis] -= h;
                    let ep = evaluate(
                        &plus,
                        format!("force-{atom}-{axis}-plus"),
                        h,
                        false,
                        None,
                        None,
                    )?
                    .energy;
                    evaluate(
                        &minus,
                        format!("force-{atom}-{axis}-minus"),
                        h,
                        false,
                        Some(central.forces[atom][axis]),
                        Some(ep),
                    )?;
                }
            }
            let mut energies = [0.; 2];
            for (i, lambda) in [1. + h, 1. - h].into_iter().enumerate() {
                let mut strained = s.clone();
                strained.cell = s.cell.map(|x| x * lambda);
                for r in &mut strained.positions {
                    for x in r {
                        *x *= lambda;
                    }
                }
                energies[i] = evaluate(
                    &strained,
                    format!("strain-{}", if i == 0 { "plus" } else { "minus" }),
                    h,
                    false,
                    if i == 1 { Some(central.strain_w) } else { None },
                    if i == 1 { Some(energies[0]) } else { None },
                )?
                .energy;
            }
        }
        b.reference.checkpoint()?;
        out.model = Some(ValidatedProfileModel {
            candidate_id: BASELINE_PROFILE_ID,
            parameters,
            snapshot: s.clone(),
            prepared_manifest_sha256: digest,
        });
        Ok(())
    })();
    if let Err(f) = result {
        out.failure = Some(f);
        out.model = None;
    }
    Ok(out)
}

#[derive(Clone, Debug, PartialEq, serde::Serialize, serde::Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum ProfileStage {
    Cold,
    Warmup,
    Batch,
}
#[derive(Clone, Debug, PartialEq, serde::Serialize, serde::Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ProfileSample {
    pub stage: ProfileStage,
    pub batch: Option<usize>,
    pub call: usize,
    pub record: PhaseRecord,
}
#[derive(serde::Serialize)]
pub struct ProfileOutcome {
    pub samples: Vec<ProfileSample>,
    pub failure: Option<CalibrationFailure>,
    pub profile_calls: u64,
    pub family_elapsed_ns: u64,
}
/// Both events are outside compute. Hard allocator exits cannot recover an
/// in-flight sample; Task7 journals admissions and completed raw records.
#[derive(serde::Serialize)]
pub enum ProfileEvent<'a> {
    Admitted {
        candidate_id: usize,
        stage: ProfileStage,
        batch: Option<usize>,
        call: usize,
        profile_calls: u64,
        usage: super::candidates::CandidateUsage,
    },
    Completed(&'a ProfileSample),
}
fn durations(p: &SpmeProfile) -> Result<[u64; 12], CalibrationFailure> {
    let fields = [
        p.neighbor,
        p.real_space,
        p.pair_correction,
        p.workspace_setup,
        p.influence_setup,
        p.spread,
        p.fft_forward,
        p.influence,
        p.fft_inverse,
        p.gather,
        p.self_background,
        p.total,
    ];
    let mut ns = [0; 12];
    for (dst, d) in ns.iter_mut().zip(fields) {
        *dst = u64::try_from(d.as_nanos()).map_err(|_| cap("phase duration overflow"))?;
    }
    Ok(ns)
}
fn profile_call(e: &mut SpmeEvaluator) -> Result<(SpmeProfile, bool), CalibrationFailure> {
    e.system.zero_forces();
    let mut virial = Virial::default();
    let (energy, profile) = e
        .model
        .compute_profiled(&mut e.system, &mut virial, &mut e.neighbors)
        .map_err(|e| invalid(e.to_string()))?;
    std::hint::black_box((energy, virial.w, e.system.forces.as_slice()));
    let finite =
        energy.is_finite() && virial.w.is_finite() && e.system.forces.iter().all(|f| f.is_finite());
    Ok((profile, finite))
}
/// Exactly106calls/model; all raw warmups and5x20batch calls retained. No
/// observer cost enters TimingRecord/ranking. Budgets continue across phases.
pub fn measure_profiles_with_observer<F>(
    prepared: &Path,
    models: &[ValidatedProfileModel],
    timing: &mut ProfileBudget,
    b: &mut CandidateBudget,
    mut observer: F,
) -> Result<ProfileOutcome, CalibrationFailure>
where
    F: FnMut(ProfileEvent<'_>) -> Result<(), CalibrationFailure>,
{
    if cfg!(debug_assertions) {
        return Err(invalid("phase diagnostics require release build"));
    }
    let (p, snapshots, _) = validate_prepared(prepared)?;
    let snapshot = representative(&snapshots)?;
    let digest = sha256(&prepared.join("manifest.json"))?;
    if models.is_empty() || (p.warmups, p.batches, p.batch_calls) != (5, 5, 20) {
        return Err(invalid("profile model roster/frozen schedule"));
    }
    for (i, m) in models.iter().enumerate() {
        if m.prepared_manifest_sha256 != digest
            || m.snapshot != *snapshot
            || models[..i]
                .iter()
                .any(|old| old.candidate_id == m.candidate_id)
        {
            return Err(invalid("profile token source/fixture/duplicate mismatch"));
        }
        candidates::checked_mesh(snapshot.cell, m.parameters.kspace_cutoff, &p.limits)?;
    }
    let calls = product(models.len() as u64, CALLS_PER_MODEL)?;
    project(
        calls,
        timing.calls,
        p.limits.max_profile_calls,
        snapshot.charges.len(),
        &p,
        b,
    )?;
    let retained = product(calls, plus(512, snapshot.id.len() as u64)?)?;
    if plus(
        plus(b.external_owned_heap_bytes, b.reference.owned_heap_baseline)?,
        retained,
    )? > p.limits.max_owned_heap_bytes
    {
        return Err(cap("profile retained record heap"));
    }
    let mut out = ProfileOutcome {
        samples: Vec::with_capacity(calls as usize),
        failure: None,
        profile_calls: timing.calls,
        family_elapsed_ns: 0,
    };
    let result = (|| {
        for m in models {
            let mut evaluator = None;
            for index in 0..CALLS_PER_MODEL as usize {
                let (stage, batch, call) = if index == 0 {
                    (ProfileStage::Cold, None, 0)
                } else if index <= 5 {
                    (ProfileStage::Warmup, None, index - 1)
                } else {
                    (
                        ProfileStage::Batch,
                        Some((index - 6) / 20),
                        (index - 6) % 20,
                    )
                };
                timing.reserve(snapshot.charges.len(), &p, b)?;
                observer(ProfileEvent::Admitted {
                    candidate_id: m.candidate_id,
                    stage: stage.clone(),
                    batch,
                    call,
                    profile_calls: timing.calls,
                    usage: b.usage(),
                })?;
                if evaluator.is_none() {
                    evaluator = Some(SpmeEvaluator::with_protocol(
                        snapshot,
                        m.parameters.clone(),
                        &p,
                    )?);
                }
                let (profile, finite) = profile_call(
                    evaluator
                        .as_mut()
                        .ok_or_else(|| invalid("profile evaluator"))?,
                )?;
                out.samples.push(ProfileSample {
                    stage: stage.clone(),
                    batch,
                    call,
                    record: PhaseRecord {
                        candidate_id: m.candidate_id,
                        snapshot_id: snapshot.id.clone(),
                        threads: p.threads,
                        cold: stage == ProfileStage::Cold,
                        calls: 1,
                        phases_ns: vec![durations(&profile)?],
                    },
                });
                observer(ProfileEvent::Completed(
                    out.samples
                        .last()
                        .ok_or_else(|| invalid("profile record"))?,
                ))?;
                if !finite {
                    return Err(functional_failure("nonfinite profiled evaluation"));
                }
                b.reference.checkpoint()?;
            }
        }
        Ok(())
    })();
    if let Err(f) = result {
        out.failure = Some(f);
    }
    out.profile_calls = timing.calls;
    out.family_elapsed_ns = u64::try_from(b.reference.elapsed().as_nanos())
        .map_err(|_| cap("profile elapsed overflow"))?;
    Ok(out)
}
/// Create-new per-call artifacts. Flush/sync each sample; partial artifacts never
/// imply a complete phase. Task7 retains the separate source/decision journal.
pub struct ProfileWriter {
    json: std::fs::File,
    csv: std::fs::File,
}
impl ProfileWriter {
    pub fn create(out: &Path) -> Result<Self, CalibrationFailure> {
        std::fs::create_dir(out).map_err(io)?;
        let open = |name| {
            std::fs::OpenOptions::new()
                .create_new(true)
                .write(true)
                .open(out.join(name))
                .map_err(io)
        };
        let mut writer = Self {
            json: open("phase-times.jsonl")?,
            csv: open("phase-times.csv")?,
        };
        writeln!(
            writer.csv,
            "candidate_id,snapshot_id,threads,stage,batch,call,{}",
            PHASE_NAMES.map(|n| format!("{n}_ns")).join(",")
        )
        .map_err(io)?;
        writer.csv.sync_data().map_err(io)?;
        Ok(writer)
    }
    pub fn completed(&mut self, s: &ProfileSample) -> Result<(), CalibrationFailure> {
        if s.record.calls != 1 || s.record.phases_ns.len() != 1 {
            return Err(invalid("writer requires individual profile sample"));
        }
        serde_json::to_writer(&mut self.json, s).map_err(io)?;
        writeln!(self.json).map_err(io)?;
        self.json.sync_data().map_err(io)?;
        let stage = match s.stage {
            ProfileStage::Cold => "cold",
            ProfileStage::Warmup => "warmup",
            ProfileStage::Batch => "batch",
        };
        let id = s.record.snapshot_id.replace('"', "\"\"");
        write!(
            self.csv,
            "{},\"{}\",{},{},{},{},",
            s.record.candidate_id,
            id,
            s.record.threads,
            stage,
            s.batch.map(|b| b.to_string()).unwrap_or_default(),
            s.call
        )
        .map_err(io)?;
        for (i, n) in s.record.phases_ns[0].iter().enumerate() {
            if i > 0 {
                write!(self.csv, ",").map_err(io)?;
            }
            write!(self.csv, "{n}").map_err(io)?;
        }
        writeln!(self.csv).map_err(io)?;
        self.csv.sync_data().map_err(io)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn phase_duration_order_is_protocol_order() {
        let d = std::time::Duration::from_nanos;
        let p = SpmeProfile {
            neighbor: d(1),
            real_space: d(2),
            pair_correction: d(3),
            workspace_setup: d(4),
            influence_setup: d(5),
            spread: d(6),
            fft_forward: d(7),
            influence: d(8),
            fft_inverse: d(9),
            gather: d(10),
            self_background: d(11),
            total: d(12),
        };
        assert_eq!(
            durations(&p).unwrap(),
            [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
        );
        assert_eq!(PHASE_NAMES[4], "influence_setup");
        assert_eq!(PHASE_NAMES[7], "influence");
    }
    #[test]
    fn profile_projection_checks_continuing_calls_work_and_overflow() {
        let p = frozen_protocol("special", 2, 1, true).unwrap();
        let mut b = CandidateBudget::new(&p.limits, ReferenceBudget::new(&p.limits));
        b.accuracy_calls = 13;
        b.candidate_real_pair_upper_bound = 29;
        b.reference.work.real_terms = 7;
        assert_eq!(project(106, 9, 1000, 3, &p, &mut b).unwrap(), (115, 347));
        assert!(project(106, 900, 1000, 3, &p, &mut b).is_err());
        b.candidate_real_pair_upper_bound = p.limits.max_real_terms_per_family - 7;
        assert!(project(1, 0, 1000, 3, &p, &mut b).is_err());
        b.candidate_real_pair_upper_bound = u64::MAX;
        assert!(project(1, 0, 1000, 3, &p, &mut b).is_err());
        assert!(project(1, u64::MAX, u64::MAX, 3, &p, &mut b).is_err());
        assert_eq!(b.accuracy_calls, 13);
        assert_eq!(b.reference.calls, 0);
        let mut timing = ProfileBudget { calls: 1000 };
        assert!(timing.reserve(3, &p, &mut b).is_err());
        assert_eq!(timing.calls, 1000);
    }
    #[test]
    fn phase_writer_keeps_individual_order_and_escaping() {
        let out = std::env::temp_dir().join(format!(
            "modena-profile-writer-{}-{}",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        let mut writer = ProfileWriter::create(&out).unwrap();
        let sample = ProfileSample {
            stage: ProfileStage::Batch,
            batch: Some(2),
            call: 3,
            record: PhaseRecord {
                candidate_id: BASELINE_PROFILE_ID,
                snapshot_id: "fixture,\"quoted\"".into(),
                threads: 1,
                cold: false,
                calls: 1,
                phases_ns: vec![[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]],
            },
        };
        writer.completed(&sample).unwrap();
        assert_eq!(
            serde_json::from_str::<ProfileSample>(
                std::fs::read_to_string(out.join("phase-times.jsonl"))
                    .unwrap()
                    .trim()
            )
            .unwrap(),
            sample
        );
        let csv = std::fs::read_to_string(out.join("phase-times.csv")).unwrap();
        assert!(csv.contains("\"fixture,\"\"quoted\"\"\""));
        assert!(csv.ends_with(",1,2,3,4,5,6,7,8,9,10,11,12\n"));
        assert!(ProfileWriter::create(&out).is_err());
        std::fs::remove_dir_all(out).unwrap();
    }
    /// Explicit registered release-only experiment. No implicit fixture preparation
    /// or source relabeling: caller supplies the fresh prepared epoch and new output.
    #[test]
    #[ignore]
    fn registered_profile_smoke() {
        let prepared =
            PathBuf::from(std::env::var("MODENA_B05_PROFILE_PREPARED").expect("prepared path"));
        let out = PathBuf::from(std::env::var("MODENA_B05_PROFILE_OUT").expect("new output path"));
        let (p, _, _) = validate_prepared(&prepared).unwrap();
        assert_eq!(p.family, "special-3");
        assert_eq!(p.threads, 1);
        assert_eq!(p.reference_ladder.len(), 6);
        std::fs::create_dir(&out).unwrap();
        let mut b = CandidateBudget::new(&p.limits, ReferenceBudget::new(&p.limits));
        let mut control_events = std::fs::OpenOptions::new()
            .create_new(true)
            .write(true)
            .open(out.join("baseline-events.jsonl"))
            .unwrap();
        let validation = validate_baseline_with_observer(&prepared, &mut b, |event| {
            serde_json::to_writer(&mut control_events, &event).map_err(io)?;
            writeln!(control_events).map_err(io)?;
            control_events.sync_data().map_err(io)
        })
        .unwrap();
        super::super::write_new(
            &out.join("baseline-validation.json"),
            &json_bytes(&validation).unwrap(),
        )
        .unwrap();
        assert!(validation.failure.is_none(), "{:?}", validation.failure);
        assert_eq!(validation.controls.len(), 42);
        assert_eq!(b.accuracy_calls, 42);
        let model = validation.model.as_ref().unwrap();
        assert_eq!(model.candidate_id(), BASELINE_PROFILE_ID);
        let models = std::slice::from_ref(model);
        // Upfront cap rejection performs no profile call and changes no work.
        let before = b.candidate_real_pair_upper_bound;
        let mut capped = ProfileBudget { calls: 900 };
        assert!(measure_profiles_with_observer(
            &prepared,
            models,
            &mut capped,
            &mut b,
            |_| panic!("rejected profile observed work")
        )
        .is_err());
        assert_eq!(capped.calls, 900);
        assert_eq!(b.candidate_real_pair_upper_bound, before);
        let stale = ValidatedProfileModel {
            candidate_id: model.candidate_id,
            parameters: model.parameters.clone(),
            snapshot: model.snapshot.clone(),
            prepared_manifest_sha256: "stale".into(),
        };
        assert!(measure_profiles_with_observer(
            &prepared,
            &[stale],
            &mut ProfileBudget::default(),
            &mut b,
            |_| panic!("stale token computed")
        )
        .is_err());
        assert!(measure_profiles_with_observer(
            &prepared,
            &[],
            &mut ProfileBudget::default(),
            &mut b,
            |_| panic!("empty roster computed")
        )
        .is_err());
        let duplicates = (0..2)
            .map(|_| ValidatedProfileModel {
                candidate_id: model.candidate_id,
                parameters: model.parameters.clone(),
                snapshot: model.snapshot.clone(),
                prepared_manifest_sha256: model.prepared_manifest_sha256.clone(),
            })
            .collect::<Vec<_>>();
        assert!(measure_profiles_with_observer(
            &prepared,
            &duplicates,
            &mut ProfileBudget::default(),
            &mut b,
            |_| panic!("duplicate token computed")
        )
        .is_err());
        assert_eq!(b.candidate_real_pair_upper_bound, before);
        let mut writer = ProfileWriter::create(&out.join("phases")).unwrap();
        let mut events = std::fs::OpenOptions::new()
            .create_new(true)
            .write(true)
            .open(out.join("profile-events.jsonl"))
            .unwrap();
        let mut timing = ProfileBudget::default();
        let measured =
            measure_profiles_with_observer(&prepared, models, &mut timing, &mut b, |event| {
                serde_json::to_writer(&mut events, &event).map_err(io)?;
                writeln!(events).map_err(io)?;
                events.sync_data().map_err(io)?;
                if let ProfileEvent::Completed(sample) = event {
                    writer.completed(sample)?;
                }
                Ok(())
            })
            .unwrap();
        super::super::write_new(
            &out.join("profile-outcome.json"),
            &json_bytes(&measured).unwrap(),
        )
        .unwrap();
        assert!(measured.failure.is_none(), "{:?}", measured.failure);
        assert_eq!(measured.samples.len(), 106);
        assert_eq!(timing.calls, 106);
        assert_eq!(
            measured
                .samples
                .iter()
                .filter(|s| s.stage == ProfileStage::Cold)
                .count(),
            1
        );
        assert_eq!(
            measured
                .samples
                .iter()
                .filter(|s| s.stage == ProfileStage::Warmup)
                .count(),
            5
        );
        for batch in 0..5 {
            let rows: Vec<_> = measured
                .samples
                .iter()
                .filter(|s| s.batch == Some(batch))
                .collect();
            assert_eq!(rows.len(), 20);
            assert!(rows.iter().enumerate().all(|(i, s)| s.call == i));
        }
        assert_eq!(
            std::fs::read_to_string(out.join("phases/phase-times.jsonl"))
                .unwrap()
                .lines()
                .count(),
            106
        );
        assert_eq!(
            std::fs::read_to_string(out.join("phases/phase-times.csv"))
                .unwrap()
                .lines()
                .count(),
            107
        );
        let failed =
            measure_profiles_with_observer(&prepared, models, &mut timing, &mut b, |event| {
                if matches!(event, ProfileEvent::Completed(_)) {
                    Err(CalibrationFailure::Io(
                        "injected completed-profile observer failure".into(),
                    ))
                } else {
                    Ok(())
                }
            })
            .unwrap();
        super::super::write_new(
            &out.join("observer-failure.json"),
            &json_bytes(&failed).unwrap(),
        )
        .unwrap();
        assert!(matches!(failed.failure, Some(CalibrationFailure::Io(_))));
        assert_eq!(failed.samples.len(), 1);
        assert_eq!(timing.calls, 107);
        assert_eq!(b.accuracy_calls, 42);
        assert_eq!(b.reference.calls, 0);
        assert_eq!(b.candidate_real_pair_upper_bound, 149 * 3);
        super::super::write_new(&out.join("usage.json"), &json_bytes(&b.usage()).unwrap()).unwrap();
    }
}

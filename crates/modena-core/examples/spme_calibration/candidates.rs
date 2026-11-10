#![allow(clippy::ptr_arg)] // Vec capacities are part of retained-heap accounting.
use super::*;
use modena_core::{spme::mesh_for_box, system::SimBox};

pub(super) fn checked_mesh(
    cell: [f64; 3],
    kcut: f64,
    limits: &Limits,
) -> Result<[usize; 3], CalibrationFailure> {
    // Bound production next-5-smooth's input before its unchecked integer loop.
    if cell.iter().any(|x| !x.is_finite() || *x <= 0.) || !kcut.is_finite() || kcut <= 0. {
        return Err(invalid("mesh needs finite positive cell and kcut"));
    }
    for edge in cell {
        let raw = (edge * kcut / std::f64::consts::PI).ceil();
        if !raw.is_finite() || raw > limits.max_mesh_points.min(1_048_576) as f64 {
            return Err(CalibrationFailure::BudgetExceeded(
                "raw mesh axis exceeds cap".into(),
            ));
        }
    }
    let mesh = mesh_for_box(&SimBox::from_array(cell), kcut);
    check_work(
        mesh,
        &WorkCounts {
            real_terms: 0,
            atom_k_terms: 0,
        },
        0,
        limits,
    )?;
    Ok(mesh)
}
/// Exact 126-tuple grid, canonical kcut aliases are actual reduced inverse lengths.
pub fn candidate_grid(
    p: &Protocol,
    cell: [f64; 3],
) -> Result<Vec<CandidateReport>, CalibrationFailure> {
    if p.kappa_a != [1.2, 1.6, 2., 2.4, 3., 3.6]
        || p.cutoff_lmin != [0.30, 0.40, 0.49]
        || p.kcut_a != [8., 12., 16., 24., 32., 48., 64.]
        || p.a != 2.2
        || p.limits.max_candidates < 126
        || p.limits.max_mesh_points > 1_048_576
    {
        return Err(invalid(
            "candidate grid/units/caps differ from frozen protocol",
        ));
    }
    let lmin = cell.into_iter().fold(f64::INFINITY, f64::min);
    let mut result: Vec<CandidateReport> = Vec::with_capacity(126);
    for ka in &p.kappa_a {
        for ratio in &p.cutoff_lmin {
            for kca in &p.kcut_a {
                let parameters = Parameters {
                    kappa: ka / p.a,
                    real_cutoff: ratio * lmin,
                    kspace_cutoff: kca / p.a,
                };
                validate_parameters(&parameters, cell)?;
                let actual_mesh = checked_mesh(cell, parameters.kspace_cutoff, &p.limits)?;
                if let Some(existing) = result.iter_mut().find(|row| {
                    row.parameters.kappa == parameters.kappa
                        && row.parameters.real_cutoff == parameters.real_cutoff
                        && row.actual_mesh == actual_mesh
                }) {
                    existing.aliases.push(parameters.kspace_cutoff);
                    continue;
                }
                result.push(CandidateReport {
                    id: result.len(),
                    family: p.family.clone(),
                    parameters: parameters.clone(),
                    actual_mesh,
                    aliases: vec![parameters.kspace_cutoff],
                    rows: vec![],
                    timing: None,
                    status: Status::Pending,
                });
            }
        }
    }
    Ok(result)
}
/// A complete every-snapshot/every-shift reduction; never ranks or selects.
pub fn feasible_targets(
    c: &CandidateReport,
    snapshot_ids: &[&str],
    p: &Protocol,
) -> Result<Vec<f64>, CalibrationFailure> {
    if c.family != p.family
        || snapshot_ids.is_empty()
        || snapshot_ids
            .iter()
            .collect::<std::collections::BTreeSet<_>>()
            .len()
            != snapshot_ids.len()
        || c.rows.len() != snapshot_ids.len() * p.translations.len()
    {
        return Err(invalid("incomplete/mixed-family candidate rows"));
    }
    for id in snapshot_ids {
        for shift in 0..p.translations.len() {
            if c.rows
                .iter()
                .filter(|r| r.snapshot_id == *id && r.shift_index == shift)
                .count()
                != 1
            {
                return Err(invalid("missing/duplicate candidate snapshot shift"));
            }
        }
    }
    let mut passing = vec![];
    for target in &p.targets {
        let mut pass = c.status == Status::Passed;
        for row in &c.rows {
            pass &= row.status == Status::Passed && passes_target(&row.metrics, *target, p, true)?;
        }
        if pass {
            passing.push(*target);
        }
    }
    Ok(passing)
}

#[derive(Clone, Debug, serde::Serialize, serde::Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ShiftReference {
    pub shift_index: usize,
    pub evaluation: Evaluation,
    pub delta: ErrorMetrics,
    pub audit: ReferenceAudit,
}
#[derive(Clone, Debug, serde::Serialize, serde::Deserialize)]
#[serde(deny_unknown_fields)]
pub struct TranslationCertificate {
    pub prepared_manifest_sha256: String,
    pub snapshot_sha256: String,
    pub actual_mesh: [usize; 3],
    pub settings: ReferenceSettings,
    pub shifts: Vec<ShiftReference>,
    pub status: Status,
    pub failure: Option<CalibrationFailure>,
    pub usage: ReferenceUsage,
}
#[derive(Clone, Debug, serde::Serialize, serde::Deserialize)]
#[serde(deny_unknown_fields)]
pub struct RawCandidateEvaluation {
    pub snapshot_id: String,
    pub shift_index: usize,
    pub evaluation: Evaluation,
}
#[derive(Clone, Debug, serde::Serialize, serde::Deserialize)]
#[serde(deny_unknown_fields)]
pub struct TranslationVariation {
    pub snapshot_id: String,
    #[serde(
        serialize_with = "serialize_finite",
        deserialize_with = "deserialize_finite"
    )]
    pub max_force_rms_from_first: f64,
    #[serde(
        serialize_with = "serialize_finite",
        deserialize_with = "deserialize_finite"
    )]
    pub max_atom_force_from_first: f64,
    #[serde(
        serialize_with = "serialize_finite",
        deserialize_with = "deserialize_finite"
    )]
    pub energy_range: f64,
    #[serde(
        serialize_with = "serialize_finite",
        deserialize_with = "deserialize_finite"
    )]
    pub w_range: f64,
}
#[derive(Clone, Debug, serde::Serialize, serde::Deserialize)]
#[serde(deny_unknown_fields)]
pub struct CandidateUsage {
    pub accuracy_calls: u64,
    pub candidate_real_pair_upper_bound: u64,
    pub reference: ReferenceUsage,
    pub process_requested_rust_heap_bytes: usize,
    pub peak_process_requested_rust_heap_bytes: usize,
}
#[derive(Clone, Debug, serde::Serialize, serde::Deserialize)]
#[serde(deny_unknown_fields)]
pub struct CandidateOutcome {
    pub report: CandidateReport,
    pub raw: Vec<RawCandidateEvaluation>,
    pub translation_variation: Vec<TranslationVariation>,
    pub failure: Option<CalibrationFailure>,
    pub usage: CandidateUsage,
}
/// One context across every candidate. Reference terms/calls retain the original
/// ReferenceBudget meaning; candidate calls and real pair bounds are separate.
/// External baseline charges OTHER retained vector/string capacities (including
/// prior CandidateOutcomes and outer containers), never global allocator totals.
pub struct CandidateBudget {
    pub reference: ReferenceBudget,
    pub accuracy_calls: u64,
    pub candidate_real_pair_upper_bound: u64,
    pub external_owned_heap_bytes: u64,
    limits: Limits,
}
impl CandidateBudget {
    pub fn new(limits: &Limits, reference: ReferenceBudget) -> Self {
        Self {
            reference,
            accuracy_calls: 0,
            candidate_real_pair_upper_bound: 0,
            external_owned_heap_bytes: 0,
            limits: limits.clone(),
        }
    }
    pub fn usage(&self) -> CandidateUsage {
        CandidateUsage {
            accuracy_calls: self.accuracy_calls,
            candidate_real_pair_upper_bound: self.candidate_real_pair_upper_bound,
            reference: usage(&self.reference),
            process_requested_rust_heap_bytes: allocator::requested_heap(),
            peak_process_requested_rust_heap_bytes: allocator::peak_requested_heap(),
        }
    }
}
fn usage(b: &ReferenceBudget) -> ReferenceUsage {
    ReferenceUsage {
        work: b.work.clone(),
        calls: b.calls,
        elapsed_ns: b.elapsed().as_nanos().min(u64::MAX as u128) as u64,
    }
}
fn cap(message: &str) -> CalibrationFailure {
    CalibrationFailure::BudgetExceeded(message.into())
}
fn add(a: u64, b: u64) -> Result<u64, CalibrationFailure> {
    a.checked_add(b)
        .ok_or_else(|| cap("candidate accounting overflow"))
}
fn mul(a: u64, b: u64) -> Result<u64, CalibrationFailure> {
    a.checked_mul(b)
        .ok_or_else(|| cap("candidate accounting overflow"))
}
fn bytes<T>(v: &Vec<T>) -> u64 {
    (v.capacity() * std::mem::size_of::<T>()) as u64
}
fn snapshot_bytes(s: &Snapshot) -> u64 {
    bytes(&s.positions)
        + bytes(&s.charges)
        + bytes(&s.scales)
        + (s.id.capacity() + s.family.capacity() + s.units.capacity()) as u64
}
fn protocol_bytes(p: &Protocol) -> u64 {
    (p.schema.capacity() + p.family.capacity() + p.fixture_set.capacity()) as u64
        + bytes(&p.targets)
        + bytes(&p.kappa_a)
        + bytes(&p.cutoff_lmin)
        + bytes(&p.kcut_a)
        + bytes(&p.reference_ladder)
        + bytes(&p.translations)
}
fn failure_bytes(f: &CalibrationFailure) -> u64 {
    match f {
        CalibrationFailure::InvalidInput(s)
        | CalibrationFailure::UnsupportedDomain(s)
        | CalibrationFailure::BudgetExceeded(s)
        | CalibrationFailure::ReferenceUnconverged(s)
        | CalibrationFailure::ReferenceStagnated(s)
        | CalibrationFailure::ReferencePrecisionLimited(s)
        | CalibrationFailure::RelativeUndefined(s)
        | CalibrationFailure::Io(s) => s.capacity() as u64,
        CalibrationFailure::TargetNotReached(_) | CalibrationFailure::ConfirmationFailed(_) => 0,
    }
}
impl CandidateOutcome {
    pub fn owned_heap_bytes(&self) -> Result<u64, CalibrationFailure> {
        let r = &self.report;
        let mut total = (r.family.capacity()) as u64
            + bytes(&r.aliases)
            + bytes(&r.rows)
            + bytes(&self.raw)
            + bytes(&self.translation_variation);
        for row in &r.rows {
            total = add(
                total,
                row.snapshot_id.capacity() as u64 + bytes(&row.passing_targets),
            )?;
        }
        for raw in &self.raw {
            total = add(
                total,
                raw.snapshot_id.capacity() as u64 + bytes(&raw.evaluation.forces),
            )?;
        }
        for v in &self.translation_variation {
            total = add(total, v.snapshot_id.capacity() as u64)?;
        }
        if let Some(f) = &self.failure {
            total = add(total, failure_bytes(f))?;
        }
        Ok(total)
    }
}
impl TranslationCertificate {
    pub fn owned_heap_bytes(&self) -> Result<u64, CalibrationFailure> {
        let mut total = (self.prepared_manifest_sha256.capacity() + self.snapshot_sha256.capacity())
            as u64
            + bytes(&self.shifts);
        for row in &self.shifts {
            total = add(total, bytes(&row.evaluation.forces))?;
        }
        if let Some(f) = &self.failure {
            total = add(total, failure_bytes(f))?;
        }
        Ok(total)
    }
}
fn shifted(s: &Snapshot, mesh: [usize; 3], fraction: [f64; 3]) -> Snapshot {
    let mut result = s.clone();
    for pos in &mut result.positions {
        for axis in 0..3 {
            pos[axis] = (pos[axis] + s.cell[axis] / mesh[axis] as f64 * fraction[axis])
                .rem_euclid(s.cell[axis]);
        }
    }
    result
}
fn valid_certificate(
    c: &TranslationCertificate,
    epoch: &str,
    hash: &str,
    mesh: [usize; 3],
    level: &ReferenceLevel,
    p: &Protocol,
    b: &ReferenceBudget,
) -> bool {
    c.prepared_manifest_sha256 == epoch
        && c.snapshot_sha256 == hash
        && c.actual_mesh == mesh
        && c.settings == level.settings
        && c.status == Status::Passed
        && c.failure.is_none()
        && c.shifts.len() == 6
        && c.usage.calls <= b.calls
        && c.usage.work.real_terms <= b.work.real_terms
        && c.usage.work.atom_k_terms <= b.work.atom_k_terms
        && c.shifts.iter().enumerate().all(|(i, row)| {
            row.shift_index == i + 1
                && metrics(
                    &level.evaluation,
                    &row.evaluation,
                    level.evaluation.forces.len(),
                    p.tolerances.relative_denominator_floor / (p.a * p.a),
                )
                .is_ok_and(|m| m == row.delta && invariant(&m, p))
        })
}
fn reference_passes(record: &ReferenceRecord) -> bool {
    record.chosen_level.is_some()
        && matches!(record.status, Status::Passed | Status::RelativeUndefined)
}
fn invariant(delta: &ErrorMetrics, p: &Protocol) -> bool {
    let f = p.q_unit * p.q_unit / (p.a * p.a);
    let e = p.q_unit * p.q_unit / p.a;
    let t = &p.tolerances;
    let rms = if delta.reference_force_rms > t.relative_denominator_floor * f {
        (t.reference_force_rms_abs_units * f)
            .min(t.reference_force_target_fraction * 1e-4 * delta.reference_force_rms)
    } else {
        t.reference_force_rms_abs_units * f
    };
    delta.force_rms_abs <= rms.min(1e-9 * f)
        && delta.max_atom_force_error <= (t.reference_force_max_units * f).min(1e-9 * f)
        && delta.energy_per_atom_abs <= (t.reference_energy_per_atom_units * e).min(1e-9 * e)
        && delta.w_per_atom_abs <= (t.reference_w_per_atom_units * e).min(1e-9 * e)
}
/// Conservative checked enumeration work; spherical filtering can only reduce it.
fn reference_upper(
    s: &Snapshot,
    settings: &ReferenceSettings,
) -> Result<WorkCounts, CalibrationFailure> {
    let mut images = 1;
    let mut vectors = 1;
    for l in s.cell {
        let r = (settings.real_cutoff / l).ceil() + 1.;
        let k = (settings.kspace_cutoff * l / (2. * std::f64::consts::PI)).ceil();
        if !r.is_finite() || !k.is_finite() || r > i32::MAX as f64 || k > i32::MAX as f64 {
            return Err(cap("reference enumeration bound"));
        }
        images = mul(images, add(mul(r as u64, 2)?, 1)?)?;
        vectors = mul(vectors, add(mul(k as u64, 2)?, 1)?)?;
    }
    let n = s.charges.len() as u64;
    Ok(WorkCounts {
        real_terms: mul(mul(n, n)?, images)?,
        atom_k_terms: mul(n, vectors - 1)?,
    })
}
/// Whole selection+worst-case confirmation call preflight and distinct-mesh
/// selection invariance work preflight. Invoke before a full-grid sweep.
pub fn preflight_selection_sweep(
    p: &Protocol,
    snapshots: &[Snapshot],
    b: &CandidateBudget,
) -> Result<WorkCounts, CalibrationFailure> {
    validate_frozen(p, snapshots)?;
    let selection: Vec<_> = snapshots
        .iter()
        .filter(|s| s.role == FixtureRole::Selection)
        .collect();
    let first = selection
        .first()
        .ok_or_else(|| invalid("empty selection roster"))?;
    let grid = candidate_grid(p, first.cell)?;
    let count = snapshots
        .iter()
        .filter(|s| s.role != FixtureRole::RegressionOnly)
        .count() as u64;
    let calls = mul(mul(grid.len() as u64, count)?, 7)?;
    if add(b.accuracy_calls, calls)? > p.limits.max_accuracy_calls {
        return Err(cap("projected whole campaign accuracy call cap"));
    }
    let meshes: std::collections::BTreeSet<_> = grid.iter().map(|c| c.actual_mesh).collect();
    let final_level = p
        .reference_ladder
        .last()
        .ok_or_else(|| invalid("reference ladder"))?;
    let settings = ReferenceSettings {
        kappa: p.reference_kappa_a / p.a,
        real_cutoff: final_level[0] * p.a,
        kspace_cutoff: final_level[1] / p.a,
    };
    let mut work = b.reference.work.clone();
    for s in selection {
        let upper = reference_upper(s, &settings)?;
        check_work([1; 3], &upper, 0, &p.limits)?;
        let copies = mul(meshes.len() as u64, 6)?;
        work.real_terms = add(work.real_terms, mul(copies, upper.real_terms)?)?;
        work.atom_k_terms = add(work.atom_k_terms, mul(copies, upper.atom_k_terms)?)?;
    }
    if work.real_terms > p.limits.max_real_terms_per_family
        || work.atom_k_terms > p.limits.max_atom_k_terms_per_family
    {
        return Err(cap(
            "projected distinct-mesh oracle invariance family work cap",
        ));
    }
    Ok(work)
}

/// Evaluates one grid candidate against ALL selection snapshots and seven phases.
/// No holdout evaluation, cost ranking, timings or emitted configuration.
/// Prepared source/spec/fixture bytes are revalidated before any cache reuse.
/// The caller retains ReferenceOutcomes from one continuing ReferenceBudget.
/// Streaming completed evidence only, called outside any timed region and before
/// later work. A journal should durably commit each event before returning Ok.
/// An allocator hard exit cannot recover in-flight heap data; parents retain
/// already committed events. I/O errors stop the evaluator and preserve its partial outcome.
#[derive(serde::Serialize)]
pub enum CandidateEvent<'a> {
    ReferenceShift {
        certificate: &'a TranslationCertificate,
        shift: &'a ShiftReference,
        usage: CandidateUsage,
    },
    CandidateRow {
        row: &'a MetricRow,
        raw: &'a RawCandidateEvaluation,
        usage: CandidateUsage,
    },
    Complete {
        outcome: &'a CandidateOutcome,
    },
}

pub fn evaluate_prepared_candidate(
    prepared: &Path,
    candidate: &CandidateReport,
    references: &Vec<ReferenceOutcome>,
    b: &mut CandidateBudget,
    certificates: &mut Vec<TranslationCertificate>,
) -> Result<CandidateOutcome, CalibrationFailure> {
    evaluate_prepared_candidate_with_observer(
        prepared,
        candidate,
        references,
        b,
        certificates,
        |_| Ok(()),
    )
}

pub fn evaluate_prepared_candidate_with_observer<F>(
    prepared: &Path,
    candidate: &CandidateReport,
    references: &Vec<ReferenceOutcome>,
    b: &mut CandidateBudget,
    certificates: &mut Vec<TranslationCertificate>,
    mut observer: F,
) -> Result<CandidateOutcome, CalibrationFailure>
where
    F: FnMut(CandidateEvent<'_>) -> Result<(), CalibrationFailure>,
{
    let (p, snapshots, manifest) = validate_prepared(prepared)?;
    if b.limits != p.limits || p.limits.max_owned_heap_bytes != allocator::PROCESS_HEAP_CAP as u64 {
        return Err(invalid(
            "candidate context requires matching limits and registered process heap cap",
        ));
    }
    let epoch = sha256(&prepared.join("manifest.json"))?;
    let hashes: Vec<String> = snapshots
        .iter()
        .map(|s| manifest.fixture_sha256[&s.id].clone())
        .collect();
    drop(manifest); // no opaque BTreeMap storage is present during reference accounting
    let first = snapshots
        .iter()
        .find(|s| s.role == FixtureRole::Selection)
        .ok_or_else(|| invalid("selection snapshots"))?;
    let grid = candidate_grid(&p, first.cell)?;
    if !grid.iter().any(|c| {
        c.id == candidate.id
            && c.family == candidate.family
            && c.parameters == candidate.parameters
            && c.actual_mesh == candidate.actual_mesh
            && c.aliases == candidate.aliases
    }) || !candidate.rows.is_empty()
        || candidate.timing.is_some()
    {
        return Err(invalid(
            "candidate must be an untouched exact prepared grid row",
        ));
    }
    drop(grid);
    let mut out = CandidateOutcome {
        report: candidate.clone(),
        raw: vec![],
        translation_variation: vec![],
        failure: None,
        usage: b.usage(),
    };
    let original_baseline = b.reference.owned_heap_baseline;
    let result = (|| -> Result<(), CalibrationFailure> {
        b.reference.checkpoint()?;
        let selection_count = snapshots
            .iter()
            .filter(|s| s.role == FixtureRole::Selection)
            .count();
        let projected = mul(selection_count as u64, 7)?;
        if add(b.accuracy_calls, projected)? > p.limits.max_accuracy_calls {
            return Err(cap("candidate accuracy calls before work"));
        }
        if references.len() != selection_count {
            return Err(invalid("exact selection reference roster required"));
        }
        // Validate every prerequisite before computing the first candidate row.
        for (index, s) in snapshots
            .iter()
            .enumerate()
            .filter(|(_, s)| s.role == FixtureRole::Selection)
        {
            let matching: Vec<_> = references
                .iter()
                .filter(|r| {
                    r.record.snapshot_id == s.id && r.record.snapshot_sha256 == hashes[index]
                })
                .collect();
            if matching.len() != 1 {
                return Err(invalid("reference hash/ID roster mismatch"));
            }
            let r = matching[0];
            if r.failure.is_some() || !reference_passes(&r.record) {
                return Err(CalibrationFailure::ReferenceUnconverged(format!(
                    "{}: reference has no strict converged level",
                    s.id
                )));
            }
            let chosen = r
                .record
                .chosen_level
                .ok_or_else(|| invalid("reference level"))?;
            let level = r
                .record
                .levels
                .get(chosen)
                .ok_or_else(|| invalid("reference level bounds"))?;
            let last = p
                .reference_ladder
                .last()
                .ok_or_else(|| invalid("reference ladder"))?;
            if chosen + 1 != p.reference_ladder.len()
                || level.settings
                    != (ReferenceSettings {
                        kappa: p.reference_kappa_a / p.a,
                        real_cutoff: last[0] * p.a,
                        kspace_cutoff: last[1] / p.a,
                    })
                || r.usage.calls > b.reference.calls
                || r.usage.work.real_terms > b.reference.work.real_terms
                || r.usage.work.atom_k_terms > b.reference.work.atom_k_terms
                || r.usage.elapsed_ns > usage(&b.reference).elapsed_ns
            {
                return Err(invalid(
                    "reference settings or retained family usage mismatch",
                ));
            }
        }
        // Preflight all missing mesh certificates before any candidate evaluation.
        let mut projected_work = b.reference.work.clone();
        let mut projected_reference_calls = b.reference.calls;
        let mut projected_pairs = b.candidate_real_pair_upper_bound;
        for (index, s) in snapshots
            .iter()
            .enumerate()
            .filter(|(_, s)| s.role == FixtureRole::Selection)
        {
            let r = references
                .iter()
                .find(|r| r.record.snapshot_id == s.id)
                .ok_or_else(|| invalid("reference"))?;
            let level = &r.record.levels[r
                .record
                .chosen_level
                .ok_or_else(|| invalid("reference level"))?];
            let cached = certificates.iter().any(|c| {
                valid_certificate(
                    c,
                    &epoch,
                    &hashes[index],
                    candidate.actual_mesh,
                    level,
                    &p,
                    &b.reference,
                )
            });
            let n = s.charges.len() as u64;
            let pairs = mul(n, n - 1)? / 2;
            if pairs > p.limits.max_real_terms_per_call {
                return Err(cap("candidate real pair per-call cap"));
            }
            projected_pairs = add(projected_pairs, mul(7, pairs)?)?;
            if !cached {
                projected_reference_calls = add(projected_reference_calls, 6)?;
                let upper = reference_upper(s, &level.settings)?;
                check_work([1; 3], &upper, 0, &p.limits)?;
                projected_work.real_terms =
                    add(projected_work.real_terms, mul(6, upper.real_terms)?)?;
                projected_work.atom_k_terms =
                    add(projected_work.atom_k_terms, mul(6, upper.atom_k_terms)?)?;
            }
        }
        if add(projected_work.real_terms, projected_pairs)? > p.limits.max_real_terms_per_family
            || projected_work.atom_k_terms > p.limits.max_atom_k_terms_per_family
            || projected_reference_calls > p.limits.max_accuracy_calls
        {
            return Err(cap("candidate invariance projected family work/call cap"));
        }
        let mut fixed = add(b.external_owned_heap_bytes, original_baseline)?;
        fixed = add(
            fixed,
            protocol_bytes(&p)
                + bytes(&snapshots)
                + epoch.capacity() as u64
                + bytes(&hashes)
                + bytes(references),
        )?;
        for s in &snapshots {
            fixed = add(fixed, snapshot_bytes(s))?;
        }
        for h in &hashes {
            fixed = add(fixed, h.capacity() as u64)?;
        }
        for r in references {
            fixed = add(fixed, r.owned_heap_bytes()?)?;
        }
        for (index, s) in snapshots
            .iter()
            .enumerate()
            .filter(|(_, s)| s.role == FixtureRole::Selection)
        {
            let r = references
                .iter()
                .find(|r| r.record.snapshot_id == s.id)
                .ok_or_else(|| invalid("reference"))?;
            let level = &r.record.levels[r
                .record
                .chosen_level
                .ok_or_else(|| invalid("reference level"))?];
            let cached = certificates.iter().any(|c| {
                valid_certificate(
                    c,
                    &epoch,
                    &hashes[index],
                    candidate.actual_mesh,
                    level,
                    &p,
                    &b.reference,
                )
            });
            if !cached {
                certificates.push(TranslationCertificate {
                    prepared_manifest_sha256: epoch.clone(),
                    snapshot_sha256: hashes[index].clone(),
                    actual_mesh: candidate.actual_mesh,
                    settings: level.settings.clone(),
                    shifts: vec![],
                    status: Status::Pending,
                    failure: None,
                    usage: usage(&b.reference),
                });
                let current = certificates.len() - 1;
                for shift_index in 1..7 {
                    b.reference.checkpoint()?;
                    let translated = shifted(s, candidate.actual_mesh, p.translations[shift_index]);
                    let mut baseline = add(fixed, out.owned_heap_bytes()?)?;
                    baseline = add(baseline, bytes(certificates))?;
                    for c in certificates.iter() {
                        baseline = add(baseline, c.owned_heap_bytes()?)?;
                    }
                    b.reference.owned_heap_baseline = baseline;
                    let attempt = reference_with_budget(
                        &translated,
                        level.settings.clone(),
                        &p.limits,
                        &mut b.reference,
                    );
                    let (evaluation, audit) = match attempt {
                        Ok(value) => value,
                        Err(f) => {
                            certificates[current].failure = Some(f.clone());
                            certificates[current].status = Status::Failed;
                            certificates[current].usage = usage(&b.reference);
                            return Err(f);
                        }
                    };
                    let delta = metrics(
                        &level.evaluation,
                        &evaluation,
                        s.charges.len(),
                        p.tolerances.relative_denominator_floor / (p.a * p.a),
                    )?;
                    let passes = invariant(&delta, &p);
                    certificates[current].shifts.push(ShiftReference {
                        shift_index,
                        evaluation,
                        delta,
                        audit,
                    });
                    certificates[current].usage = usage(&b.reference);
                    observer(CandidateEvent::ReferenceShift {
                        certificate: &certificates[current],
                        shift: certificates[current]
                            .shifts
                            .last()
                            .ok_or_else(|| invalid("reference shift checkpoint"))?,
                        usage: b.usage(),
                    })?;
                    if !passes {
                        let f = CalibrationFailure::ReferenceUnconverged(format!(
                            "{}: oracle translation invariance failed at shift {}",
                            s.id, shift_index
                        ));
                        certificates[current].status = Status::Failed;
                        certificates[current].failure = Some(f.clone());
                        return Err(f);
                    }
                }
                certificates[current].status = Status::Passed;
            }
        }
        b.reference.owned_heap_baseline = original_baseline;
        let mut evaluator = SpmeEvaluator::with_protocol(first, candidate.parameters.clone(), &p)?;
        for s in snapshots
            .iter()
            .filter(|s| s.role == FixtureRole::Selection)
        {
            let r = references
                .iter()
                .find(|r| r.record.snapshot_id == s.id)
                .ok_or_else(|| invalid("reference"))?;
            let reference = &r.record.levels[r
                .record
                .chosen_level
                .ok_or_else(|| invalid("reference level"))?]
            .evaluation;
            let start = out.raw.len();
            for shift_index in 0..7 {
                b.reference.checkpoint()?;
                let translated = shifted(s, candidate.actual_mesh, p.translations[shift_index]);
                evaluator.set_snapshot(&translated)?;
                if evaluator.actual_mesh() != candidate.actual_mesh {
                    return Err(invalid("selection cell changed actual mesh"));
                }
                b.accuracy_calls = add(b.accuracy_calls, 1)?; // count attempted work, including a failed compute
                let n = s.charges.len() as u64;
                let pairs = mul(n, n - 1)? / 2;
                b.candidate_real_pair_upper_bound = add(b.candidate_real_pair_upper_bound, pairs)?;
                let evaluation = evaluator.evaluate()?;
                out.raw.push(RawCandidateEvaluation {
                    snapshot_id: s.id.clone(),
                    shift_index,
                    evaluation,
                });
                let evaluation = &out
                    .raw
                    .last()
                    .ok_or_else(|| invalid("raw candidate evaluation"))?
                    .evaluation;
                let m = metrics(
                    reference,
                    evaluation,
                    s.charges.len(),
                    p.tolerances.relative_denominator_floor * p.q_unit * p.q_unit / (p.a * p.a),
                )?;
                let mut passing_targets = vec![];
                for t in &p.targets {
                    if passes_target(&m, *t, &p, r.record.status == Status::Passed)? {
                        passing_targets.push(*t);
                    }
                }
                let status = if m.force_rms_rel.is_none() {
                    Status::RelativeUndefined
                } else if passing_targets.is_empty() {
                    Status::Failed
                } else {
                    Status::Passed
                };
                out.report.rows.push(MetricRow {
                    snapshot_id: s.id.clone(),
                    shift_index,
                    metrics: m,
                    passing_targets,
                    status,
                });
                observer(CandidateEvent::CandidateRow {
                    row: out
                        .report
                        .rows
                        .last()
                        .ok_or_else(|| invalid("row checkpoint"))?,
                    raw: out.raw.last().ok_or_else(|| invalid("raw checkpoint"))?,
                    usage: b.usage(),
                })?;
                b.reference.checkpoint()?;
            }
            let group = &out.raw[start..];
            let mut v = TranslationVariation {
                snapshot_id: s.id.clone(),
                max_force_rms_from_first: 0.,
                max_atom_force_from_first: 0.,
                energy_range: 0.,
                w_range: 0.,
            };
            let mut emin = group[0].evaluation.energy;
            let mut emax = emin;
            let mut wmin = group[0].evaluation.strain_w;
            let mut wmax = wmin;
            for row in group.iter().skip(1) {
                let m = metrics(&group[0].evaluation, &row.evaluation, s.charges.len(), 0.)?;
                v.max_force_rms_from_first = v.max_force_rms_from_first.max(m.force_rms_abs);
                v.max_atom_force_from_first =
                    v.max_atom_force_from_first.max(m.max_atom_force_error);
                emin = emin.min(row.evaluation.energy);
                emax = emax.max(row.evaluation.energy);
                wmin = wmin.min(row.evaluation.strain_w);
                wmax = wmax.max(row.evaluation.strain_w);
            }
            v.energy_range = emax - emin;
            v.w_range = wmax - wmin;
            out.translation_variation.push(v);
        }
        out.report.status = Status::Passed;
        let ids: Vec<_> = snapshots
            .iter()
            .filter(|s| s.role == FixtureRole::Selection)
            .map(|s| s.id.as_str())
            .collect();
        if feasible_targets(&out.report, &ids, &p)?.is_empty() {
            out.report.status = if out
                .report
                .rows
                .iter()
                .any(|r| r.status == Status::RelativeUndefined)
            {
                Status::RelativeUndefined
            } else {
                Status::TargetNotReached
            };
        }
        Ok(())
    })();
    b.reference.owned_heap_baseline = original_baseline;
    if let Err(f) = result {
        out.report.status = match f {
            CalibrationFailure::BudgetExceeded(_) => Status::BudgetExceeded,
            CalibrationFailure::ReferenceUnconverged(_) => Status::ReferenceUnconverged,
            _ => Status::Failed,
        };
        out.failure = Some(f);
    }
    out.usage = b.usage();
    if let Err(f) = observer(CandidateEvent::Complete { outcome: &out }) {
        out.failure = Some(f);
        out.report.status = Status::Failed;
    }
    Ok(out)
}

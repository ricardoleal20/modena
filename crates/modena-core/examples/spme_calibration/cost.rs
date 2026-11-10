//! Bounded ordinary costs, source-bound reference handoff and frozen confirmation.
use super::*;
use candidates::{bytes, snapshot_bytes};

/// Origin cannot be attached to arbitrary ReferenceOutcomes by a public constructor.
/// Persisted reload requires an authoritative journal/phase digest supplied by the
/// orchestrator, not a digest read from the untrusted reference file itself.
#[derive(serde::Serialize)]
pub struct PreparedReferences {
    prepared_manifest_sha256: String,
    role: FixtureRole,
    records: Vec<ReferenceOutcome>,
}
#[derive(serde::Deserialize)]
#[serde(deny_unknown_fields)]
struct ReferenceWire {
    prepared_manifest_sha256: String,
    role: FixtureRole,
    records: Vec<ReferenceOutcome>,
}
impl PreparedReferences {
    pub fn records(&self) -> &[ReferenceOutcome] {
        &self.records
    }
    pub fn failure(&self) -> Option<&CalibrationFailure> {
        self.records.iter().find_map(|r| r.failure.as_ref())
    }
    pub fn write_new(&self, path: &Path) -> Result<String, CalibrationFailure> {
        super::write_new(path, &json_bytes(self)?)?;
        sha256(path)
    }
    pub fn load_selection(
        prepared: &Path,
        path: &Path,
        authoritative_phase_sha256: &str,
        b: &mut CandidateBudget,
    ) -> Result<Self, CalibrationFailure> {
        if authoritative_phase_sha256.len() != 64 || sha256(path)? != authoritative_phase_sha256 {
            return Err(invalid(
                "reference phase digest differs from authoritative journal",
            ));
        }
        let wire: ReferenceWire =
            serde_json::from_slice(&std::fs::read(path).map_err(io)?).map_err(io)?;
        let result = Self {
            prepared_manifest_sha256: wire.prepared_manifest_sha256,
            role: wire.role,
            records: wire.records,
        };
        if result.role != FixtureRole::Selection {
            return Err(invalid("selection reference role"));
        }
        result.validate(prepared, b)?;
        let (p, s, _) = validate_prepared(prepared)?;
        preflight_selection_sweep(&p, &s, b)?;
        Ok(result)
    }
    fn validate(&self, prepared: &Path, b: &mut CandidateBudget) -> Result<(), CalibrationFailure> {
        let (_, snapshots, manifest) = validate_prepared(prepared)?;
        b.reference.checkpoint()?;
        if self.prepared_manifest_sha256 != sha256(&prepared.join("manifest.json"))? {
            return Err(invalid(
                "reference origin differs from current prepared source/protocol",
            ));
        }
        let roster: Vec<_> = snapshots.iter().filter(|s| s.role == self.role).collect();
        if roster.len() != self.records.len() {
            return Err(invalid("source-bound reference roster count"));
        }
        for s in roster {
            let matches: Vec<_> = self
                .records
                .iter()
                .filter(|r| {
                    r.record.snapshot_id == s.id
                        && r.record.snapshot_sha256 == manifest.fixture_sha256[&s.id]
                })
                .collect();
            if matches.len() != 1 {
                return Err(invalid("source-bound reference fixture hash/ID"));
            }
            let r = matches[0];
            if r.usage.calls > b.reference.calls
                || r.usage.work.real_terms > b.reference.work.real_terms
                || r.usage.work.atom_k_terms > b.reference.work.atom_k_terms
                || u128::from(r.usage.elapsed_ns) > b.reference.elapsed().as_nanos()
            {
                return Err(invalid(
                    "reference family calls/work/elapsed were not resumed",
                ));
            }
        }
        Ok(())
    }
    /// Strict selection engine with origin checked before any oracle/candidate call.
    pub fn evaluate_candidate(
        &self,
        prepared: &Path,
        candidate: &CandidateReport,
        b: &mut CandidateBudget,
        certificates: &mut Vec<TranslationCertificate>,
    ) -> Result<CandidateOutcome, CalibrationFailure> {
        self.evaluate_candidate_with_observer(prepared, candidate, b, certificates, |_| Ok(()))
    }
    pub fn evaluate_candidate_with_observer<F>(
        &self,
        prepared: &Path,
        candidate: &CandidateReport,
        b: &mut CandidateBudget,
        certificates: &mut Vec<TranslationCertificate>,
        observer: F,
    ) -> Result<CandidateOutcome, CalibrationFailure>
    where
        F: FnMut(CandidateEvent<'_>) -> Result<(), CalibrationFailure>,
    {
        if self.role != FixtureRole::Selection {
            return Err(invalid("selection references required"));
        }
        self.validate(prepared, b)?;
        let original = b.external_owned_heap_bytes;
        b.external_owned_heap_bytes =
            heap_sum([original, self.prepared_manifest_sha256.capacity() as u64])?;
        let result = evaluate_prepared_candidate_with_observer(
            prepared,
            candidate,
            &self.records,
            b,
            certificates,
            observer,
        );
        b.external_owned_heap_bytes = original;
        result
    }
}
fn heap_sum(values: impl IntoIterator<Item = u64>) -> Result<u64, CalibrationFailure> {
    values.into_iter().try_fold(0u64, |sum, n| {
        sum.checked_add(n).ok_or_else(|| {
            CalibrationFailure::BudgetExceeded("retained heap accounting overflow".into())
        })
    })
}
fn generate_references(
    prepared: &Path,
    role: FixtureRole,
    b: &mut CandidateBudget,
) -> Result<PreparedReferences, CalibrationFailure> {
    let (p, snapshots, _) = validate_prepared(prepared)?;
    let epoch = sha256(&prepared.join("manifest.json"))?;
    let mut records: Vec<ReferenceOutcome> = Vec::new();
    let original = b.reference.owned_heap_baseline;
    let result = (|| {
        for s in snapshots.iter().filter(|s| s.role == role) {
            let retained = records
                .iter()
                .try_fold(0u64, |sum, r| heap_sum([sum, r.owned_heap_bytes()?]))?;
            b.reference.owned_heap_baseline = heap_sum([
                original,
                b.external_owned_heap_bytes,
                bytes(&snapshots),
                bytes(&records),
                epoch.capacity() as u64,
                heap_sum(
                    snapshots
                        .iter()
                        .filter(|other| other.id != s.id)
                        .map(snapshot_bytes),
                )?,
                retained,
            ])?;
            let r = converged_reference_with_budget(s, &p, &mut b.reference)?;
            let failure = r.failure.clone();
            records.push(r);
            if failure.is_some() {
                break;
            }
        }
        Ok(())
    })();
    b.reference.owned_heap_baseline = original;
    result?;
    Ok(PreparedReferences {
        prepared_manifest_sha256: epoch,
        role,
        records,
    })
}
/// Required beginning of a real selection campaign: project the whole grid and
/// confirmation BEFORE partial accuracy work, then generate current-source refs.
/// Task7 journals this wrapper and retains the SAME CandidateBudget thereafter.
pub fn begin_selection(
    prepared: &Path,
    b: &mut CandidateBudget,
) -> Result<PreparedReferences, CalibrationFailure> {
    begin_selection_with_generator(prepared, b, |prepared, b| {
        generate_references(prepared, FixtureRole::Selection, b)
    })
}
fn begin_selection_with_generator<F>(
    prepared: &Path,
    b: &mut CandidateBudget,
    generate: F,
) -> Result<PreparedReferences, CalibrationFailure>
where
    F: FnOnce(&Path, &mut CandidateBudget) -> Result<PreparedReferences, CalibrationFailure>,
{
    let (p, s, _) = validate_prepared(prepared)?;
    preflight_selection_sweep(&p, &s, b)?;
    generate(prepared, b)
}

#[derive(Clone, Debug, Default, serde::Serialize, serde::Deserialize)]
#[serde(deny_unknown_fields)]
pub struct TimingBudget {
    pub calls: u64,
}
impl TimingBudget {
    fn reserve(
        &mut self,
        n: u64,
        atoms: usize,
        p: &Protocol,
        b: &mut CandidateBudget,
    ) -> Result<(), CalibrationFailure> {
        b.reference.checkpoint()?;
        let total = self
            .calls
            .checked_add(n)
            .ok_or_else(|| CalibrationFailure::BudgetExceeded("timing call overflow".into()))?;
        if total > p.limits.max_timing_calls {
            return Err(CalibrationFailure::BudgetExceeded("timing call cap".into()));
        }
        let pairs = (atoms as u64)
            .checked_mul(atoms.saturating_sub(1) as u64)
            .and_then(|n| n.checked_div(2))
            .and_then(|pairs| pairs.checked_mul(n))
            .and_then(|pairs| pairs.checked_add(b.candidate_real_pair_upper_bound))
            .ok_or_else(|| {
                CalibrationFailure::BudgetExceeded("timing pair work overflow".into())
            })?;
        if pairs
            .checked_add(b.reference.work.real_terms)
            .is_none_or(|sum| sum > p.limits.max_real_terms_per_family)
        {
            return Err(CalibrationFailure::BudgetExceeded(
                "timing family real work cap".into(),
            ));
        }
        b.candidate_real_pair_upper_bound = pairs;
        self.calls = total;
        Ok(())
    }
}
#[derive(Clone, Debug, serde::Serialize)]
pub struct OrdinaryCostOutcome {
    pub cold_ns: Option<u64>,
    pub rebuild_ns: Option<u64>,
    pub warmups_completed: usize,
    pub batch_total_ns: Vec<u64>,
    pub timing: Option<TimingRecord>,
    pub failure: Option<CalibrationFailure>,
    pub timing_calls: u64,
    pub family_elapsed_ns: u64,
}
#[derive(Clone, Copy, Debug, PartialEq, serde::Serialize)]
pub enum CostStage {
    ColdStarted,
    RebuildStarted,
    WarmupStarted,
    BatchStarted,
    Cold,
    Rebuild,
    Warmup,
    Batch,
    Complete,
}
#[derive(serde::Serialize)]
pub struct CostEvent<'a> {
    pub stage: CostStage,
    pub outcome: &'a OrdinaryCostOutcome,
}
pub fn measure_ordinary_with_observer<F>(
    snapshot: &Snapshot,
    parameters: Parameters,
    p: &Protocol,
    timing: &mut TimingBudget,
    b: &mut CandidateBudget,
    mut observer: F,
) -> Result<OrdinaryCostOutcome, CalibrationFailure>
where
    F: FnMut(CostEvent<'_>) -> Result<(), CalibrationFailure>,
{
    if cfg!(debug_assertions) {
        return Err(invalid("ordinary campaign timing requires release build"));
    }
    if (p.warmups, p.batches, p.batch_calls) != (5, 5, 20) {
        return Err(invalid("frozen timing protocol"));
    }
    let mut out = OrdinaryCostOutcome {
        cold_ns: None,
        rebuild_ns: None,
        warmups_completed: 0,
        batch_total_ns: Vec::with_capacity(5),
        timing: None,
        failure: None,
        timing_calls: timing.calls,
        family_elapsed_ns: 0,
    };
    let result = (|| {
        timing.reserve(1, snapshot.charges.len(), p, b)?;
        cost_event(CostStage::ColdStarted, &mut out, timing, b, &mut observer)?;
        let started = std::time::Instant::now();
        let mut evaluator = SpmeEvaluator::with_protocol(snapshot, parameters, p)?;
        evaluator
            .compute_without_output_copy()
            .map_err(|e| invalid(e.to_string()))?;
        out.cold_ns = Some(
            u64::try_from(started.elapsed().as_nanos())
                .map_err(|_| CalibrationFailure::BudgetExceeded("cold timing overflow".into()))?,
        );
        cost_event(CostStage::Cold, &mut out, timing, b, &mut observer)?;
        timing.reserve(1, snapshot.charges.len(), p, b)?;
        cost_event(
            CostStage::RebuildStarted,
            &mut out,
            timing,
            b,
            &mut observer,
        )?;
        out.rebuild_ns = Some(evaluator.forced_rebuild_call()?);
        cost_event(CostStage::Rebuild, &mut out, timing, b, &mut observer)?;
        for _ in 0..5 {
            timing.reserve(1, snapshot.charges.len(), p, b)?;
            cost_event(CostStage::WarmupStarted, &mut out, timing, b, &mut observer)?;
            evaluator
                .compute_without_output_copy()
                .map_err(|e| invalid(e.to_string()))?;
            out.warmups_completed += 1;
            cost_event(CostStage::Warmup, &mut out, timing, b, &mut observer)?;
        }
        for _ in 0..5 {
            timing.reserve(20, snapshot.charges.len(), p, b)?;
            cost_event(CostStage::BatchStarted, &mut out, timing, b, &mut observer)?;
            let ns = evaluator.ordinary_batch(20)?;
            out.batch_total_ns.push(ns);
            // Durable observation comes before the next deadline/admission check.
            cost_event(CostStage::Batch, &mut out, timing, b, &mut observer)?;
            b.reference.checkpoint()?;
        }
        out.timing = Some(ordinary_timing_summary(
            out.cold_ns.ok_or_else(|| invalid("cold cost"))?,
            out.rebuild_ns.ok_or_else(|| invalid("rebuild cost"))?,
            out.batch_total_ns.clone(),
            p,
        )?);
        Ok(())
    })();
    if let Err(f) = result {
        out.failure = Some(f);
        out.timing = None;
    }
    if let Err(f) = cost_event(CostStage::Complete, &mut out, timing, b, &mut observer) {
        out.failure = Some(f);
        out.timing = None;
    }
    Ok(out)
}
fn cost_event<F>(
    stage: CostStage,
    out: &mut OrdinaryCostOutcome,
    timing: &TimingBudget,
    b: &CandidateBudget,
    observer: &mut F,
) -> Result<(), CalibrationFailure>
where
    F: FnMut(CostEvent<'_>) -> Result<(), CalibrationFailure>,
{
    out.timing_calls = timing.calls;
    out.family_elapsed_ns = u64::try_from(b.reference.elapsed().as_nanos())
        .map_err(|_| CalibrationFailure::BudgetExceeded("family elapsed overflow".into()))?;
    observer(CostEvent {
        stage,
        outcome: out,
    })
}
/// Convenience diagnostic. Journaling callers use the observer API to retain
/// partial costs even on errors or a later hard allocator exit.
pub fn measure_ordinary(
    snapshot: &Snapshot,
    parameters: Parameters,
    p: &Protocol,
    timing: &mut TimingBudget,
    b: &mut CandidateBudget,
) -> Result<TimingRecord, CalibrationFailure> {
    let out = measure_ordinary_with_observer(snapshot, parameters, p, timing, b, |_| Ok(()))?;
    if let Some(f) = out.failure {
        return Err(f);
    }
    out.timing
        .ok_or_else(|| invalid("incomplete ordinary timing"))
}
fn exact_grid(
    p: &Protocol,
    s: &[Snapshot],
    c: &[CandidateReport],
) -> Result<(), CalibrationFailure> {
    let first = s
        .iter()
        .find(|s| s.role == FixtureRole::Selection)
        .ok_or_else(|| invalid("selection roster"))?;
    let grid = candidate_grid(p, first.cell)?;
    if c.len() != grid.len() {
        return Err(invalid(
            "complete exact grid required; no measured shortlist",
        ));
    }
    for expected in grid {
        let matching: Vec<_> = c.iter().filter(|row| row.id == expected.id).collect();
        if matching.len() != 1 {
            return Err(invalid("missing/duplicate exact grid ID"));
        }
        let row = matching[0];
        if row.family != expected.family
            || row.parameters != expected.parameters
            || row.actual_mesh != expected.actual_mesh
            || row.aliases != expected.aliases
        {
            return Err(invalid("candidate differs from exact prepared grid"));
        }
    }
    Ok(())
}
/// All feasible candidates, never a mesh shortlist. A preflight cap failure
/// occurs before first timing; partial accepted costs survive any later error in
/// the caller's reports. Caller must not freeze on an Err.
pub fn measure_feasible_candidates(
    prepared: &Path,
    candidates: &mut [CandidateReport],
    timing: &mut TimingBudget,
    b: &mut CandidateBudget,
) -> Result<(), CalibrationFailure> {
    measure_feasible_candidates_with_observer(prepared, candidates, timing, b, |_, _| Ok(()))
}
pub fn measure_feasible_candidates_with_observer<F>(
    prepared: &Path,
    candidates: &mut [CandidateReport],
    timing: &mut TimingBudget,
    b: &mut CandidateBudget,
    mut observer: F,
) -> Result<(), CalibrationFailure>
where
    F: FnMut(usize, CostEvent<'_>) -> Result<(), CalibrationFailure>,
{
    let (p, s, _) = validate_prepared(prepared)?;
    exact_grid(&p, &s, candidates)?;
    let ids: Vec<_> = s
        .iter()
        .filter(|s| s.role == FixtureRole::Selection)
        .map(|s| s.id.as_str())
        .collect();
    let mut indices = Vec::new();
    for (i, c) in candidates.iter().enumerate() {
        if !feasible_targets(c, &ids, &p)?.is_empty() {
            indices.push(i);
        }
    }
    let projected = (indices.len() as u64)
        .checked_mul(107)
        .and_then(|n| timing.calls.checked_add(n))
        .ok_or_else(|| CalibrationFailure::BudgetExceeded("projected timing overflow".into()))?;
    if projected > p.limits.max_timing_calls {
        return Err(CalibrationFailure::BudgetExceeded(
            "cannot time every feasible candidate within cap".into(),
        ));
    }
    if candidates.iter().any(|c| c.timing.is_some()) {
        return Err(invalid("ordinary sweep requires untouched timing records"));
    }
    let snapshot = s
        .iter()
        .find(|s| s.role == FixtureRole::Selection)
        .ok_or_else(|| invalid("timing fixture"))?;
    for i in indices {
        let id = candidates[i].id;
        let out = measure_ordinary_with_observer(
            snapshot,
            candidates[i].parameters.clone(),
            &p,
            timing,
            b,
            |event| observer(id, event),
        )?;
        if let Some(f) = out.failure {
            return Err(f);
        }
        candidates[i].timing = out.timing;
    }
    Ok(())
}

fn rows_heap(rows: &Vec<MetricRow>) -> Result<u64, CalibrationFailure> {
    heap_sum([
        bytes(rows),
        heap_sum(
            rows.iter()
                .map(|r| r.snapshot_id.capacity() as u64 + bytes(&r.passing_targets)),
        )?,
    ])
}
fn report_heap(c: &CandidateReport) -> Result<u64, CalibrationFailure> {
    heap_sum([
        c.family.capacity() as u64,
        bytes(&c.aliases),
        rows_heap(&c.rows)?,
        c.timing.as_ref().map_or(0, |t| bytes(&t.batch_total_ns)),
    ])
}
fn decisions_heap(d: &Vec<TargetDecision>) -> Result<u64, CalibrationFailure> {
    heap_sum([
        bytes(d),
        d.iter().try_fold(0, |sum, d| {
            heap_sum([sum, rows_heap(&d.confirmation_rows)?])
        })?,
    ])
}
/// Private construction and no Deserialize: heldout compute can only consume
/// the all-target decision set frozen from the complete measured grid.
#[derive(serde::Serialize)]
pub struct FrozenSelection {
    protocol: Protocol,
    prepared_manifest_sha256: String,
    candidates: Vec<CandidateReport>,
    decisions: Vec<TargetDecision>,
}
pub fn freeze_prepared_selection(
    prepared: &Path,
    candidates: &[CandidateReport],
) -> Result<FrozenSelection, CalibrationFailure> {
    let (p, s, _) = validate_prepared(prepared)?;
    exact_grid(&p, &s, candidates)?;
    let ids: Vec<_> = s
        .iter()
        .filter(|s| s.role == FixtureRole::Selection)
        .map(|s| s.id.as_str())
        .collect();
    let decisions = p
        .targets
        .iter()
        .map(|target| select(candidates, &ids, *target, &p))
        .collect::<Result<Vec<_>, _>>()?;
    Ok(FrozenSelection {
        protocol: p,
        prepared_manifest_sha256: sha256(&prepared.join("manifest.json"))?,
        candidates: candidates.to_vec(),
        decisions,
    })
}
#[derive(serde::Serialize)]
pub struct ConfirmationOutcome {
    pub decisions: Vec<TargetDecision>,
    pub references: Option<PreparedReferences>,
    pub candidates: Vec<CandidateOutcome>,
    pub failure: Option<CalibrationFailure>,
}
#[derive(serde::Deserialize)]
#[serde(deny_unknown_fields)]
struct FrozenWire {
    protocol: Protocol,
    prepared_manifest_sha256: String,
    candidates: Vec<CandidateReport>,
    decisions: Vec<TargetDecision>,
}
impl FrozenSelection {
    pub fn write_new(&self, path: &Path) -> Result<String, CalibrationFailure> {
        super::write_new(path, &json_bytes(self)?)?;
        sha256(path)
    }
    /// Validate authoritative journal hash, exact epoch/grid and every frozen
    /// decision. This performs pure consistency checks, never heldout tuning.
    pub fn load(
        prepared: &Path,
        path: &Path,
        authoritative_phase_sha256: &str,
    ) -> Result<Self, CalibrationFailure> {
        if authoritative_phase_sha256.len() != 64 || sha256(path)? != authoritative_phase_sha256 {
            return Err(invalid(
                "frozen phase digest differs from authoritative journal",
            ));
        }
        let wire: FrozenWire =
            serde_json::from_slice(&std::fs::read(path).map_err(io)?).map_err(io)?;
        let frozen = freeze_prepared_selection(prepared, &wire.candidates)?;
        if frozen.protocol != wire.protocol
            || frozen.prepared_manifest_sha256 != wire.prepared_manifest_sha256
            || frozen.decisions != wire.decisions
        {
            return Err(invalid("frozen decision/source consistency mismatch"));
        }
        Ok(frozen)
    }
    pub fn decisions(&self) -> &[TargetDecision] {
        &self.decisions
    }
    pub fn confirm(
        &self,
        prepared: &Path,
        b: &mut CandidateBudget,
        certificates: &mut Vec<TranslationCertificate>,
    ) -> Result<ConfirmationOutcome, CalibrationFailure> {
        self.confirm_with_observer(prepared, b, certificates, |_| Ok(()))
    }
    pub fn confirm_with_observer<F>(
        &self,
        prepared: &Path,
        b: &mut CandidateBudget,
        certificates: &mut Vec<TranslationCertificate>,
        mut observer: F,
    ) -> Result<ConfirmationOutcome, CalibrationFailure>
    where
        F: FnMut(CandidateEvent<'_>) -> Result<(), CalibrationFailure>,
    {
        let (p, s, _) = validate_prepared(prepared)?;
        if self.prepared_manifest_sha256 != sha256(&prepared.join("manifest.json"))? {
            return Err(invalid("frozen selection origin mismatch"));
        }
        let original_external = b.external_owned_heap_bytes;
        let retained = heap_sum([
            original_external,
            candidates::protocol_bytes(&self.protocol),
            self.prepared_manifest_sha256.capacity() as u64,
            bytes(&self.candidates),
            self.candidates
                .iter()
                .try_fold(0, |sum, c| heap_sum([sum, report_heap(c)?]))?,
            decisions_heap(&self.decisions)?,
            candidates::protocol_bytes(&p),
            bytes(&s),
            heap_sum(s.iter().map(snapshot_bytes))?,
        ])?;
        let mut out = ConfirmationOutcome {
            decisions: self.decisions.clone(),
            references: None,
            candidates: Vec::new(),
            failure: None,
        };
        let mut selected: Vec<_> = self
            .decisions
            .iter()
            .filter_map(|d| d.selected_id)
            .collect();
        selected.sort_unstable();
        selected.dedup();
        if selected.is_empty() {
            return Ok(out);
        }
        // Every target ID is already frozen above, BEFORE any heldout reference or candidate compute.
        let result = (|| {
            // References generated here must also account for all retained
            // selection certificates; the row engine charges them itself later.
            b.external_owned_heap_bytes = heap_sum([
                retained,
                bytes(&selected),
                decisions_heap(&out.decisions)?,
                bytes(certificates),
                certificates
                    .iter()
                    .try_fold(0, |sum, c| heap_sum([sum, c.owned_heap_bytes()?]))?,
            ])?;
            out.references = Some(generate_references(prepared, FixtureRole::Confirmation, b)?);
            let refs = out
                .references
                .as_ref()
                .ok_or_else(|| invalid("confirmation references"))?;
            if let Some(f) = refs.failure() {
                return Err(f.clone());
            }
            refs.validate(prepared, b)?;
            let ids: Vec<_> = s
                .iter()
                .filter(|s| s.role == FixtureRole::Confirmation)
                .map(|s| s.id.as_str())
                .collect();
            let grid = candidate_grid(&p, s[0].cell)?;
            for id in selected.iter().copied() {
                let candidate = grid
                    .iter()
                    .find(|c| c.id == id)
                    .ok_or_else(|| invalid("frozen grid ID"))?;
                b.external_owned_heap_bytes = heap_sum([
                    retained,
                    bytes(&selected),
                    decisions_heap(&out.decisions)?,
                    bytes(&grid),
                    grid.iter()
                        .try_fold(0, |sum, c| heap_sum([sum, report_heap(c)?]))?,
                    bytes(&out.candidates),
                    out.candidates
                        .iter()
                        .try_fold(0, |sum, c| heap_sum([sum, c.owned_heap_bytes()?]))?,
                ])?;
                let evaluated = candidates::evaluate_prepared_role_with_observer(
                    prepared,
                    candidate,
                    &refs.records,
                    b,
                    certificates,
                    FixtureRole::Confirmation,
                    &mut observer,
                )?;
                let failure = evaluated.failure.clone();
                out.candidates.push(evaluated);
                let c = &out
                    .candidates
                    .last()
                    .ok_or_else(|| invalid("confirmation outcome"))?
                    .report;
                for d in out
                    .decisions
                    .iter_mut()
                    .filter(|d| d.selected_id == Some(id))
                {
                    if failure.is_some() {
                        d.confirmation_rows = c.rows.clone();
                        d.confirmation_status = Status::ConfirmationFailed;
                    } else {
                        finish_confirmation(d, c, &ids, &p)?;
                    }
                }
                if let Some(f) = failure {
                    return Err(f);
                }
            }
            Ok(())
        })();
        b.external_owned_heap_bytes = original_external;
        if let Err(f) = result {
            // A resource/reference/I/O failure prevents every new fragment. IDs stay frozen.
            for d in out.decisions.iter_mut().filter(|d| d.selected_id.is_some()) {
                d.confirmation_status = Status::ConfirmationFailed;
            }
            out.failure = Some(f);
        }
        Ok(out)
    }
    pub fn fragment(
        &self,
        out: &ConfirmationOutcome,
        target: f64,
        p: &Protocol,
    ) -> Result<Option<String>, CalibrationFailure> {
        if out.failure.is_some() {
            return Ok(None);
        }
        if p != &self.protocol {
            return Err(invalid("fragment protocol differs from frozen selection"));
        }
        let d = out
            .decisions
            .iter()
            .find(|d| d.target == target)
            .ok_or_else(|| invalid("confirmation target"))?;
        let original = self
            .decisions
            .iter()
            .find(|d| d.target == target)
            .ok_or_else(|| invalid("frozen target"))?;
        if d.selected_id != original.selected_id {
            return Err(invalid("confirmation changed frozen ID"));
        }
        let Some(id) = d.selected_id else {
            return Ok(None);
        };
        if d.confirmation_status != Status::Passed {
            return Ok(None);
        }
        let confirmation = out
            .candidates
            .iter()
            .find(|c| c.report.id == id)
            .ok_or_else(|| invalid("confirmed candidate outcome missing"))?;
        if confirmation.failure.is_some() || confirmation.report.rows != d.confirmation_rows {
            return Err(invalid("confirmation decision/raw outcome mismatch"));
        }
        let snapshots = materialize(p)?;
        let ids: Vec<_> = snapshots
            .iter()
            .filter(|s| s.role == FixtureRole::Confirmation)
            .map(|s| s.id.as_str())
            .collect();
        if !feasible_targets(&confirmation.report, &ids, p)?.contains(&target) {
            return Err(invalid("complete confirmation gates missing"));
        }
        let c = self
            .candidates
            .iter()
            .find(|c| c.id == id)
            .ok_or_else(|| invalid("frozen candidate"))?;
        selected_fragment(d, c, p)
    }
}

/// One registered three-ion functional probe. Diagnostic ordinary costs are
/// retained separately even when this single candidate is infeasible; no ranking,
/// heldout tuning, selected fragment or full-campaign claim is made.
pub fn cost_smoke(prepared: &Path, out: &Path) -> Result<(), CalibrationFailure> {
    let (p, snapshots, _) = validate_prepared(prepared)?;
    if p.family != "special-3" || p.threads != 1 || p.reference_ladder.len() != 6 {
        return Err(invalid(
            "cost-smoke requires six-level special-3, one-thread prepared protocol",
        ));
    }
    std::fs::create_dir(out).map_err(io)?;
    let mut b = CandidateBudget::new(&p.limits, ReferenceBudget::new(&p.limits));
    b.external_owned_heap_bytes = heap_sum([
        candidates::protocol_bytes(&p),
        bytes(&snapshots),
        heap_sum(snapshots.iter().map(snapshot_bytes))?,
    ])?;
    let refs = begin_selection(prepared, &mut b)?;
    let reference_path = out.join("selection-reference.json");
    let digest = refs.write_new(&reference_path)?;
    if let Some(f) = refs.failure() {
        return Err(f.clone());
    }
    // Round-trip through the explicit trusted phase hash. The hash was obtained
    // from the newly written current-source phase, never from an imported file.
    let refs = PreparedReferences::load_selection(prepared, &reference_path, &digest, &mut b)?;
    let fixture = snapshots
        .iter()
        .find(|s| s.id == "special-3-selection-charged-scaled")
        .ok_or_else(|| invalid("registered smoke fixture"))?;
    let grid = candidate_grid(&p, fixture.cell)?;
    let candidate = grid
        .iter()
        .find(|c| {
            c.parameters.kappa == 1.6 / p.a
                && c.parameters.real_cutoff
                    == 0.40 * fixture.cell.into_iter().fold(f64::INFINITY, f64::min)
                && c.parameters.kspace_cutoff == 24. / p.a
        })
        .ok_or_else(|| invalid("registered smoke candidate"))?;
    b.external_owned_heap_bytes = heap_sum([
        b.external_owned_heap_bytes,
        bytes(&grid),
        grid.iter()
            .try_fold(0, |sum, c| heap_sum([sum, report_heap(c)?]))?,
    ])?;
    let mut certificates = Vec::new();
    let outcome = refs.evaluate_candidate(prepared, candidate, &mut b, &mut certificates)?;
    super::write_new(&out.join("candidate.json"), &json_bytes(&outcome)?)?;
    super::write_new(
        &out.join("translation-certificates.json"),
        &json_bytes(&certificates)?,
    )?;
    if let Some(f) = outcome.failure {
        return Err(f);
    }
    let mut timing = TimingBudget::default();
    let costs = measure_ordinary(
        fixture,
        candidate.parameters.clone(),
        &p,
        &mut timing,
        &mut b,
    )?;
    #[derive(serde::Serialize)]
    struct Smoke<'a> {
        label: &'a str,
        prepared_manifest_sha256: String,
        candidate_id: usize,
        timing: TimingRecord,
        timing_calls: u64,
        usage: candidates::CandidateUsage,
    }
    let record = Smoke {
        label: "bounded functional ordinary-cost diagnostic; no ranking or selected configuration",
        prepared_manifest_sha256: sha256(&prepared.join("manifest.json"))?,
        candidate_id: candidate.id,
        timing: costs,
        timing_calls: timing.calls,
        usage: b.usage(),
    };
    super::write_new(&out.join("ordinary-diagnostic.json"), &json_bytes(&record)?)?;
    Ok(())
}

#[cfg(test)]
mod projection_admission_tests {
    use super::*;
    #[test]
    fn complete_projection_rejects_before_begin_selection_generates_references() {
        let mut p = frozen_protocol("special", 2, 1, true).unwrap();
        // The old selection-certificate-only bound is45,307,584atom-k visits;
        // complete reference ladders+confirmation require68,684,928.
        p.limits.max_atom_k_terms_per_family = 50_000_000;
        let dir =
            std::env::temp_dir().join(format!("b05-task5-fix1-projection-{}", std::process::id()));
        if dir.exists() {
            std::fs::remove_dir_all(&dir).unwrap();
        }
        prepare(&dir, p.clone()).unwrap();
        let mut b = CandidateBudget::new(&p.limits, ReferenceBudget::new(&p.limits));
        let mut generated = 0;
        // Pure admission seam: replacing reference generation prevents the RED
        // from doing physics when the incomplete preflight incorrectly admits it.
        let result = begin_selection_with_generator(&dir, &mut b, |_, _| {
            generated += 1;
            Ok(PreparedReferences {
                prepared_manifest_sha256: String::new(),
                role: FixtureRole::Selection,
                records: vec![],
            })
        });
        assert!(matches!(result, Err(CalibrationFailure::BudgetExceeded(_))));
        assert_eq!(generated, 0);
        assert_eq!(b.reference.calls, 0);
        assert_eq!(b.accuracy_calls, 0);
        assert_eq!(
            b.reference.work,
            WorkCounts {
                real_terms: 0,
                atom_k_terms: 0
            }
        );
        assert_eq!(b.candidate_real_pair_upper_bound, 0);
        std::fs::remove_dir_all(dir).unwrap();
    }
}

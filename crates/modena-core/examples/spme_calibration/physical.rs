//! Source-bound runtime oracle physical controls. No target/cost certification.
use super::candidates::{bytes, protocol_bytes, snapshot_bytes, CandidateUsage};
use super::reference::{planned_work, ReferenceOrigin};
use super::*;
use std::mem::size_of;
#[path = "physical_controls.rs"]
mod controls;
const A: f64 = 2.2;
fn cap(s: &str) -> CalibrationFailure {
    CalibrationFailure::BudgetExceeded(s.into())
}
fn add(a: u64, b: u64) -> Result<u64, CalibrationFailure> {
    a.checked_add(b)
        .ok_or_else(|| cap("physical accounting overflow"))
}
fn mul(a: u64, b: u64) -> Result<u64, CalibrationFailure> {
    a.checked_mul(b)
        .ok_or_else(|| cap("physical accounting overflow"))
}
fn failed(s: impl Into<String>) -> CalibrationFailure {
    CalibrationFailure::ReferencePrecisionLimited(s.into())
}
#[derive(Clone, Debug, PartialEq, serde::Serialize)]
pub struct Membership {
    pub real: Vec<(usize, usize, [i32; 3])>,
    pub reciprocal: Vec<[i32; 3]>,
}
#[derive(Clone, Debug, serde::Serialize)]
pub struct PhysicalRecord {
    pub label: String,
    pub origin: ReferenceOrigin,
    pub snapshot: Snapshot,
    pub settings: ReferenceSettings,
    pub evaluation: Option<Evaluation>,
    pub audit: Option<ReferenceAudit>,
    pub membership: Option<Membership>,
    pub usage: CandidateUsage,
}
#[derive(Clone, Debug, serde::Serialize)]
pub struct PhysicalCheck {
    pub label: String,
    pub actual: f64,
    pub expected: f64,
    pub tolerance: f64,
    pub passed: bool,
}
/// Private construction, no Deserialize: a status cannot mint a physical proof.
#[derive(serde::Serialize)]
pub struct PhysicalOutcome {
    manifest_sha256: String,
    protocol_sha256: String,
    specification_sha256: String,
    source_sha256: std::collections::BTreeMap<String, String>,
    family: String,
    threads: usize,
    records: Vec<PhysicalRecord>,
    checks: Vec<PhysicalCheck>,
    failure: Option<CalibrationFailure>,
    finalization_failure: Option<CalibrationFailure>,
    usage: CandidateUsage,
    auxiliary_calls: u64,
    completed: bool,
}
impl PhysicalOutcome {
    pub fn records(&self) -> &[PhysicalRecord] {
        &self.records
    }
    pub fn checks(&self) -> &[PhysicalCheck] {
        &self.checks
    }
    pub fn failure(&self) -> Option<&CalibrationFailure> {
        self.failure.as_ref()
    }
    pub fn finalization_failure(&self) -> Option<&CalibrationFailure> {
        self.finalization_failure.as_ref()
    }
    pub fn usage(&self) -> &CandidateUsage {
        &self.usage
    }
    pub fn passed(&self) -> bool {
        self.completed
            && self.failure.is_none()
            && self.finalization_failure.is_none()
            && self.checks.iter().all(|c| c.passed)
    }
    /// Cross-family oracle proof checks sources and thread identity; it does not
    /// transfer/reset the special family's continuing usage into an ionic family.
    pub fn validate_for(&self, prepared: &Path) -> Result<(), CalibrationFailure> {
        let (p, _, m) = validate_prepared(prepared)?;
        if !self.passed()
            || self.family != "special-3"
            || self.threads != p.threads
            || self.source_sha256 != m.source_sha256
            || self.specification_sha256 != m.specification_sha256
        {
            return Err(invalid("physical proof source/thread/status mismatch"));
        }
        Ok(())
    }
    /// Restore actual usage through ReferenceBudget::resume and CandidateBudget
    /// fields before reusing this proof in the later special-family phase.
    pub fn validate_special_resume(
        &self,
        prepared: &Path,
        b: &CandidateBudget,
    ) -> Result<(), CalibrationFailure> {
        self.validate_for(prepared)?;
        let used = b.usage();
        if self.manifest_sha256 != sha256(&prepared.join("manifest.json"))?
            || self.usage.reference.calls > used.reference.calls
            || self.usage.reference.work.real_terms > used.reference.work.real_terms
            || self.usage.reference.work.atom_k_terms > used.reference.work.atom_k_terms
            || self.usage.reference.elapsed_ns > used.reference.elapsed_ns
            || self.usage.accuracy_calls > used.accuracy_calls
            || self.usage.candidate_real_pair_upper_bound > used.candidate_real_pair_upper_bound
        {
            return Err(invalid(
                "physical prerequisite usage was not resumed in its own family",
            ));
        }
        Ok(())
    }
    pub fn owned_heap_bytes(&self) -> Result<u64, CalibrationFailure> {
        let mut n = add(bytes(&self.records), bytes(&self.checks))?;
        n = add(
            n,
            (self.manifest_sha256.capacity()
                + self.protocol_sha256.capacity()
                + self.specification_sha256.capacity()
                + self.family.capacity()) as u64,
        )?;
        // BTreeMap node capacity is opaque; charge a conservative node allowance.
        for (k, v) in &self.source_sha256 {
            n = add(n, add((k.capacity() + v.capacity()) as u64, 256)?)?;
        }
        for r in &self.records {
            n = add(n, record_heap(r)?)?;
        }
        for c in &self.checks {
            n = add(n, c.label.capacity() as u64)?;
        }
        n = add(n, failure_heap(self.failure.as_ref()))?;
        add(n, failure_heap(self.finalization_failure.as_ref()))
    }
}
fn failure_heap(failure: Option<&CalibrationFailure>) -> u64 {
    match failure {
        Some(
            CalibrationFailure::InvalidInput(s)
            | CalibrationFailure::UnsupportedDomain(s)
            | CalibrationFailure::BudgetExceeded(s)
            | CalibrationFailure::ReferenceUnconverged(s)
            | CalibrationFailure::ReferenceStagnated(s)
            | CalibrationFailure::ReferencePrecisionLimited(s)
            | CalibrationFailure::RelativeUndefined(s)
            | CalibrationFailure::Io(s),
        ) => s.capacity() as u64,
        _ => 0,
    }
}
fn retain_finalization_failure(
    primary: &mut Option<CalibrationFailure>,
    finalization: &mut Option<CalibrationFailure>,
    error: CalibrationFailure,
) {
    if primary.is_none() {
        *primary = Some(error.clone());
    }
    *finalization = Some(error);
}
fn origin_heap(o: &ReferenceOrigin) -> u64 {
    (o.manifest_sha256.capacity()
        + o.protocol_sha256.capacity()
        + o.snapshot_sha256.capacity()
        + o.snapshot_id.capacity()) as u64
}
fn record_heap(r: &PhysicalRecord) -> Result<u64, CalibrationFailure> {
    let mut n = add(
        snapshot_bytes(&r.snapshot),
        add(r.label.capacity() as u64, origin_heap(&r.origin))?,
    )?;
    if let Some(e) = &r.evaluation {
        n = add(n, bytes(&e.forces))?;
    }
    if let Some(m) = &r.membership {
        n = add(n, add(bytes(&m.real), bytes(&m.reciprocal))?)?;
    }
    Ok(n)
}
#[derive(serde::Serialize)]
pub enum PhysicalEvent<'a> {
    Admitted {
        origin: &'a ReferenceOrigin,
        label: &'a str,
        reserved: CandidateUsage,
    },
    AuxiliaryCompleted {
        record: &'a PhysicalRecord,
    },
    Completed {
        record: &'a PhysicalRecord,
    },
    Checked {
        check: &'a PhysicalCheck,
        usage: CandidateUsage,
    },
    Complete {
        outcome: &'a PhysicalOutcome,
    },
}
struct Runner<'a, F> {
    p: Protocol,
    budget: &'a mut CandidateBudget,
    out: PhysicalOutcome,
    observer: F,
    external: u64,
    retained: u64,
}
impl<F: FnMut(PhysicalEvent<'_>) -> Result<(), CalibrationFailure>> Runner<'_, F> {
    fn reserve_record(&mut self) -> Result<(), CalibrationFailure> {
        self.heap(size_of::<PhysicalRecord>() as u64)?;
        self.out.records.try_reserve_exact(1).map_err(io)
    }
    fn heap(&mut self, extra: u64) -> Result<(), CalibrationFailure> {
        let used = add(
            add(self.external, self.retained)?,
            add(self.out.owned_heap_bytes()?, extra)?,
        )?;
        if used > self.p.limits.max_owned_heap_bytes {
            return Err(cap("physical retained capacity cap"));
        }
        self.budget.reference.owned_heap_baseline = used;
        Ok(())
    }
    fn admit(
        &mut self,
        origin: &ReferenceOrigin,
        label: &str,
        w: &WorkCounts,
        reference_calls: u64,
        core: bool,
    ) -> Result<(), CalibrationFailure> {
        self.budget.reference.checkpoint()?;
        check_work(
            [1; 3],
            w,
            self.budget.reference.owned_heap_baseline,
            &self.p.limits,
        )?;
        let mut reserved = self.budget.usage();
        reserved.reference.calls = add(reserved.reference.calls, reference_calls)?;
        reserved.reference.work.real_terms = add(reserved.reference.work.real_terms, w.real_terms)?;
        reserved.reference.work.atom_k_terms =
            add(reserved.reference.work.atom_k_terms, w.atom_k_terms)?;
        if core {
            reserved.accuracy_calls = add(reserved.accuracy_calls, 1)?;
            reserved.candidate_real_pair_upper_bound =
                add(reserved.candidate_real_pair_upper_bound, 3)?;
        }
        if reserved.reference.calls > self.p.limits.max_accuracy_calls
            || reserved.accuracy_calls > self.p.limits.max_accuracy_calls
            || add(
                reserved.reference.work.real_terms,
                reserved.candidate_real_pair_upper_bound,
            )? > self.p.limits.max_real_terms_per_family
            || reserved.reference.work.atom_k_terms > self.p.limits.max_atom_k_terms_per_family
        {
            return Err(cap("physical continuing call/work cap"));
        }
        (self.observer)(PhysicalEvent::Admitted {
            origin,
            label,
            reserved,
        })
    }
    fn store(&mut self, record: PhysicalRecord, aux: bool) -> Result<usize, CalibrationFailure> {
        let index = self.out.records.len();
        self.out.records.push(record);
        self.heap(0)?;
        let r = &self.out.records[index];
        if aux {
            (self.observer)(PhysicalEvent::AuxiliaryCompleted { record: r })?;
        } else {
            (self.observer)(PhysicalEvent::Completed { record: r })?;
        }
        self.budget.reference.checkpoint()?;
        Ok(index)
    }
    fn oracle(
        &mut self,
        o: &ReferenceOrigin,
        s: &Snapshot,
        p: &ReferenceSettings,
        label: &str,
    ) -> Result<usize, CalibrationFailure> {
        self.heap(add(snapshot_bytes(s), origin_heap(o))?)?;
        self.reserve_record()?;
        let w = planned_work(s, p)?;
        self.admit(o, label, &w, 1, false)?;
        let (e, a) =
            reference_with_budget(s, p.clone(), &self.p.limits, &mut self.budget.reference)?;
        self.store(
            PhysicalRecord {
                label: label.into(),
                origin: o.clone(),
                snapshot: s.clone(),
                settings: p.clone(),
                evaluation: Some(e),
                audit: Some(a),
                membership: None,
                usage: self.budget.usage(),
            },
            false,
        )
    }
    fn aux_admit(
        &mut self,
        o: &ReferenceOrigin,
        s: &Snapshot,
        p: &ReferenceSettings,
        label: &str,
        half: bool,
    ) -> Result<WorkCounts, CalibrationFailure> {
        self.reserve_record()?;
        let mut w = planned_work(s, p)?;
        if half {
            w.real_terms = 0;
            w.atom_k_terms = mul(w.atom_k_terms, s.positions.len() as u64)?;
        } else {
            w.real_terms = add(w.real_terms, s.positions.len() as u64)?;
        }
        self.heap(add(
            snapshot_bytes(s),
            if half {
                mul(s.positions.len() as u64, 24)?
            } else {
                add(
                    mul(w.real_terms, size_of::<(usize, usize, [i32; 3])>() as u64)?,
                    mul(w.atom_k_terms, size_of::<[i32; 3]>() as u64)?,
                )?
            },
        )?)?;
        self.admit(o, label, &w, 1, false)?;
        // Auxiliary loops conservatively charge their complete rectangular bound,
        // including excluded tuples. Oracle calls retain their actual counters.
        self.budget.reference.work.real_terms =
            add(self.budget.reference.work.real_terms, w.real_terms)?;
        self.budget.reference.work.atom_k_terms =
            add(self.budget.reference.work.atom_k_terms, w.atom_k_terms)?;
        self.budget.reference.calls = add(self.budget.reference.calls, 1)?;
        self.out.auxiliary_calls = add(self.out.auxiliary_calls, 1)?;
        Ok(w)
    }
    fn value(&self, id: usize) -> Result<&Evaluation, CalibrationFailure> {
        self.out
            .records
            .get(id)
            .and_then(|r| r.evaluation.as_ref())
            .ok_or_else(|| invalid("physical evaluation index"))
    }
    fn check(
        &mut self,
        label: impl Into<String>,
        actual: f64,
        expected: f64,
        tolerance: f64,
    ) -> Result<(), CalibrationFailure> {
        if !actual.is_finite() || !expected.is_finite() || !tolerance.is_finite() {
            return Err(failed("nonfinite physical comparison"));
        }
        let passed = (actual - expected).abs() <= tolerance;
        let label = label.into();
        self.heap(add(
            size_of::<PhysicalCheck>() as u64,
            label.capacity() as u64,
        )?)?;
        self.out.checks.try_reserve_exact(1).map_err(io)?;
        self.out.checks.push(PhysicalCheck {
            label,
            actual,
            expected,
            tolerance,
            passed,
        });
        self.heap(0)?;
        let check = self
            .out
            .checks
            .last()
            .ok_or_else(|| invalid("physical check"))?;
        (self.observer)(PhysicalEvent::Checked {
            check,
            usage: self.budget.usage(),
        })?;
        if !passed {
            return Err(failed(check.label.clone()));
        }
        Ok(())
    }
    fn compare(
        &mut self,
        label: &str,
        a: &Evaluation,
        b: &Evaluation,
        tol: f64,
    ) -> Result<(), CalibrationFailure> {
        self.check(format!("{label}: energy"), a.energy, b.energy, tol / A)?;
        self.check(format!("{label}: W"), a.strain_w, b.strain_w, tol / A)?;
        if a.forces.len() != b.forces.len() {
            return Err(invalid("physical force shape"));
        }
        for (i, (x, y)) in a.forces.iter().zip(&b.forces).enumerate() {
            for axis in 0..3 {
                self.check(
                    format!("{label}: force {i}/{axis}"),
                    x[axis],
                    y[axis],
                    tol / A.powi(2),
                )?;
            }
        }
        Ok(())
    }
    fn strict(
        &mut self,
        label: &str,
        left: &Evaluation,
        right: &Evaluation,
        final_value: &Evaluation,
    ) -> Result<(), CalibrationFailure> {
        let m = metrics(right, left, left.forces.len(), 1e-8 / A.powi(2))?;
        let d = metrics(
            final_value,
            final_value,
            final_value.forces.len(),
            1e-8 / A.powi(2),
        )?
        .reference_force_rms;
        let absolute = 1e-8 / A.powi(2);
        let rms = if d > absolute {
            absolute.min(0.01 * 1e-4 * d)
        } else {
            absolute
        };
        self.check(format!("{label}: RMS"), m.force_rms_abs, 0., rms)?;
        self.check(
            format!("{label}: max atom"),
            m.max_atom_force_error,
            0.,
            1e-7 / A.powi(2),
        )?;
        self.check(format!("{label}: E/N"), m.energy_per_atom_abs, 0., 1e-9 / A)?;
        self.check(format!("{label}: W/N"), m.w_per_atom_abs, 0., 1e-9 / A)
    }
}
fn settings() -> ReferenceSettings {
    ReferenceSettings {
        kappa: 1.1 / A,
        real_cutoff: 8.73,
        kspace_cutoff: 3.75,
    }
}
/// Preflight exact registered call roster with conservative two-pass bounds.
fn special_plan(
    p: &Protocol,
    ss: &[Snapshot],
    b: &CandidateBudget,
) -> Result<(), CalibrationFailure> {
    let mut work = b.reference.work.clone();
    let mut calls = b.reference.calls;
    let mut model = b.accuracy_calls;
    let mut heap = 0;
    let mut core_scratch = 0;
    let mut project = |s: &Snapshot,
                       r: ReferenceSettings,
                       count: u64,
                       aux: bool,
                       half: bool|
     -> Result<(), CalibrationFailure> {
        let mut w = planned_work(s, &r)?;
        if half {
            w.real_terms = 0;
            w.atom_k_terms = mul(w.atom_k_terms, s.positions.len() as u64)?;
        } else if aux {
            w.real_terms = add(w.real_terms, s.positions.len() as u64)?;
        }
        check_work([1; 3], &w, 0, &p.limits)?;
        work.real_terms = add(work.real_terms, mul(w.real_terms, count)?)?;
        work.atom_k_terms = add(work.atom_k_terms, mul(w.atom_k_terms, count)?)?;
        calls = add(calls, count)?;
        let payload = if aux && !half {
            add(mul(w.real_terms, 32)?, mul(w.atom_k_terms, 12)?)?
        } else {
            mul(s.positions.len() as u64, 24)?
        };
        // Per-record structs, copies, origins/labels; checks <=32 per oracle.
        heap = add(heap, mul(count, add(payload, 8192)?)?)?;
        Ok(())
    };
    for s in ss {
        // Strain lambda=1+h has a possibly larger rectangular bound.
        let mut strained = s.clone();
        strained.cell = strained.cell.map(|x| x * 1.00001);
        project(&strained, settings(), 42, false, false)?; // force37 + strain5
        project(&strained, settings(), 42, true, false)?;
        project(s, settings(), 1, false, false)?; // full k
        project(
            s,
            ReferenceSettings {
                kspace_cutoff: 0.1,
                ..settings()
            },
            1,
            false,
            false,
        )?;
        project(s, settings(), 1, true, true)?;
        project(s, settings(), 1, false, false)?; // background
        project(s, settings(), 8 + 19, false, false)?; // translation/wrap incl bases
        if s.scales.is_empty() {
            project(
                s,
                ReferenceSettings {
                    kappa: 0.5,
                    real_cutoff: 3.7,
                    kspace_cutoff: 3.71,
                },
                9,
                false,
                false,
            )?;
        }
        for (r, k) in [(8.8, 22.), (11., 24.), (13.2, 26.)] {
            project(
                s,
                ReferenceSettings {
                    kappa: 2.,
                    real_cutoff: r,
                    kspace_cutoff: k,
                },
                1,
                false,
                false,
            )?;
        }
        project(
            s,
            ReferenceSettings {
                kappa: 2.,
                real_cutoff: 13.2,
                kspace_cutoff: 26.,
            },
            2,
            false,
            false,
        )?;
        project(
            s,
            ReferenceSettings {
                kappa: 2.,
                real_cutoff: 4.41,
                kspace_cutoff: 26.,
            },
            1,
            false,
            false,
        )?;
        // Core's bounded direct reciprocal work charged conservatively as an
        // auxiliary two-pass*n bound. No oracle call is charged for this core call.
        project(
            s,
            ReferenceSettings {
                kappa: 2.,
                real_cutoff: 4.41,
                kspace_cutoff: 26.,
            },
            1,
            true,
            true,
        )?;
        let core_work = planned_work(
            s,
            &ReferenceSettings {
                kappa: 2.,
                real_cutoff: 4.41,
                kspace_cutoff: 26.,
            },
        )?;
        core_scratch = core_scratch.max(add(
            mul(core_work.atom_k_terms / (2 * s.positions.len() as u64), 64)?,
            65536,
        )?);
        model = add(model, 1)?;
    }
    let mut one = ss
        .first()
        .ok_or_else(|| invalid("special control roster"))?
        .clone();
    one.cell = [A; 3];
    one.positions = vec![[0.4, 0.5, 0.6]];
    one.charges = vec![1.];
    one.scales.clear();
    project(
        &one,
        ReferenceSettings {
            kappa: 0.5,
            real_cutoff: 2.42,
            kspace_cutoff: 0.1,
        },
        1,
        false,
        false,
    )?;
    if calls > p.limits.max_accuracy_calls
        || model > p.limits.max_accuracy_calls
        || add(
            work.real_terms,
            add(b.candidate_real_pair_upper_bound, mul(ss.len() as u64, 3)?)?,
        )? > p.limits.max_real_terms_per_family
        || work.atom_k_terms > p.limits.max_atom_k_terms_per_family
    {
        return Err(cap("complete physical control work/call cap"));
    }
    if add(
        add(b.reference.owned_heap_baseline, b.external_owned_heap_bytes)?,
        add(heap, core_scratch)?,
    )? > p.limits.max_owned_heap_bytes
    {
        return Err(cap("complete physical control heap cap"));
    }
    Ok(())
}
pub fn verify_special_with_observer<F>(
    prepared: &Path,
    b: &mut CandidateBudget,
    observer: F,
) -> Result<PhysicalOutcome, CalibrationFailure>
where
    F: FnMut(PhysicalEvent<'_>) -> Result<(), CalibrationFailure>,
{
    let (p, snapshots, m) = validate_prepared(prepared)?;
    if !b.matches_limits(&p.limits) {
        return Err(invalid("physical continuing budget limits mismatch"));
    }
    if p.family != "special-3" {
        return Err(invalid(
            "special physical proof requires special prepared family",
        ));
    }
    b.reference.checkpoint()?;
    let ss: Vec<_> = snapshots
        .iter()
        .filter(|s| {
            s.role == FixtureRole::Selection
                && (s.id.ends_with("-bare") || s.id.ends_with("-scaled"))
        })
        .cloned()
        .collect();
    special_plan(&p, &ss, b)?;
    let original = b.reference.owned_heap_baseline;
    let external = add(original, b.external_owned_heap_bytes)?;
    let snapshot_heap = snapshots
        .iter()
        .try_fold(bytes(&snapshots), |n, s| add(n, snapshot_bytes(s)))?;
    let special_heap = ss
        .iter()
        .try_fold(bytes(&ss), |n, s| add(n, snapshot_bytes(s)))?;
    // Registered3atom recipes have at most two perturbed snapshots, eight
    // copied evaluations and short labels outside retained records at one time.
    // Reserve their conservative8KiB scratch bound; never global allocator totals.
    let retained = add(
        8192,
        add(protocol_bytes(&p), add(snapshot_heap, special_heap)?)?,
    )?;
    let out = PhysicalOutcome {
        manifest_sha256: sha256(&prepared.join("manifest.json"))?,
        protocol_sha256: m.protocol_sha256,
        specification_sha256: m.specification_sha256,
        source_sha256: m.source_sha256,
        family: p.family.clone(),
        threads: p.threads,
        records: vec![],
        checks: vec![],
        failure: None,
        finalization_failure: None,
        usage: b.usage(),
        auxiliary_calls: 0,
        completed: false,
    };
    let mut r = Runner {
        p,
        budget: b,
        out,
        observer,
        external,
        retained,
    };
    let result = (|| {
        for s in &ss {
            let o = ReferenceOrigin::prepared(prepared, &s.id)?;
            controls::special(&mut r, &o, s)?;
        }
        let s = ss.first().ok_or_else(|| invalid("special roster"))?;
        controls::self_images(&mut r, &ReferenceOrigin::prepared(prepared, &s.id)?, s)?;
        validate_prepared(prepared)?; // refuse source edits during work
        Ok(())
    })();
    r.out.failure = result.err();
    r.out.completed = r.out.failure.is_none();
    r.out.usage = r.budget.usage();
    r.budget.reference.owned_heap_baseline = original;
    if let Err(e) = (r.observer)(PhysicalEvent::Complete { outcome: &r.out }) {
        retain_finalization_failure(&mut r.out.failure, &mut r.out.finalization_failure, e);
        r.out.completed = false;
    }
    Ok(r.out)
}

/// Absolute crystal sanity only. RelativeUndefined is valid after strict ladder
/// convergence; no relative target decision is represented by this type.
#[derive(serde::Serialize)]
pub struct MadelungOutcome {
    manifest_sha256: String,
    source_sha256: std::collections::BTreeMap<String, String>,
    specification_sha256: String,
    threads: usize,
    replicas: u32,
    reference: ReferenceOutcome,
    sanity: Option<PhysicalCheck>,
    failure: Option<CalibrationFailure>,
    finalization_failure: Option<CalibrationFailure>,
    usage: CandidateUsage,
}
impl MadelungOutcome {
    pub fn reference(&self) -> &ReferenceOutcome {
        &self.reference
    }
    pub fn sanity(&self) -> Option<&PhysicalCheck> {
        self.sanity.as_ref()
    }
    pub fn failure(&self) -> Option<&CalibrationFailure> {
        self.failure.as_ref()
    }
    pub fn finalization_failure(&self) -> Option<&CalibrationFailure> {
        self.finalization_failure.as_ref()
    }
    pub fn usage(&self) -> &CandidateUsage {
        &self.usage
    }
    pub fn passed(&self) -> bool {
        self.validate_success().is_ok()
    }
    /// Shared proof invariant; future checked persistence must validate this AND
    /// the authoritative source/phase chain. A successful scalar sanity alone is
    /// never a substitute for the strict absolute reference ladder.
    pub fn validate_success(&self) -> Result<(), CalibrationFailure> {
        let r = &self.reference;
        let levels = &r.record.levels;
        let sanity = self
            .sanity
            .as_ref()
            .ok_or_else(|| invalid("missing Madelung sanity"))?;
        let n = match self.replicas {
            2 => 64,
            3 => 216,
            _ => return Err(invalid("Madelung replicas")),
        };
        if self.failure.is_some()
            || self.finalization_failure.is_some()
            || r.failure.is_some()
            || r.finalization_failure.is_some()
            || r.failed_attempt.is_some()
            || !matches!(r.record.status, Status::Passed | Status::RelativeUndefined)
            || !matches!(levels.len(), 4 | 6)
            || r.record.chosen_level != levels.len().checked_sub(1)
            || r.record.crosschecks.len() != 2
            || r.audits.len() != levels.len() + 2
            || r.refinement_scores.len() + 1 != levels.len()
            || r.crosscheck_scores.len() != 2
            || r.refinement_scores
                .iter()
                .chain(&r.crosscheck_scores)
                .any(|x| !x.is_finite() || *x < 0.)
            || r.refinement_scores.iter().rev().take(2).any(|x| *x > 1.)
            || r.crosscheck_scores.iter().any(|x| *x > 1.)
            || !sanity.passed
            || !sanity.actual.is_finite()
            || sanity.expected != -1.747558775 / (A / 2.)
            || sanity.tolerance != 2e-5 / A
            || (sanity.actual - sanity.expected).abs() > sanity.tolerance
        {
            return Err(invalid(
                "Madelung receipt lacks successful absolute physics proof",
            ));
        }
        let e = &levels
            .last()
            .ok_or_else(|| invalid("Madelung final level"))?
            .evaluation;
        if e.forces.len() != n
            || !e.energy.is_finite()
            || sanity.actual != e.energy / (n as f64 / 2.)
        {
            return Err(invalid("Madelung sanity differs from strict reference"));
        }
        Ok(())
    }
    /// Dynamic capacities retained by this receipt. Charge the outer container
    /// and path separately; retaining a comparison receipt transfers no work/calls.
    pub fn owned_heap_bytes(&self) -> Result<u64, CalibrationFailure> {
        let mut n = add(
            self.reference.owned_heap_bytes()?,
            (self.manifest_sha256.capacity() + self.specification_sha256.capacity()) as u64,
        )?;
        for (key, value) in &self.source_sha256 {
            // Conservative allowance for opaque BTreeMap nodes, as in the
            // special proof; actual string capacities are counted separately.
            n = add(n, add((key.capacity() + value.capacity()) as u64, 256)?)?;
        }
        if let Some(check) = &self.sanity {
            n = add(n, check.label.capacity() as u64)?;
        }
        n = add(n, failure_heap(self.failure.as_ref()))?;
        n = add(n, failure_heap(self.finalization_failure.as_ref()))?;
        Ok(n)
    }
}
#[derive(serde::Serialize)]
pub enum MadelungEvent<'a> {
    Reference(ReferenceEvent<'a>),
    Checked {
        check: &'a PhysicalCheck,
        usage: CandidateUsage,
    },
    Complete {
        outcome: &'a MadelungOutcome,
    },
}
pub fn verify_madelung_with_observer<F>(
    prepared: &Path,
    b: &mut CandidateBudget,
    mut observer: F,
) -> Result<MadelungOutcome, CalibrationFailure>
where
    F: FnMut(MadelungEvent<'_>) -> Result<(), CalibrationFailure>,
{
    let (p, ss, m) = validate_prepared(prepared)?;
    if !b.matches_limits(&p.limits) {
        return Err(invalid("Madelung continuing budget limits mismatch"));
    }
    if !p.family.starts_with("cubic-") {
        return Err(invalid("Madelung requires cubic family"));
    }
    let s = ss
        .iter()
        .find(|s| s.id.ends_with("-crystal"))
        .ok_or_else(|| invalid("crystal fixture"))?;
    let mut work = b.reference.work.clone();
    let mut calls = b.reference.calls;
    let last = p
        .reference_ladder
        .last()
        .ok_or_else(|| invalid("Madelung ladder"))?;
    for (kappa, [real, k]) in p
        .reference_ladder
        .iter()
        .map(|v| (p.reference_kappa_a, *v))
        .chain(p.crosscheck_kappa_a.iter().map(|ka| (*ka, *last)))
    {
        let w = planned_work(
            s,
            &ReferenceSettings {
                kappa: kappa / p.a,
                real_cutoff: real * p.a,
                kspace_cutoff: k / p.a,
            },
        )?;
        check_work([1; 3], &w, 0, &p.limits)?;
        work.real_terms = add(work.real_terms, w.real_terms)?;
        work.atom_k_terms = add(work.atom_k_terms, w.atom_k_terms)?;
        calls = add(calls, 1)?;
    }
    if calls > p.limits.max_accuracy_calls
        || add(work.real_terms, b.candidate_real_pair_upper_bound)?
            > p.limits.max_real_terms_per_family
        || work.atom_k_terms > p.limits.max_atom_k_terms_per_family
    {
        return Err(cap("complete Madelung work/call cap"));
    }
    let original = b.reference.owned_heap_baseline;
    let retained = add(
        protocol_bytes(&p),
        add(
            bytes(&ss),
            ss.iter().try_fold(0, |n, s| add(n, snapshot_bytes(s)))?,
        )?,
    )?;
    let reserve = add(
        retained,
        mul(
            (p.reference_ladder.len() + 2) as u64,
            add(mul(s.positions.len() as u64, 24)?, 4096)?,
        )?,
    )?;
    let baseline = add(add(original, b.external_owned_heap_bytes)?, reserve)?;
    if baseline > p.limits.max_owned_heap_bytes {
        return Err(cap("Madelung planned owned heap cap"));
    }
    b.reference.owned_heap_baseline = add(add(original, b.external_owned_heap_bytes)?, retained)?;
    let reference =
        converged_prepared_reference_with_observer(prepared, &s.id, &mut b.reference, |e| {
            observer(MadelungEvent::Reference(e))
        });
    b.reference.owned_heap_baseline = original;
    let reference = reference?;
    let mut out = MadelungOutcome {
        manifest_sha256: sha256(&prepared.join("manifest.json"))?,
        source_sha256: m.source_sha256,
        specification_sha256: m.specification_sha256,
        threads: p.threads,
        replicas: p.replicas,
        failure: reference
            .failure
            .clone()
            .or_else(|| reference.finalization_failure.clone()),
        finalization_failure: None,
        reference,
        sanity: None,
        usage: b.usage(),
    };
    if out.failure.is_none() {
        let result = (|| {
            validate_prepared(prepared)?;
            if out.reference.record.chosen_level != Some(p.reference_ladder.len() - 1)
                || !matches!(
                    out.reference.record.status,
                    Status::Passed | Status::RelativeUndefined
                )
            {
                return Err(failed("Madelung lacks strict absolute convergence"));
            }
            let e = &out
                .reference
                .record
                .levels
                .last()
                .ok_or_else(|| invalid("Madelung final level"))?
                .evaluation;
            let actual = e.energy / (s.positions.len() as f64 / 2.);
            let expected = -1.747558775 / (p.a / 2.);
            let tolerance = 2e-5 / p.a;
            let check = PhysicalCheck {
                label: "Madelung sign/normalization sanity only".into(),
                actual,
                expected,
                tolerance,
                passed: actual.is_finite() && (actual - expected).abs() <= tolerance,
            };
            let passed = check.passed;
            out.sanity = Some(check);
            observer(MadelungEvent::Checked {
                check: out
                    .sanity
                    .as_ref()
                    .ok_or_else(|| invalid("Madelung sanity"))?,
                usage: b.usage(),
            })?;
            b.reference.checkpoint()?;
            if !passed {
                return Err(failed("Madelung sign/normalization"));
            }
            Ok(())
        })();
        out.failure = result.err();
    }
    out.usage = b.usage();
    if let Err(e) = observer(MadelungEvent::Complete { outcome: &out }) {
        retain_finalization_failure(&mut out.failure, &mut out.finalization_failure, e);
    }
    Ok(out)
}
pub fn compare_madelung_replication(
    left: &MadelungOutcome,
    left_prepared: &Path,
    right: &MadelungOutcome,
    right_prepared: &Path,
) -> Result<PhysicalCheck, CalibrationFailure> {
    left.validate_success()?;
    right.validate_success()?;
    validate_prepared(left_prepared)?;
    validate_prepared(right_prepared)?;
    if left.manifest_sha256 != sha256(&left_prepared.join("manifest.json"))?
        || right.manifest_sha256 != sha256(&right_prepared.join("manifest.json"))?
        || left.source_sha256 != right.source_sha256
        || left.specification_sha256 != right.specification_sha256
        || left.threads != right.threads
        || left.replicas != 2
        || right.replicas != 3
        || left.failure.is_some()
        || right.failure.is_some()
    {
        return Err(invalid("Madelung replication proof identity/status"));
    }
    let a = left
        .sanity
        .as_ref()
        .ok_or_else(|| invalid("left Madelung sanity"))?
        .actual;
    let b = right
        .sanity
        .as_ref()
        .ok_or_else(|| invalid("right Madelung sanity"))?
        .actual;
    let tolerance = 1e-9 / A;
    Ok(PhysicalCheck {
        label: "Madelung finite-cell replication sanity".into(),
        actual: a,
        expected: b,
        tolerance,
        passed: (a - b).abs() <= tolerance,
    })
}

#[cfg(test)]
mod ownership_tests {
    use super::*;
    #[test]
    fn physical_madelung_failed_receipt_cannot_authorize_replication() {
        let temp = tempfile::tempdir().unwrap();
        let mut receipts = Vec::new();
        for replicas in [2, 3] {
            let dir = temp.path().join(format!("cubic-{replicas}"));
            let p = frozen_protocol("cubic", replicas, 1, true).unwrap();
            prepare(&dir, p.clone()).unwrap();
            let mut b = CandidateBudget::new(&p.limits, ReferenceBudget::new(&p.limits));
            let mut out = verify_madelung_with_observer(&dir, &mut b, |event| {
                if matches!(
                    event,
                    MadelungEvent::Reference(ReferenceEvent::Admitted { .. })
                ) {
                    Err(CalibrationFailure::Io("no oracle work authorized".into()))
                } else {
                    Ok(())
                }
            })
            .unwrap();
            assert_eq!(b.reference.calls, 0);
            assert!(out.reference.record.levels.is_empty());
            assert!(out.validate_success().is_err());
            assert!(!out.passed());
            // Internal corruption simulates the former public-field exploit and
            // an invalid future import. External callers must have no setters.
            out.failure = None;
            out.reference.failure = None;
            out.reference.failed_attempt = None;
            out.reference.record.status = Status::RelativeUndefined;
            out.reference.record.chosen_level = Some(5);
            let expected = -1.747558775 / (A / 2.);
            out.sanity = Some(PhysicalCheck {
                label: "invented without any reference evaluations".into(),
                actual: expected,
                expected,
                tolerance: 2e-5 / A,
                passed: true,
            });
            receipts.push((dir, out));
        }
        assert!(compare_madelung_replication(
            &receipts[0].1,
            &receipts[0].0,
            &receipts[1].1,
            &receipts[1].0
        )
        .is_err());
    }
    #[test]
    fn physical_special_failure_capacities_include_both_causes() {
        let temp = tempfile::tempdir().unwrap();
        let dir = temp.path().join("special");
        let p = frozen_protocol("special", 2, 1, true).unwrap();
        prepare(&dir, p.clone()).unwrap();
        let mut b = CandidateBudget::new(&p.limits, ReferenceBudget::new(&p.limits));
        let mut out = verify_special_with_observer(&dir, &mut b, |event| match event {
            PhysicalEvent::Admitted { .. } => Err(CalibrationFailure::Io("primary".into())),
            PhysicalEvent::Complete { .. } => Err(CalibrationFailure::Io("finalization".into())),
            _ => Ok(()),
        })
        .unwrap();
        assert_eq!(b.reference.calls, 0);
        let before = out.owned_heap_bytes().unwrap();
        let old =
            failure_heap(out.failure.as_ref()) + failure_heap(out.finalization_failure.as_ref());
        if let Some(CalibrationFailure::Io(s)) = &mut out.failure {
            s.reserve_exact(4096);
        } else {
            panic!("missing primary");
        }
        if let Some(CalibrationFailure::Io(s)) = &mut out.finalization_failure {
            s.reserve_exact(8192);
        } else {
            panic!("missing finalization");
        }
        let new =
            failure_heap(out.failure.as_ref()) + failure_heap(out.finalization_failure.as_ref());
        assert_eq!(out.owned_heap_bytes().unwrap() - before, new - old);
        assert!(!out.passed());
    }
    #[test]
    fn physical_madelung_owned_capacity_tracks_retained_growth() {
        let limits = frozen_protocol("cubic", 2, 1, true).unwrap().limits;
        let budget = CandidateBudget::new(&limits, ReferenceBudget::new(&limits));
        let mut out = MadelungOutcome {
            manifest_sha256: String::new(),
            source_sha256: Default::default(),
            specification_sha256: String::new(),
            threads: 1,
            replicas: 2,
            reference: ReferenceOutcome {
                record: ReferenceRecord {
                    snapshot_id: String::new(),
                    snapshot_sha256: String::new(),
                    levels: vec![],
                    crosschecks: vec![],
                    chosen_level: None,
                    status: Status::Failed,
                },
                audits: vec![],
                refinement_scores: vec![],
                crosscheck_scores: vec![],
                failure: None,
                finalization_failure: None,
                failed_attempt: None,
                usage: ReferenceUsage {
                    work: WorkCounts {
                        real_terms: 0,
                        atom_k_terms: 0,
                    },
                    calls: 0,
                    elapsed_ns: 0,
                },
            },
            sanity: None,
            failure: None,
            finalization_failure: None,
            usage: budget.usage(),
        };
        let before = out.owned_heap_bytes().unwrap();
        out.manifest_sha256.reserve_exact(123);
        out.reference.record.snapshot_id.reserve_exact(157);
        assert_eq!(
            out.owned_heap_bytes().unwrap() - before,
            (out.manifest_sha256.capacity() + out.reference.record.snapshot_id.capacity()) as u64
        );
        let before = out.owned_heap_bytes().unwrap();
        let key = "core source".to_string();
        let value = "digest".to_string();
        let expected = (key.capacity() + value.capacity()) as u64 + 256;
        out.source_sha256.insert(key, value);
        assert_eq!(out.owned_heap_bytes().unwrap() - before, expected);
        let before = out.owned_heap_bytes().unwrap();
        let mut primary = String::from("primary");
        primary.reserve_exact(4096);
        let mut finalization = String::from("finalization");
        finalization.reserve_exact(8192);
        let extra = (primary.capacity() + finalization.capacity()) as u64;
        out.failure = Some(CalibrationFailure::Io(primary));
        out.finalization_failure = Some(CalibrationFailure::Io(finalization));
        assert_eq!(out.owned_heap_bytes().unwrap() - before, extra);
        assert_eq!(budget.reference.calls, 0);
    }
}

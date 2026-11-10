//! Pure measured-cost decisions. Held-out outcomes never enter ranking.
use super::*;

pub fn ordinary_timing_summary(
    cold_ns: u64,
    rebuild_ns: u64,
    batches: Vec<u64>,
    p: &Protocol,
) -> Result<TimingRecord, CalibrationFailure> {
    if p.warmups != 5
        || p.batches != 5
        || p.batch_calls != 20
        || batches.len() != 5
        || cold_ns == 0
        || rebuild_ns == 0
        || batches.contains(&0)
    {
        return Err(invalid("complete positive ordinary 5x20 timing required"));
    }
    let mut sorted = batches.clone();
    sorted.sort_unstable();
    Ok(TimingRecord {
        reset_included: true,
        cold_ns,
        rebuild_ns,
        batch_total_ns: batches,
        calls_per_batch: 20,
        median_ns_per_call: sorted[2] as f64 / 20.,
        q1_ns_per_call: sorted[1] as f64 / 20.,
        q3_ns_per_call: sorted[3] as f64 / 20.,
    })
}
fn valid_timing(t: &TimingRecord, p: &Protocol) -> bool {
    t.reset_included
        && ordinary_timing_summary(t.cold_ns, t.rebuild_ns, t.batch_total_ns.clone(), p)
            .is_ok_and(|expected| expected == *t)
}
/// All feasible rows must have accepted ordinary cost. Missing/invalid timing
/// prevents a fastest-feasible claim even when another candidate was measured.
pub fn select(
    candidates: &[CandidateReport],
    ids: &[&str],
    target: f64,
    p: &Protocol,
) -> Result<TargetDecision, CalibrationFailure> {
    if !p.targets.contains(&target)
        || !target.is_finite()
        || p.tolerances.cost_tie_fraction != 0.05
        || candidates
            .iter()
            .map(|c| c.id)
            .collect::<std::collections::BTreeSet<_>>()
            .len()
            != candidates.len()
    {
        return Err(invalid("target/tie protocol or duplicate candidate IDs"));
    }
    let mut feasible = Vec::new();
    for c in candidates {
        if feasible_targets(c, ids, p)?.contains(&target) {
            if [
                c.parameters.kappa,
                c.parameters.real_cutoff,
                c.parameters.kspace_cutoff,
            ]
            .iter()
            .any(|x| !x.is_finite() || *x <= 0.)
            {
                return Err(invalid("finite positive cost parameters required"));
            }
            check_work(
                c.actual_mesh,
                &WorkCounts {
                    real_terms: 0,
                    atom_k_terms: 0,
                },
                0,
                &p.limits,
            )?;
            let timing = c
                .timing
                .as_ref()
                .filter(|t| valid_timing(t, p))
                .ok_or_else(|| {
                    invalid(
                        "every feasible candidate needs complete finite accepted ordinary timing",
                    )
                })?;
            feasible.push((c, timing.median_ns_per_call));
        }
    }
    let best = feasible
        .iter()
        .map(|(_, ns)| *ns)
        .fold(f64::INFINITY, f64::min);
    feasible.retain(|(_, ns)| *ns <= best * 1.05);
    feasible.sort_by(|(a, _), (b, _)| {
        let product = |c: &CandidateReport| c.actual_mesh.iter().product::<usize>();
        product(a)
            .cmp(&product(b))
            .then_with(|| {
                a.parameters
                    .real_cutoff
                    .total_cmp(&b.parameters.real_cutoff)
            })
            .then_with(|| a.parameters.kappa.total_cmp(&b.parameters.kappa))
            .then_with(|| {
                a.parameters
                    .kspace_cutoff
                    .total_cmp(&b.parameters.kspace_cutoff)
            })
            .then_with(|| a.id.cmp(&b.id))
    });
    let selected_id = feasible.first().map(|(c, _)| c.id);
    Ok(TargetDecision {
        target,
        selected_id,
        selection_status: if selected_id.is_some() {
            Status::Passed
        } else {
            Status::TargetNotReached
        },
        confirmation_rows: vec![],
        confirmation_status: Status::Pending,
    })
}
/// Complete one frozen decision, consuming exactly that ID's complete holdout
/// roster. No candidate list or cost data exists here to support a fallback.
pub fn finish_confirmation(
    d: &mut TargetDecision,
    c: &CandidateReport,
    ids: &[&str],
    p: &Protocol,
) -> Result<(), CalibrationFailure> {
    if d.selected_id != Some(c.id)
        || d.selection_status != Status::Passed
        || d.confirmation_status != Status::Pending
        || !d.confirmation_rows.is_empty()
    {
        return Err(invalid(
            "confirmation requires untouched frozen selected ID",
        ));
    }
    let passes = feasible_targets(c, ids, p)?.contains(&d.target);
    d.confirmation_rows = c.rows.clone();
    d.confirmation_status = if passes {
        Status::Passed
    } else {
        Status::ConfirmationFailed
    };
    Ok(())
}
/// Explicit reduced-unit fragment for a confirmed target only. The caller must
/// append its immutable report/manifest identities when writing the artifact.
pub fn selected_fragment(
    d: &TargetDecision,
    c: &CandidateReport,
    p: &Protocol,
) -> Result<Option<String>, CalibrationFailure> {
    if d.selected_id.is_none()
        || d.selection_status != Status::Passed
        || d.confirmation_status != Status::Passed
    {
        return Ok(None);
    }
    if d.selected_id != Some(c.id)
        || c.family != p.family
        || d.confirmation_rows.is_empty()
        || !p.targets.contains(&d.target)
    {
        return Err(invalid("fragment candidate/confirmed decision mismatch"));
    }
    for row in &d.confirmation_rows {
        if row.status != Status::Passed || !passes_target(&row.metrics, d.target, p, true)? {
            return Err(invalid("fragment confirmation gate mismatch"));
        }
    }
    let q = &c.parameters;
    if [q.kappa, q.real_cutoff, q.kspace_cutoff]
        .iter()
        .any(|x| !x.is_finite() || *x <= 0.)
    {
        return Err(invalid("fragment parameters"));
    }
    Ok(Some(format!("# empirical helper calibration; phase profiling pending\n# Merge into an existing units = \"reduced\" deck.\n[electrostatics]\nmethod = \"spme\"\nkappa = {:?}\nreal_cutoff = {:?}\nkspace_cutoff = {:?}\n", q.kappa,q.real_cutoff,q.kspace_cutoff)))
}

//! Deterministic retained evidence presentation; no compute or reselection.
use super::*;
use journal::{Artifact, Receipt};
use std::io::{BufRead, BufReader};
fn csv(s: &str) -> String {
    format!("\"{}\"", s.replace('"', "\"\""))
}
fn xml(s: &str) -> String {
    s.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
        .replace('"', "&quot;")
        .replace('\'', "&apos;")
}
pub fn tables(candidates: &[CandidateReport], decisions: &[TargetDecision]) -> (String, String) {
    let mut rows = String::from("candidate_id,fixture,shift,status,relative_rms,absolute_rms,max_atom,max_component,E_per_atom_signed,E_per_atom_abs,W_per_atom_signed,W_per_atom_abs,net_x,net_y,net_z,net_norm,cold_ns,rebuild_ns,warm_median_ns,warm_q1_ns,warm_q3_ns,raw_batch_ns\n");
    let mut sorted: Vec<_> = candidates.iter().collect();
    sorted.sort_by_key(|c| c.id);
    for c in sorted {
        let costs = c
            .timing
            .as_ref()
            .map(|t| {
                format!(
                    "{},{},{},{},{},{}",
                    t.cold_ns,
                    t.rebuild_ns,
                    t.median_ns_per_call,
                    t.q1_ns_per_call,
                    t.q3_ns_per_call,
                    csv(&t
                        .batch_total_ns
                        .iter()
                        .map(u64::to_string)
                        .collect::<Vec<_>>()
                        .join(";"))
                )
            })
            .unwrap_or_else(|| ",,,,,".into());
        if c.rows.is_empty() {
            rows.push_str(&format!(
                "{},{},,{:?},,,,,,,,,,,,,,{costs}\n",
                c.id,
                csv("no completed accuracy rows"),
                c.status
            ));
        }
        for row in &c.rows {
            let m = &row.metrics;
            rows.push_str(&format!(
                "{},{},{},{:?},{},{},{},{},{},{},{},{},{},{},{},{},{costs}\n",
                c.id,
                csv(&row.snapshot_id),
                row.shift_index,
                row.status,
                m.force_rms_rel.map(|x| x.to_string()).unwrap_or_default(),
                m.force_rms_abs,
                m.max_atom_force_error,
                m.max_component_force_error,
                m.energy_per_atom_signed,
                m.energy_per_atom_abs,
                m.w_per_atom_signed,
                m.w_per_atom_abs,
                m.net_force[0],
                m.net_force[1],
                m.net_force[2],
                m.net_force_norm
            ));
        }
    }
    let mut targets = String::from("target,frozen_id,selection_status,confirmation_status\n");
    for d in decisions {
        targets.push_str(&format!(
            "{},{},{:?},{:?}\n",
            d.target,
            d.selected_id.map(|x| x.to_string()).unwrap_or_default(),
            d.selection_status,
            d.confirmation_status
        ));
    }
    (rows, targets)
}
pub fn svg(candidates: &[CandidateReport], decisions: &[TargetDecision], label: &str) -> String {
    let mut points: Vec<_> = candidates
        .iter()
        .filter_map(|c| {
            let t = c.timing.as_ref()?;
            let error = c
                .rows
                .iter()
                .filter_map(|r| r.metrics.force_rms_rel)
                .fold(0., f64::max);
            (t.median_ns_per_call.is_finite()
                && t.median_ns_per_call > 0.
                && error.is_finite()
                && error > 0.)
                .then_some((c, t.median_ns_per_call.log10(), error.log10()))
        })
        .collect();
    points.sort_by_key(|(c, _, _)| c.id);
    let xmin = points.iter().map(|(_, x, _)| *x).fold(2., f64::min).floor();
    let xmax = points
        .iter()
        .map(|(_, x, _)| *x)
        .fold(xmin + 1., f64::max)
        .ceil()
        .max(xmin + 1.);
    let ymin = points
        .iter()
        .map(|(_, _, y)| *y)
        .fold(-5., f64::min)
        .floor();
    let ymax = points
        .iter()
        .map(|(_, _, y)| *y)
        .fold(-2., f64::max)
        .ceil()
        .max(ymin + 1.);
    let x = |v| 90. + (v - xmin) / (xmax - xmin) * 700.;
    let y = |v| 400. - (v - ymin) / (ymax - ymin) * 310.;
    let height = 510 + decisions.len() * 24;
    let mut out = format!("<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"880\" height=\"{height}\" viewBox=\"0 0 880 {height}\"><rect width=\"100%\" height=\"100%\" fill=\"white\"/><g font-family=\"sans-serif\" font-size=\"13\" fill=\"#17212b\"><text x=\"30\" y=\"28\">{}</text><text x=\"30\" y=\"50\">Observed exact-snapshot error; ordinary uninstrumented cost. No speedup claim.</text><path d=\"M90 80 V400 H800\" stroke=\"#17212b\" fill=\"none\"/><text x=\"310\" y=\"450\">Ordinary ns/call (log10 scale)</text><text x=\"15\" y=\"76\">Relative RMS (log10)</text>", xml(label));
    for i in (xmin as i32)..=(xmax as i32) {
        out.push_str(&format!(
            "<text x=\"{:.2}\" y=\"423\">1e{i}</text>",
            x(f64::from(i))
        ));
    }
    for i in (ymin as i32)..=(ymax as i32) {
        out.push_str(&format!(
            "<text x=\"34\" y=\"{:.2}\">1e{i}</text>",
            y(f64::from(i))
        ));
    }
    for t in [1e-3_f64, 1e-4_f64] {
        let py = y(t.log10());
        out.push_str(&format!("<path d=\"M90 {py:.2} H800\" stroke=\"#9ba5ad\" stroke-dasharray=\"5 4\"/><text x=\"700\" y=\"{:.2}\">target {t}</text>", py-5.));
    }
    for (c, px, py) in points {
        let color = if c.id == 126 {
            "#b06716"
        } else if c.status == Status::Passed {
            "#137c66"
        } else {
            "#b13a43"
        };
        out.push_str(&format!(
            "<circle cx=\"{:.2}\" cy=\"{:.2}\" r=\"4\" fill=\"{color}\"><title>{}</title></circle>",
            x(px),
            y(py),
            xml(&format!("ID{} {:?}; {}", c.id, c.status, c.family))
        ));
    }
    out.push_str("<text x=\"30\" y=\"477\">Green: passing selection; red: failed; ochre: baseline126. Untimed rows remain in CSV.</text>");
    for (i, d) in decisions.iter().enumerate() {
        out.push_str(&format!(
            "<text x=\"30\" y=\"{}\">{}</text>",
            505 + 24 * i,
            xml(&format!(
                "Target {}: frozen ID {:?}; selection {:?}; held-out {:?}",
                d.target, d.selected_id, d.selection_status, d.confirmation_status
            ))
        ));
    }
    out.push_str("</g></svg>\n");
    out
}
#[derive(serde::Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Memory {
    pub schema: String,
    pub manifest_sha256: String,
    pub children: Vec<ChildMemory>,
}
#[derive(serde::Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ChildMemory {
    pub phase: String,
    pub actual_child_exit: Option<i32>,
    pub monitor_exit: Option<i32>,
    pub wrapper_exit: i32,
    pub peak_rss_bytes: Option<u64>,
    pub active_ns: u64,
    pub timed_out: bool,
    pub command: Vec<String>,
    pub stdout: String,
    pub stderr: String,
}
pub fn accepted_memory(m: &Memory, epoch: &str) -> bool {
    m.schema == "modena-spme-memory-v1"
        && m.manifest_sha256 == epoch
        && ["verify", "sweep", "confirm", "profile"]
            .iter()
            .all(|phase| m.children.iter().any(|c| c.phase == *phase))
        && m.children.iter().all(|c| {
            !c.timed_out
                && c.monitor_exit == Some(0)
                && c.wrapper_exit == 0
                && matches!(c.actual_child_exit, Some(0 | 4))
                && c.peak_rss_bytes.is_some_and(|n| n > 0 && n <= 536_870_912)
        })
        && m.children
            .iter()
            .filter(|c| ["verify", "sweep", "confirm", "profile"].contains(&c.phase.as_str()))
            .try_fold(0_u64, |sum, c| sum.checked_add(c.active_ns))
            .is_some_and(|sum| sum < 1_800_000_000_000)
}
fn chain(r: &Receipt, rows: &mut Vec<Receipt>) -> Result<(), CalibrationFailure> {
    if let Some(p) = &r.prior {
        chain(&Receipt::load(&p.path, &p.sha256)?, rows)?;
    }
    rows.push(r.clone());
    Ok(())
}
fn decode<T: serde::de::DeserializeOwned>(
    value: &serde_json::Value,
) -> Result<T, CalibrationFailure> {
    serde_json::from_value(value.clone()).map_err(|e| invalid(e.to_string()))
}
pub fn run(args: &cli::Arguments) -> Result<(), CalibrationFailure> {
    let receipt = Receipt::load(args.path("--campaign")?, args.value("--expected-sha256")?)?;
    let (p, snapshots, manifest) = validate_prepared(&receipt.prepared)?;
    let memory_path = args.path("--memory")?;
    journal::check_hash(memory_path, args.value("--memory-sha256")?)?;
    let memory: Memory = journal::read_json(memory_path)?;
    let mut phases = Vec::new();
    chain(&receipt, &mut phases)?;
    let grid = candidate_grid(&p, snapshots[0].cell)?;
    let mut candidates = BTreeMap::<usize, CandidateReport>::new();
    let mut decisions: Vec<TargetDecision> = p
        .targets
        .iter()
        .map(|t| TargetDecision {
            target: *t,
            selected_id: None,
            selection_status: Status::Pending,
            confirmation_rows: vec![],
            confirmation_status: Status::Pending,
        })
        .collect();
    let mut ordinary = String::from("phase,candidate_id,stage,cold_ns,rebuild_ns,warmups_completed,raw_batches_ns,timing_calls,failure,finalization_failure\n");
    let mut references = Vec::new();
    let mut variation = Vec::<serde_json::Value>::new();
    let mut baseline_confirmation = false;
    let mut complete_grid = false;
    let mut baseline_heldout = Vec::new();
    for phase in &phases {
        if let Some(a) = phase.artifacts.get("decisions") {
            decisions = journal::read_json(&a.path)?;
        }
        if let Some(a) = phase.artifacts.get("candidates") {
            let rows: Vec<CandidateReport> = journal::read_json(&a.path)?;
            complete_grid = rows.len() == grid.len();
            for c in rows {
                candidates.insert(c.id, c);
            }
        }
        if let Some(a) = phase.artifacts.get("baseline-selection") {
            let row: Option<CandidateOutcome> = journal::read_json(&a.path)?;
            if let Some(o) = row {
                candidates.insert(126, o.report);
            }
        }
        if let Some(a) = phase.artifacts.get("baseline-confirmation") {
            let o: CandidateOutcome = journal::read_json(&a.path)?;
            baseline_confirmation = o.failure.is_none()
                && o.finalization_failure.is_none()
                && o.report.rows.len()
                    == snapshots
                        .iter()
                        .filter(|s| s.role == FixtureRole::Confirmation)
                        .count()
                        * 7;
            baseline_heldout.push(o.report);
        }
        if let Some(a) = phase.artifacts.get("events") {
            let mut current = None;
            journal::visit_events(&a.path, &a.sha256, &receipt.manifest_sha256, |event| {
                if event.kind == "ordinary" || event.kind == "baseline-ordinary" {
                    let (id, value) = if event.kind == "ordinary" {
                        (event.data["candidate_id"].to_string(), &event.data["event"])
                    } else {
                        ("126".into(), &event.data)
                    };
                    let o = &value["outcome"];
                    let field = |key: &str| {
                        if o[key].is_null() {
                            String::new()
                        } else {
                            o[key].to_string()
                        }
                    };
                    ordinary.push_str(&format!(
                        "{},{},{},{},{},{},{},{},{},{}\n",
                        phase.phase,
                        id,
                        csv(&value["stage"].to_string()),
                        field("cold_ns"),
                        field("rebuild_ns"),
                        field("warmups_completed"),
                        csv(&field("batch_total_ns")),
                        field("timing_calls"),
                        csv(&field("failure")),
                        csv(&field("finalization_failure"))
                    ));
                }
                if event.kind == "reference" {
                    if let Some(v) = event.data.get("Complete").and_then(|v| v.get("outcome")) {
                        let o: ReferenceOutcome = decode(v)?;
                        references.push(o.record);
                    }
                }
                if ["candidate", "baseline-selection", "baseline-confirmation"]
                    .contains(&event.kind.as_str())
                {
                    if let Some(a) = event.data.get("Admitted") {
                        current = a
                            .get("candidate_id")
                            .and_then(serde_json::Value::as_u64)
                            .map(|id| id as usize);
                    }
                    if let Some(v) = event.data.get("Complete").and_then(|v| v.get("outcome")) {
                        let o: CandidateOutcome = decode(v)?;
                        variation.push(serde_json::json!({"phase":phase.phase,"id":o.report.id,"variation":o.translation_variation,"failure":o.failure,"finalization_failure":o.finalization_failure}));
                    }
                    if phase.phase != "confirm" {
                        if let (Some(id), Some(v)) = (
                            current,
                            event.data.get("CandidateRow").and_then(|v| v.get("row")),
                        ) {
                            let row: MetricRow = decode(v)?;
                            if let std::collections::btree_map::Entry::Vacant(entry) =
                                candidates.entry(id)
                            {
                                entry.insert(if id == 126 {
                                    candidates::baseline_report(&p, snapshots[0].cell)?
                                } else {
                                    grid.iter()
                                        .find(|c| c.id == id)
                                        .cloned()
                                        .ok_or_else(|| invalid("journal ID"))?
                                });
                            }
                            let c = candidates
                                .get_mut(&id)
                                .ok_or_else(|| invalid("journal candidate"))?;
                            if !c.rows.iter().any(|old| {
                                old.snapshot_id == row.snapshot_id
                                    && old.shift_index == row.shift_index
                            }) {
                                c.rows.push(row);
                            }
                        }
                    }
                }
                Ok(())
            })?;
        }
    }
    let mut phase_records = Vec::new();
    let mut profile_complete = false;
    if let Some(path) = args.values.get("--profiles") {
        let path = Path::new(path);
        journal::check_hash(path, args.value("--profiles-sha256")?)?;
        let pr = Receipt::load(
            &path
                .parent()
                .ok_or_else(|| invalid("profile parent"))?
                .join("summary.json"),
            args.value("--profile-receipt-sha256")?,
        )?;
        let parent = Artifact::capture(args.path("--campaign")?)?;
        if pr.phase != "profile"
            || pr.manifest_sha256 != receipt.manifest_sha256
            || pr.prior.as_ref().map(|a| &a.sha256) != Some(&parent.sha256)
            || pr.artifact("profiles")?.sha256 != args.value("--profiles-sha256")?
        {
            return Err(invalid("profile confirmation chain mismatch"));
        }
        let mut schedules = BTreeMap::<usize, Vec<profile::ProfileSample>>::new();
        for line in BufReader::new(std::fs::File::open(path).map_err(io)?).lines() {
            let s: profile::ProfileSample =
                serde_json::from_str(&line.map_err(io)?).map_err(|e| invalid(e.to_string()))?;
            phase_records.push(s.record.clone());
            schedules.entry(s.record.candidate_id).or_default().push(s);
        }
        let mut ids = vec![126];
        for d in &decisions {
            if d.confirmation_status == Status::Passed {
                if let Some(id) = d.selected_id {
                    if !ids.contains(&id) {
                        ids.push(id);
                    }
                }
            }
        }
        let representative = snapshots
            .iter()
            .find(|s| s.role == FixtureRole::Selection)
            .ok_or_else(|| invalid("profile representative"))?;
        profile_complete = pr.require_success().is_ok()
            && schedules.len() == ids.len()
            && ids.iter().all(|id| {
                schedules.get(id).is_some_and(|ss| {
                    ss.len() == 106
                        && ss.iter().enumerate().all(|(i, s)| {
                            s.record.threads == p.threads
                                && s.record.calls == 1
                                && s.record.phases_ns.len() == 1
                                && s.record.snapshot_id == representative.id
                                && if i == 0 {
                                    s.stage == profile::ProfileStage::Cold
                                        && s.call == 0
                                        && s.batch.is_none()
                                } else if i <= 5 {
                                    s.stage == profile::ProfileStage::Warmup
                                        && s.call == i - 1
                                        && s.batch.is_none()
                                } else {
                                    s.stage == profile::ProfileStage::Batch
                                        && s.call == (i - 6) % 20
                                        && s.batch == Some((i - 6) / 20)
                                }
                        })
                })
            });
    }
    let all: Vec<_> = candidates.into_values().collect();
    let selection_count = snapshots
        .iter()
        .filter(|s| s.role == FixtureRole::Selection)
        .count()
        * 7;
    let full = complete_grid
        && baseline_confirmation
        && all.iter().all(|c| c.rows.len() == selection_count)
        && all
            .iter()
            .find(|c| c.id == 126)
            .is_some_and(|c| c.timing.is_some())
        && decisions.iter().all(|d| {
            matches!(
                d.selection_status,
                Status::Passed | Status::TargetNotReached
            ) && (d.selected_id.is_none()
                || matches!(
                    d.confirmation_status,
                    Status::Passed | Status::ConfirmationFailed
                ))
        })
        && phases.iter().all(|r| r.require_success().is_ok())
        && profile_complete
        && accepted_memory(&memory, &receipt.manifest_sha256);
    let status = if full {
        Status::FullScopeMeasured
    } else if phases.iter().any(|r| r.require_success().is_err()) {
        Status::Failed
    } else if !complete_grid || receipt.phase != "confirm" {
        Status::NotCalibrated
    } else {
        Status::OfflineHelperCompletePhaseProfilingPending
    };
    let capture = |prog: &str, argv: &[&str]| -> String {
        Command::new(prog)
            .args(argv)
            .output()
            .ok()
            .filter(|o| o.status.success())
            .map(|o| String::from_utf8_lossy(&o.stdout).trim().into())
            .unwrap_or_else(|| "unavailable".into())
    };
    let report = CampaignReport {
        protocol: p.clone(),
        provenance: Provenance {
            git_head: manifest.git_head,
            source_sha256: manifest.source_sha256,
            protocol_sha256: manifest.protocol_sha256,
            snapshot_sha256: manifest.snapshots_sha256,
            reference_algorithm: "all-image-full-pm-k-v1".into(),
            rustc_verbose: capture("rustc", &["-Vv"]),
            cargo_version: capture("cargo", &["-V"]),
            os: std::env::consts::OS.into(),
            cpu: capture("sysctl", &["-n", "machdep.cpu.brand_string"]),
            threads: p.threads,
            profile: "release".into(),
            features: vec![],
            command: std::env::args().collect(),
            utc_started: "see launcher registration.json".into(),
            family_elapsed_ns: receipt.usage.candidate.reference.elapsed_ns,
            measured_peak_rss_bytes: memory
                .children
                .iter()
                .filter_map(|c| c.peak_rss_bytes)
                .max(),
        },
        snapshot_ids: snapshots.iter().map(|s| s.id.clone()).collect(),
        references,
        candidates: all,
        decisions,
        phase_records,
        status,
    };
    let out = args.path("--out")?;
    std::fs::create_dir(out).map_err(io)?;
    super::write_new(&out.join("summary.json"), &json_bytes(&report)?)?;
    super::write_new(
        &out.join("diagnostics.json"),
        &json_bytes(
            &serde_json::json!({"phases":phases,"translation_variation":variation,"full_scope_measured":full,"memory":memory_path,"scientific_status_separate":true}),
        )?,
    )?;
    super::write_new(&out.join("ordinary-costs.csv"), ordinary.as_bytes())?;
    let (rows, targets) = tables(&report.candidates, &report.decisions);
    super::write_new(&out.join("errors.csv"), rows.as_bytes())?;
    super::write_new(&out.join("targets.csv"), targets.as_bytes())?;
    super::write_new(
        &out.join("cost-error.svg"),
        svg(&report.candidates, &report.decisions, &p.family).as_bytes(),
    )?;
    let mut heldout = baseline_heldout;
    for d in &report.decisions {
        if let Some(id) = d.selected_id {
            if let Some(c) = report.candidates.iter().find(|c| c.id == id) {
                let mut row = c.clone();
                row.family = format!("{} held-out target {}", p.family, d.target);
                row.rows = d.confirmation_rows.clone();
                row.timing = None;
                row.status = d.confirmation_status.clone();
                heldout.push(row);
            }
        }
    }
    super::write_new(
        &out.join("held-out-errors.csv"),
        tables(&heldout, &[]).0.as_bytes(),
    )?;
    Ok(())
}

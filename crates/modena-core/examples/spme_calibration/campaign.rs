//! Six-phase orchestration. Scientific work stays in the shared bounded engines.
use super::*;
use journal::{Artifact, Journal, Receipt, Usage};
use std::cell::RefCell;

pub fn validate_fragment_schema(fragment: &str, q: &Parameters) -> Result<(), CalibrationFailure> {
    let deck = format!("units = \"reduced\"\n[structure]\ntype = \"fcc\"\nlattice_constant = 100.0\nreplicas = [1,1,1]\nelement = \"Ar\"\n[model]\ntype = \"lennard_jones\"\nsigma = 1.0\nepsilon = 1.0\ncutoff = 2.5\n[integrator]\ntype = \"velocity_verlet\"\ndt = 0.001\n[run]\nsteps = 1\n{fragment}");
    let parsed = modena_core::InputConfig::from_str(&deck).map_err(|e| invalid(e.to_string()))?;
    let e = parsed
        .electrostatics
        .ok_or_else(|| invalid("fragment electrostatics missing"))?;
    if !matches!(e.method, modena_core::EwaldMethod::Spme)
        || e.kappa.to_bits() != q.kappa.to_bits()
        || e.real_cutoff.to_bits() != q.real_cutoff.to_bits()
        || e.kspace_cutoff.to_bits() != q.kspace_cutoff.to_bits()
    {
        return Err(invalid(
            "selected fragment failed exact explicit schema roundtrip",
        ));
    }
    Ok(())
}

fn parent_elapsed(args: &cli::Arguments) -> Result<u64, CalibrationFailure> {
    args.values
        .get("--elapsed-ns")
        .map(|s| {
            s.parse()
                .map_err(|_| invalid("invalid cumulative elapsed nanoseconds"))
        })
        .unwrap_or(Ok(0))
}
fn runtime(p: &Protocol) -> Result<(), CalibrationFailure> {
    if cfg!(debug_assertions)
        || std::env::var("RAYON_NUM_THREADS").ok().as_deref()
            != Some(p.threads.to_string().as_str())
        || rayon::current_num_threads() != p.threads
    {
        return Err(invalid(
            "compute requires release and matching RAYON_NUM_THREADS/runtime pool",
        ));
    }
    Ok(())
}
fn prepared_paths(args: &cli::Arguments) -> Result<PathBuf, CalibrationFailure> {
    let protocol = args.path("--protocol")?.canonicalize().map_err(io)?;
    let snapshots = args.path("--snapshots")?.canonicalize().map_err(io)?;
    let dir = protocol
        .parent()
        .ok_or_else(|| invalid("protocol parent"))?;
    if protocol != dir.join("protocol.json") || snapshots != dir.join("snapshots.json") {
        return Err(invalid("protocol/snapshots must be one prepared directory"));
    }
    validate_prepared(dir)?;
    Ok(dir.into())
}
fn fresh_receipt(
    phase: &str,
    prepared: &Path,
    b: &CandidateBudget,
    prior: Option<Artifact>,
) -> Result<Receipt, CalibrationFailure> {
    Ok(Receipt {
        schema: "modena-spme-phase-v1".into(),
        phase: phase.into(),
        prepared: prepared.canonicalize().map_err(io)?,
        manifest_sha256: sha256(&prepared.join("manifest.json"))?,
        prior,
        inputs: BTreeMap::new(),
        artifacts: BTreeMap::new(),
        usage: Usage::capture(b, 0, 0),
        failure: None,
        finalization_failure: None,
    })
}
fn persist<T: serde::Serialize>(
    r: &mut Receipt,
    out: &Path,
    key: &str,
    value: &T,
) -> Result<(), CalibrationFailure> {
    let path = out.join(format!("{key}.json"));
    super::write_new(&path, &json_bytes(value)?)?;
    r.artifacts.insert(key.into(), Artifact::capture(&path)?);
    Ok(())
}
fn record_failure(r: &mut Receipt, f: CalibrationFailure) {
    if r.failure.as_ref() == Some(&f) || r.finalization_failure.as_ref() == Some(&f) {
        return;
    }
    if r.failure.is_none() {
        r.failure = Some(f);
    } else if let Some(previous) = r.finalization_failure.take() {
        r.finalization_failure = Some(io(format!("{previous}; additional finalization: {f}")));
    } else {
        r.finalization_failure = Some(f);
    }
}
fn complete(
    r: &mut Receipt,
    out: &Path,
    summary: &Path,
    b: &CandidateBudget,
    timing: u64,
    profiles: u64,
    result: Result<(), CalibrationFailure>,
) -> Result<(), CalibrationFailure> {
    if let Err(f) = result {
        record_failure(r, f);
    }
    r.usage = Usage::capture(b, timing, profiles);
    let final_result = (|| {
        validate_prepared(&r.prepared)?;
        r.artifacts.insert(
            "events".into(),
            Artifact::capture(&out.join("events.jsonl"))?,
        );
        r.write_new(summary)?;
        Ok(())
    })();
    if let Err(f) = final_result {
        eprintln!(
            "primary phase failure: {:?}; finalization failure: {f}",
            r.failure
        );
        return Err(f);
    }
    r.require_success()
}
pub fn run(args: cli::Arguments) -> Result<(), CalibrationFailure> {
    match args.mode.as_str() {
        "prepare" => {
            let p = frozen_protocol(
                args.value("--family")?,
                args.value("--replicas")?
                    .parse()
                    .map_err(|_| invalid("invalid replicas"))?,
                args.value("--threads")?
                    .parse()
                    .map_err(|_| invalid("invalid threads"))?,
                args.extended,
            )?;
            prepare(args.path("--out")?, p)?;
            Ok(())
        }
        "verify" => verify(&args),
        "sweep" => sweep(&args),
        "confirm" => confirm(&args),
        "profile" => profile(&args),
        "report" => super::report::run(&args),
        _ => Err(invalid("unknown mode")),
    }
}
fn verify(args: &cli::Arguments) -> Result<(), CalibrationFailure> {
    let prepared = prepared_paths(args)?;
    let (p, _, _) = validate_prepared(&prepared)?;
    runtime(&p)?;
    let reference = ReferenceBudget::resume(
        &p.limits,
        WorkCounts {
            real_terms: 0,
            atom_k_terms: 0,
        },
        0,
        std::time::Duration::from_nanos(parent_elapsed(args)?),
        0,
    )?;
    let mut b = CandidateBudget::new(&p.limits, reference);
    let out = prepared.join("verification");
    std::fs::create_dir(&out).map_err(io)?;
    let mut r = fresh_receipt("verify", &prepared, &b, None)?;
    let mut j = Journal::create(&out.join("events.jsonl"), r.manifest_sha256.clone())?;
    b.external_owned_heap_bytes =
        j.owned_heap_bytes() + r.owned_heap_bytes() + candidates::protocol_bytes(&p);
    let result = (|| {
        if p.family == "special-3" {
            if args.values.contains_key("--special") || args.values.contains_key("--madelung") {
                return Err(invalid(
                    "special verification cannot import another physical proof",
                ));
            }
            let proof = physical::verify_special_with_observer(&prepared, &mut b, |e| {
                j.record("physical", &e)
            })?;
            persist(&mut r, &out, "physical", &proof)?;
            if let Some(f) = proof.failure() {
                r.failure = Some(f.clone());
            }
            if let Some(f) = proof.finalization_failure() {
                r.finalization_failure = Some(f.clone());
            }
            proof.validate_for(&prepared)?;
        } else {
            let path = args.path("--special")?;
            let hash = args.value("--special-sha256")?;
            let proof = physical::PhysicalOutcome::load_checked(path, hash, &prepared)?;
            b.external_owned_heap_bytes += proof.owned_heap_bytes()?;
            r.inputs.insert(
                "special".into(),
                Artifact {
                    path: path.canonicalize().map_err(io)?,
                    sha256: hash.into(),
                },
            );
            if p.family.starts_with("cubic-") {
                let current = physical::verify_madelung_with_observer(&prepared, &mut b, |e| {
                    j.record("madelung", &e)
                })?;
                persist(&mut r, &out, "madelung", &current)?;
                if let Some(f) = current.failure() {
                    r.failure = Some(f.clone());
                }
                if let Some(f) = current.finalization_failure() {
                    r.finalization_failure = Some(f.clone());
                }
                current.validate_success()?;
                if p.replicas == 3 {
                    let path = args.path("--madelung")?;
                    let hash = args.value("--madelung-sha256")?;
                    let left = physical::MadelungOutcome::load_checked(path, hash)?;
                    b.external_owned_heap_bytes += left.owned_heap_bytes()?;
                    let left_receipt = Receipt::load(path, hash)?;
                    physical::compare_madelung_replication(
                        &left,
                        &left_receipt.prepared,
                        &current,
                        &prepared,
                    )?;
                    r.inputs.insert(
                        "madelung64".into(),
                        Artifact {
                            path: path.canonicalize().map_err(io)?,
                            sha256: hash.into(),
                        },
                    );
                }
            }
        }
        Ok(())
    })();
    complete(
        &mut r,
        &out,
        &prepared.join("verification.json"),
        &b,
        0,
        0,
        result,
    )
}
fn sweep(args: &cli::Arguments) -> Result<(), CalibrationFailure> {
    let prepared = prepared_paths(args)?;
    let prior_path = prepared.join("verification.json");
    let prior = Receipt::load(&prior_path, args.value("--verification-sha256")?)?;
    prior.require_success()?;
    if prior.phase != "verify" || prior.prepared != prepared {
        return Err(invalid("sweep requires current verification"));
    }
    let (p, snapshots, _) = validate_prepared(&prepared)?;
    runtime(&p)?;
    let mut b = prior.usage.resume(&p, parent_elapsed(args)?)?;
    if p.family == "special-3" {
        physical::PhysicalOutcome::load_checked(
            &prior_path,
            args.value("--verification-sha256")?,
            &prepared,
        )?
        .validate_special_resume(&prepared, &b)?;
    }
    let out = args.path("--out")?;
    std::fs::create_dir(out).map_err(io)?;
    let mut r = fresh_receipt(
        "sweep",
        &prepared,
        &b,
        Some(Artifact::capture(&prior_path)?),
    )?;
    let mut j = Journal::create(&out.join("events.jsonl"), r.manifest_sha256.clone())?;
    let mut timing = TimingBudget {
        calls: prior.usage.timing_calls,
    };
    let mut certificates = Vec::new();
    let mut reports = Vec::new();
    let mut baseline: Option<CandidateOutcome> = None;
    let result = (|| {
        // Complete selection + confirmation + explicit baseline + diagnostic
        // work bound precedes the first selection reference, without reducing grid.
        let outer = candidates::protocol_bytes(&p)
            + candidates::bytes(&snapshots)
            + snapshots
                .iter()
                .map(candidates::snapshot_bytes)
                .sum::<u64>()
            + prior.owned_heap_bytes();
        b.external_owned_heap_bytes = outer + r.owned_heap_bytes() + j.owned_heap_bytes();
        preflight_selection_sweep(&p, &snapshots, &b)?;
        let refs = begin_selection_with_reference_observer(&prepared, &mut b, |e| {
            j.record("reference", &e)
        })?;
        persist(&mut r, out, "selection-references", &refs)?;
        if let Some(f) = refs.failure() {
            return Err(f.clone());
        }
        let grid = candidate_grid(&p, snapshots[0].cell)?;
        for candidate in &grid {
            b.external_owned_heap_bytes = outer
                + r.owned_heap_bytes()
                + j.owned_heap_bytes()
                + cost::reports_owned_heap(&reports)?
                + cost::reports_owned_heap(&grid)?;
            let evaluated = refs.evaluate_candidate_with_observer(
                &prepared,
                candidate,
                &mut b,
                &mut certificates,
                |e| j.record("candidate", &e),
            )?;
            let failure = evaluated
                .failure
                .clone()
                .or_else(|| evaluated.finalization_failure.clone());
            reports.push(evaluated.report);
            if let Some(f) = failure {
                return Err(f);
            }
        }
        baseline = Some(
            refs.evaluate_baseline(&prepared, &mut b, &mut certificates, |e| {
                j.record("baseline-selection", &e)
            })?,
        );
        if let Some(f) = baseline
            .as_ref()
            .and_then(|o| o.failure.clone().or_else(|| o.finalization_failure.clone()))
        {
            return Err(f);
        }
        // Before accepting in-memory raw reports for selection, bind them to
        // the durable current-epoch journal and re-read each completed outcome.
        b.external_owned_heap_bytes = outer
            + r.owned_heap_bytes()
            + j.owned_heap_bytes()
            + cost::reports_owned_heap(&reports)?
            + cost::reports_owned_heap(&grid)?
            + refs.owned_heap_bytes()?
            + baseline
                .as_ref()
                .map(|o| o.owned_heap_bytes())
                .transpose()?
                .unwrap_or(0)
            + candidates::bytes(&certificates)
            + certificates
                .iter()
                .map(|c| c.owned_heap_bytes())
                .collect::<Result<Vec<_>, _>>()?
                .into_iter()
                .sum::<u64>();
        let events_path = out.join("events.jsonl");
        let ledger = Artifact::capture(&events_path)?;
        validate_report_origin(&reports, &ledger, &r.manifest_sha256)?;
        measure_feasible_candidates_with_observer(
            &prepared,
            &mut reports,
            &mut timing,
            &mut b,
            |id, e| {
                j.record(
                    "ordinary",
                    &serde_json::json!({"candidate_id": id, "event": e}),
                )
            },
        )?;
        let first = snapshots
            .iter()
            .find(|s| s.role == FixtureRole::Selection)
            .ok_or_else(|| invalid("representative missing"))?;
        let baseline_parameters = candidates::baseline_report(&p, first.cell)?.parameters;
        let measured = measure_ordinary_with_observer(
            first,
            baseline_parameters,
            &p,
            &mut timing,
            &mut b,
            |e| j.record("baseline-ordinary", &e),
        )?;
        persist(&mut r, out, "baseline-cost", &measured)?;
        if let Some(o) = &mut baseline {
            o.report.timing = measured.timing;
        }
        if let Some(f) = measured.failure.or(measured.finalization_failure) {
            return Err(f);
        }
        let frozen = freeze_prepared_selection(&prepared, &reports)?;
        let path = out.join("frozen.json");
        frozen.write_new(&path)?;
        r.artifacts
            .insert("frozen".into(), Artifact::capture(&path)?);
        persist(&mut r, out, "decisions", &frozen.decisions())?;
        Ok(())
    })();
    // These retained records are written on ordinary errors as well. The durable
    // journal remains authoritative after immediate allocator/process exits.
    if let Err(f) = &result {
        record_failure(&mut r, f.clone());
    }
    let retained = (|| {
        persist(&mut r, out, "candidates", &reports)?;
        persist(&mut r, out, "baseline-selection", &baseline)?;
        persist(&mut r, out, "certificates", &certificates)
    })();
    if let Err(f) = retained {
        record_failure(&mut r, f);
    }
    let decision_failure = if result.is_ok() && r.artifacts.contains_key("decisions") {
        let ds: Vec<TargetDecision> = journal::read_json(&r.artifact("decisions")?.path)?;
        ds.iter()
            .find(|d| d.selected_id.is_none())
            .map(|d| CalibrationFailure::TargetNotReached(d.target))
    } else {
        None
    };
    complete(
        &mut r,
        out,
        &out.join("summary.json"),
        &b,
        timing.calls,
        0,
        result,
    )?;
    if let Some(f) = decision_failure {
        Err(f)
    } else {
        Ok(())
    }
}
fn validate_report_origin(
    reports: &[CandidateReport],
    ledger: &Artifact,
    epoch: &str,
) -> Result<(), CalibrationFailure> {
    let mut completed = 0;
    journal::visit_events(&ledger.path, &ledger.sha256, epoch, |event| {
        if event.kind != "candidate" {
            return Ok(());
        }
        if let Some(value) = event.data.get("Complete").and_then(|v| v.get("outcome")) {
            let outcome: CandidateOutcome =
                serde_json::from_value(value.clone()).map_err(|e| invalid(e.to_string()))?;
            if outcome.failure.is_some() || outcome.finalization_failure.is_some() {
                return Err(invalid("failed journal candidate cannot freeze"));
            }
            if reports.get(completed) != Some(&outcome.report) {
                return Err(invalid(
                    "raw candidates differ from trusted current journal",
                ));
            }
            completed += 1;
        }
        Ok(())
    })?;
    if completed != reports.len() {
        return Err(invalid("incomplete raw candidate journal"));
    }
    Ok(())
}
fn load_campaign(
    args: &cli::Arguments,
    phase: &str,
) -> Result<(Receipt, Protocol, CandidateBudget), CalibrationFailure> {
    let r = Receipt::load(args.path("--campaign")?, args.value("--expected-sha256")?)?;
    r.require_success()?;
    if r.phase != phase {
        return Err(invalid(format!("requires {phase} decision manifest")));
    }
    let (p, _, _) = validate_prepared(&r.prepared)?;
    runtime(&p)?;
    let b = r.usage.resume(&p, parent_elapsed(args)?)?;
    Ok((r, p, b))
}
fn confirm(args: &cli::Arguments) -> Result<(), CalibrationFailure> {
    let (prior, p, mut b) = load_campaign(args, "sweep")?;
    let prepared = &prior.prepared;
    let frozen_artifact = prior.artifact("frozen")?;
    let frozen = FrozenSelection::load(prepared, &frozen_artifact.path, &frozen_artifact.sha256)?;
    let mut certificates: Vec<TranslationCertificate> =
        journal::read_json(&prior.artifact("certificates")?.path)?;
    let out = args.path("--out")?;
    std::fs::create_dir(out).map_err(io)?;
    let mut r = fresh_receipt(
        "confirm",
        prepared,
        &b,
        Some(Artifact::capture(args.path("--campaign")?)?),
    )?;
    let j = RefCell::new(Journal::create(
        &out.join("events.jsonl"),
        r.manifest_sha256.clone(),
    )?);
    let mut target_failure = None;
    let result = (|| {
        b.external_owned_heap_bytes = prior.owned_heap_bytes()
            + r.owned_heap_bytes()
            + j.borrow().owned_heap_bytes()
            + candidates::protocol_bytes(&p);
        let mut confirmation = frozen.confirm_with_observers(
            prepared,
            &mut b,
            &mut certificates,
            |e| j.borrow_mut().record("candidate", &e),
            |e| j.borrow_mut().record("reference", &e),
        )?;
        if confirmation.failure.is_none() {
            // Baseline holdout is after BOTH IDs have been frozen, even when
            // no target is reachable. It reuses the same reference row engine.
            if confirmation.references.is_none() {
                confirmation.references = Some(cost::generate_references_with_observer(
                    prepared,
                    FixtureRole::Confirmation,
                    &mut b,
                    |e| j.borrow_mut().record("reference", &e),
                )?);
            }
            b.external_owned_heap_bytes = prior.owned_heap_bytes()
                + r.owned_heap_bytes()
                + j.borrow().owned_heap_bytes()
                + candidates::protocol_bytes(&p)
                + frozen.owned_heap_bytes()?
                + confirmation
                    .candidates
                    .iter()
                    .map(|o| o.owned_heap_bytes())
                    .collect::<Result<Vec<_>, _>>()?
                    .into_iter()
                    .sum::<u64>();
            let refs = confirmation
                .references
                .as_ref()
                .ok_or_else(|| invalid("confirmation refs"))?;
            if let Some(f) = refs.failure() {
                confirmation.failure = Some(f.clone());
            } else {
                let baseline =
                    refs.evaluate_baseline(prepared, &mut b, &mut certificates, |e| {
                        j.borrow_mut().record("baseline-confirmation", &e)
                    })?;
                confirmation.failure = baseline
                    .failure
                    .clone()
                    .or_else(|| baseline.finalization_failure.clone());
                persist(&mut r, out, "baseline-confirmation", &baseline)?;
            }
        }
        persist(&mut r, out, "confirmation", &confirmation)?;
        persist(&mut r, out, "decisions", &confirmation.decisions)?;
        persist(&mut r, out, "certificates", &certificates)?;
        if let Some(f) = &confirmation.failure {
            return Err(f.clone());
        }
        for d in &confirmation.decisions {
            if let Some(fragment) = frozen.fragment(&confirmation, d.target, &p)? {
                let candidate = confirmation
                    .candidates
                    .iter()
                    .find(|c| Some(c.report.id) == d.selected_id)
                    .ok_or_else(|| invalid("fragment confirmed candidate"))?;
                validate_fragment_schema(&fragment, &candidate.report.parameters)?;
                let annotated = format!("# report: {}\n# prepared manifest SHA256: {}\n# source/protocol/snapshot identities: {}/manifest.json\n{fragment}", out.join("summary.json").display(), r.manifest_sha256, prepared.display());
                super::write_new(
                    &out.join(format!("selected-{}.toml", d.target)),
                    annotated.as_bytes(),
                )?;
            } else if target_failure.is_none() {
                target_failure = Some(if d.selected_id.is_none() {
                    CalibrationFailure::TargetNotReached(d.target)
                } else {
                    CalibrationFailure::ConfirmationFailed(d.target)
                });
            }
        }
        Ok(())
    })();
    complete(
        &mut r,
        out,
        &out.join("summary.json"),
        &b,
        prior.usage.timing_calls,
        0,
        result,
    )?;
    target_failure.map_or(Ok(()), Err)
}
fn profile(args: &cli::Arguments) -> Result<(), CalibrationFailure> {
    let (prior, _, mut b) = load_campaign(args, "confirm")?;
    let prepared = &prior.prepared;
    let selection_ref = prior
        .prior
        .as_ref()
        .ok_or_else(|| invalid("missing selection chain"))?;
    let selection = Receipt::load(&selection_ref.path, &selection_ref.sha256)?;
    let frozen_artifact = selection.artifact("frozen")?;
    let frozen = FrozenSelection::load(prepared, &frozen_artifact.path, &frozen_artifact.sha256)?;
    let confirmation_artifact = prior.artifact("confirmation")?;
    let confirmation = ConfirmationOutcome::load_checked(
        prepared,
        &confirmation_artifact.path,
        &confirmation_artifact.sha256,
        &mut b,
    )?;
    let out = args.path("--out")?;
    let mut writer = profile::ProfileWriter::create(out)?;
    let mut r = fresh_receipt(
        "profile",
        prepared,
        &b,
        Some(Artifact::capture(args.path("--campaign")?)?),
    )?;
    let mut j = Journal::create(&out.join("events.jsonl"), r.manifest_sha256.clone())?;
    let mut timing = profile::ProfileBudget {
        calls: prior.usage.profile_calls,
    };
    let result = (|| {
        b.external_owned_heap_bytes = prior.owned_heap_bytes()
            + selection.owned_heap_bytes()
            + r.owned_heap_bytes()
            + j.owned_heap_bytes()
            + frozen.owned_heap_bytes()?
            + confirmation.owned_heap_bytes()?;
        let mut models = profile::confirmed_profile_models(prepared, &frozen, &confirmation)?;
        let mut baseline = profile::validate_baseline_with_observer(prepared, &mut b, |e| {
            j.record("baseline-functional", &e)
        })?;
        persist(&mut r, out, "baseline-functional", &baseline)?;
        if let Some(f) = baseline.failure {
            return Err(f);
        }
        models.insert(
            0,
            baseline
                .model
                .take()
                .ok_or_else(|| invalid("baseline controls did not authorize model"))?,
        );
        drop(baseline);
        persist(&mut r, out, "profile-models", &models)?;
        b.external_owned_heap_bytes +=
            candidates::bytes(&models) + models.iter().map(|m| m.owned_heap_bytes()).sum::<u64>();
        let measured = profile::measure_profiles_with_observer(
            prepared,
            &models,
            &mut timing,
            &mut b,
            |event| {
                j.record("profile", &event)?;
                if let profile::ProfileEvent::Completed(sample) = event {
                    writer.completed(sample)?;
                }
                Ok(())
            },
        )?;
        persist(&mut r, out, "profile-outcome", &measured)?;
        if let Some(f) = measured.failure {
            return Err(f);
        }
        Ok(())
    })();
    for (key, file) in [
        ("profiles", "phase-times.jsonl"),
        ("profile-csv", "phase-times.csv"),
    ] {
        r.artifacts
            .insert(key.into(), Artifact::capture(&out.join(file))?);
    }
    complete(
        &mut r,
        out,
        &out.join("summary.json"),
        &b,
        prior.usage.timing_calls,
        timing.calls,
        result,
    )
}

/// Preregistered test-only diagnostic. It never freezes a partial grid or emits
/// selected parameters. The production CLI remains exactly six modes.
#[cfg(test)]
pub fn registered_special_smoke() -> Result<(), CalibrationFailure> {
    let env = |key| std::env::var(key).map_err(|_| invalid(format!("missing {key}")));
    let previous = PathBuf::from(env("MODENA_B05_TASK7B_VERIFY")?);
    let expected = env("MODENA_B05_TASK7B_VERIFY_SHA256")?;
    let prior = Receipt::load(&previous, &expected)?;
    prior.require_success()?;
    let prepared = &prior.prepared;
    let (p, snapshots, _) = validate_prepared(prepared)?;
    if p.family != "special-3" || p.reference_ladder.len() != 6 || p.threads != 1 {
        return Err(invalid(
            "registered special six-level serial diagnostic only",
        ));
    }
    runtime(&p)?;
    let elapsed: u64 = env("MODENA_B05_TASK7B_ACTIVE_NS")?
        .parse()
        .map_err(|_| invalid("smoke elapsed"))?;
    let mut b = prior.usage.resume(&p, elapsed)?;
    physical::PhysicalOutcome::load_checked(&previous, &expected, prepared)?
        .validate_special_resume(prepared, &b)?;
    let out = PathBuf::from(env("MODENA_B05_TASK7B_OUT")?);
    let mut writer = profile::ProfileWriter::create(&out)?;
    let mut receipt = fresh_receipt(
        "diagnostic-smoke",
        prepared,
        &b,
        Some(Artifact::capture(&previous)?),
    )?;
    let mut journal = Journal::create(&out.join("events.jsonl"), receipt.manifest_sha256.clone())?;
    let mut timing = TimingBudget {
        calls: prior.usage.timing_calls,
    };
    let mut profiling = profile::ProfileBudget {
        calls: prior.usage.profile_calls,
    };
    let result = (|| {
        b.external_owned_heap_bytes = candidates::protocol_bytes(&p)
            + candidates::bytes(&snapshots)
            + snapshots
                .iter()
                .map(candidates::snapshot_bytes)
                .sum::<u64>()
            + receipt.owned_heap_bytes()
            + prior.owned_heap_bytes()
            + journal.owned_heap_bytes();
        let refs = begin_selection_with_reference_observer(prepared, &mut b, |e| {
            journal.record("reference", &e)
        })?;
        persist(&mut receipt, &out, "selection-references", &refs)?;
        if let Some(f) = refs.failure() {
            return Err(f.clone());
        }
        let mut certificates = Vec::new();
        let mut baseline = refs.evaluate_baseline(prepared, &mut b, &mut certificates, |e| {
            journal.record("baseline-selection", &e)
        })?;
        if let Some(f) = baseline
            .failure
            .as_ref()
            .or(baseline.finalization_failure.as_ref())
        {
            persist(&mut receipt, &out, "baseline-selection", &Some(&baseline))?;
            return Err(f.clone());
        }
        let first = snapshots
            .iter()
            .find(|s| s.role == FixtureRole::Selection)
            .ok_or_else(|| invalid("smoke representative"))?;
        b.external_owned_heap_bytes += refs.owned_heap_bytes()?
            + baseline.owned_heap_bytes()?
            + candidates::bytes(&certificates)
            + certificates
                .iter()
                .map(|c| c.owned_heap_bytes())
                .collect::<Result<Vec<_>, _>>()?
                .into_iter()
                .sum::<u64>();
        let costs = measure_ordinary_with_observer(
            first,
            baseline.report.parameters.clone(),
            &p,
            &mut timing,
            &mut b,
            |e| journal.record("baseline-ordinary", &e),
        )?;
        baseline.report.timing = costs.timing.clone();
        persist(&mut receipt, &out, "baseline-selection", &Some(&baseline))?;
        persist(&mut receipt, &out, "baseline-cost", &costs)?;
        persist(&mut receipt, &out, "certificates", &certificates)?;
        if let Some(f) = costs.failure.or(costs.finalization_failure) {
            return Err(f);
        }
        let mut validated = profile::validate_baseline_with_observer(prepared, &mut b, |e| {
            journal.record("baseline-functional", &e)
        })?;
        persist(&mut receipt, &out, "baseline-functional", &validated)?;
        if let Some(f) = validated.failure {
            return Err(f);
        }
        let model = validated
            .model
            .take()
            .ok_or_else(|| invalid("smoke baseline model"))?;
        drop(validated);
        let profiled = profile::measure_profiles_with_observer(
            prepared,
            &[model],
            &mut profiling,
            &mut b,
            |e| {
                journal.record("profile", &e)?;
                if let profile::ProfileEvent::Completed(sample) = e {
                    writer.completed(sample)?;
                }
                Ok(())
            },
        )?;
        persist(&mut receipt, &out, "profile-outcome", &profiled)?;
        if let Some(f) = profiled.failure {
            return Err(f);
        }
        // Ordinary timings are unchanged by the diagnostic phase.
        if baseline.report.timing != costs.timing {
            return Err(invalid("profile rewrote primary timing"));
        }
        Ok(())
    })();
    complete(
        &mut receipt,
        &out,
        &out.join("summary.json"),
        &b,
        timing.calls,
        profiling.calls,
        result,
    )
}

#[cfg(test)]
mod integration_tests {
    use super::*;
    #[test]
    fn phase_primary_and_finalization_failures_remain_inspectable() {
        let p = frozen_protocol("special", 2, 1, true).unwrap();
        let b = CandidateBudget::new(&p.limits, ReferenceBudget::new(&p.limits));
        let mut r = Receipt {
            schema: "modena-spme-phase-v1".into(),
            phase: "sweep".into(),
            prepared: "/tmp".into(),
            manifest_sha256: "epoch".into(),
            prior: None,
            inputs: BTreeMap::new(),
            artifacts: BTreeMap::new(),
            usage: Usage::capture(&b, 0, 0),
            failure: None,
            finalization_failure: None,
        };
        record_failure(&mut r, CalibrationFailure::BudgetExceeded("primary".into()));
        let before = r.owned_heap_bytes();
        record_failure(&mut r, io("first finalization".repeat(100)));
        record_failure(&mut r, io("second finalization"));
        assert!(matches!(
            r.failure,
            Some(CalibrationFailure::BudgetExceeded(_))
        ));
        let text = r.finalization_failure.as_ref().unwrap().to_string();
        assert!(text.contains("first finalization") && text.contains("second finalization"));
        assert!(r.owned_heap_bytes() > before + 1000);
    }
    #[test]
    fn checked_special_import_rejects_status_without_controls() {
        let tmp = tempfile::tempdir().unwrap();
        let prepared = tmp.path().join("prepared");
        let p = frozen_protocol("special", 2, 1, true).unwrap();
        prepare(&prepared, p.clone()).unwrap();
        let mut b = CandidateBudget::new(&p.limits, ReferenceBudget::new(&p.limits));
        let proof = physical::verify_special_with_observer(&prepared, &mut b, |_| {
            Err(io("before first call"))
        })
        .unwrap();
        assert_eq!(b.reference.calls, 0);
        let mut value = serde_json::to_value(proof).unwrap();
        value["failure"] = serde_json::Value::Null;
        value["finalization_failure"] = serde_json::Value::Null;
        value["completed"] = serde_json::json!(true);
        let mut r = fresh_receipt("verify", &prepared, &b, None).unwrap();
        persist(&mut r, tmp.path(), "physical", &value).unwrap();
        let path = tmp.path().join("receipt.json");
        let authority = r.write_new(&path).unwrap();
        assert!(
            physical::PhysicalOutcome::load_checked(&path, &authority.sha256, &prepared).is_err()
        );
    }
}

# B05 SPME Empirical Calibration Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use `subagent-driven-development` after controller handoff. Delegation is authorized: execute Tasks 1–8 sequentially, with one implementer per task followed by independent review before advancing. The controller performs preflight/scaffolding directly, without a separate agent or review task. Implementation checkbox steps are future work; controller preflight and the fixed contract are adopted before Task1.

**Goal:** build and evaluate a bounded offline SPME selector with independently converged Ewald reference, two empirical force targets, held-out confirmation, retained cost/error evidence, and actual opt-in phase profiling.

**Architecture:** add a `modena-core` example and example-only support modules; preserve existing explicit production electrostatics configuration. Evaluate the existing analytic order-4 model against an independent all-image/full-±k reference. Share the existing SPME implementation through a compile-time no-op/profile observer so ordinary primary timing and opt-in diagnostic timing remain separate.

**Tech Stack:** existing Rust std/serde/serde_json/toml/libm/rayon/rustfft, existing core public model API, host shasum (agent commands use RTK) and macOS `/usr/bin/time -l`, Python standard library for the sequential campaign launcher. No new Rust dependency.

**Spec:** `docs/scientific-hardening/B05-SPME-PROTOCOL.md`. Binding sections2–4 and7 fix all scientific/serialized/profiling values; per-task briefs identify the applicable sections. The controller adopts this specification and plan before implementation.

## Global Constraints

- B03 is delivered. Controller may hand off independent example-only Tasks 1–3 on the verified B03 base while B04 proceeds in a separate checkout. Reconcile the unpublished B05 base with final reviewed B04 before Task 4 candidate evaluation or production profiling. No concurrent work within a B05 task sequence; no competing builds/tests/writes during another case's exclusive measurement lease.
- Preserve explicit `EwaldSpec` parameters, current defaults, unit conversion, runner construction, analytic order-4 B-spline differentiation, existing FFT-plan/workspace reuse, additive ForceModel semantics, and charged-mesh COM behavior.
- Never attribute a P3M/ik estimator or a universal error bound to this analytic SPME implementation. No hidden force-mean subtraction.
- Support only reduced units, finite orthorhombic cells with 3D periodicity, rc < half the shortest cell edge, and the exact frozen fixture families. Reject unsupported geometry/units and unvalidated inputs.
- Force targets are empirical; independent all-image real-space/full-positive-and-negative-k Ewald convergence and held-out confirmation are mandatory before emitting selected configurations.
- Absolute RMS, maximum atom-vector error, energy/atom and scalar isotropic-strain W/atom gates accompany relative RMS. Current Virial is scalar: B05 does not implement or certify the six-component stress tensor; that belongs to B07.
- Freeze fixtures, tolerance, operation/memory budgets and timing protocol before candidate measurements; preserve failed rows and raw data.
- Actual spread/FFT/gather/real-space phase measurements are required for full B05 completion. Their observer is opt-in per call, and primary candidate timing uses the ordinary uninstrumented compute path. Staging the offline-helper alternative may leave profiling explicitly pending; it cannot close full B05.
- No new Rust dependency, production autotuning, interpolation-order change, P3M replacement, tensor extension, automatic input rewrite, or unbounded refinement.
- All future shell commands use rtk; no competing timing campaign or source mutation during another stream's exclusive measurement lease.
- New branches use the user-required `ricardo/` prefix, overriding the host default `ricardoleal20/`.
- The controller owns preflight/scaffolding; Tasks 1–8 use authorized sequential subagent-driven development with one implementer per task and independent review before advancing.

## File ownership and interfaces

Create only during execution:

| File | Responsibility |
|---|---|
| `docs/scientific-hardening/B05-SPME-PROTOCOL.md` | versioned frozen spec, family identities, tolerances, limits, statuses and commands |
| `docs/superpowers/plans/2026-11-06-b05-spme-calibration.md` | executed checkbox plan copied from this draft and adjusted only by recorded decisions |
| `crates/modena-core/examples/spme_convergence.rs` | strict manual CLI, phase orchestration, files, hashes, process exits |
| `crates/modena-core/examples/spme_calibration/mod.rs` | exact serialized records from spec 7.2, validation, bounded family orchestration, reusable evaluator |
| `crates/modena-core/examples/spme_calibration/fixtures.rs` | deterministic selection/confirmation/regression-only snapshot arrays |
| `crates/modena-core/examples/spme_calibration/reference.rs` | independent full image/full ±k reference, its term counts, ladder and convergence status |
| `crates/modena-core/examples/spme_calibration/metrics.rs` | raw vector RMS, other gates, deduplication, deterministic cost selector, CSV/SVG |
| `crates/modena-core/tests/spme_calibration.rs` | import example modules by path; pure metric, oracle, selector, provenance and CLI regressions |
| `crates/modena-core/tests/spme_profile.rs` | ordinary/profile U/F/W parity, additive semantics, phase labels and cache reuse |
| `bench/spme_calibration_campaign.py` | sequential frozen family/threads campaign, `/usr/bin/time -l` peak RSS capture, raw stdout/stderr and environment manifests |

Modify only `crates/modena-core/src/spme.rs` for observer/profile interface and stale analytic-force comments, and `crates/modena-core/tests/electrostatics.rs` for its stale `ik` explanatory comment. No planned production input/runner/Ewald/Cargo changes. No source edit to B03 files. Other existing tests are run unchanged.

Example/test interfaces, all errors `CalibrationFailure`:

```rust
pub fn reference(snapshot: &Snapshot, settings: ReferenceSettings,
                 limits: &Limits) -> Result<Evaluation, CalibrationFailure>;
pub fn converged_reference(snapshot: &Snapshot, protocol: &Protocol)
    -> Result<ReferenceRecord, CalibrationFailure>;
pub fn evaluate_spme(snapshot: &Snapshot, p: Parameters)
    -> Result<Evaluation, CalibrationFailure>;
pub fn metrics(reference: &Evaluation, candidate: &Evaluation,
               n: usize, force_floor: f64) -> Result<ErrorMetrics, CalibrationFailure>;
pub fn campaign(protocol: &Protocol, snapshots: &[Snapshot])
    -> Result<CampaignReport, CalibrationFailure>;
pub fn select(rows: &[CandidateReport], target: f64)
    -> Result<Parameters, CalibrationFailure>;

// Fields are private; types are existing modena_core types.
pub struct SpmeEvaluator {
    model: SpmeForceModel, system: System, neighbors: NeighborList,
}
impl SpmeEvaluator {
    pub fn new(snapshot: &Snapshot, parameters: Parameters)
        -> Result<Self, CalibrationFailure>;
    pub fn evaluate(&mut self) -> Result<Evaluation, CalibrationFailure>;
    pub fn ordinary_batch(&mut self, calls: usize) -> Result<u64, CalibrationFailure>;
}
```

The support types are exactly those in spec 6 and 7.2; include every field, enum variant, finite-value check and unit. Integration test import:

```rust
#[path = "../examples/spme_calibration/mod.rs"]
mod calibration;
```

Core additive interface is exactly `SpmeForceModel::compute_profiled(&self,&mut System,&mut Virial,&mut NeighborList)->modena_core::Result<(f64,SpmeProfile)>` and the twelve Duration fields in spec 7.4. Keep the ordinary `ForceModel` trait signature unchanged. No full tensor type is introduced.

## Controller preflight: Bind the independent reference base and final B04 integration gate

**Files:** read current `AGENTS.md`/RTK instructions, `crates/modena-core/src/spme.rs`, `src/ewald.rs`, `src/io/input.rs`, `src/runner.rs`, `tests/pair_scales.rs`, `crates/modena-bench/Cargo.toml`, and `crates/modena-bench/src/lib.rs` if present after B03. Only the final plan/spec documents may be edited in this task.

**Owner:** controller directly; no implementer agent or separate review dispatch. **Consumes:** completed B03/source identity and released measurement lease; final reviewed B04 revision is required before Task4 and numerical campaigns. **Produces:** recorded exact source hashes, final ownership table, verified API compatibility and unchanged physical contract for Tasks 1–8.

- [ ] Confirm controller handoff and absence of another source/timing lease before any repository mutation or compute.
- [ ] Use the user-required `ricardo/` prefix for any new implementation branch; do not use the host default `ricardoleal20/`.
- [ ] Read the listed stable source and benchmark public interfaces once. B03 is expected to have added a benchmark library plus clap/rayon/serde_json dependencies. Record whether reusable timing/provenance conventions exist; consume them only if they are public and do not create a core->bench dependency or edit B03-owned files.
- [ ] Compare all interfaces in this plan with the actual base. If names/signatures changed, update only call-site adapters and record the change; if the physical model/virial/background/pair policy changed, refreeze spec and fixtures before implementation measurements. Do not assume this preparation's old line numbers or source hashes are current.
- [ ] Create final repository spec/plan from these drafts, retaining `Adopted contract, no B05 measurements` status. Freeze physical protocol before any candidate evaluation.
- [ ] Commit only the two documents: `:memo: Add: B05 SPME calibration specification and plan`.

This is controller preflight/scaffolding rather than a separate implementation task; it has no agent/review cycle. It is not permission to repeat broad research, implement B03, or redesign benchmark architecture. Tasks 1–8 retain their task boundaries and each receives one sequential implementer followed by independent review under `subagent-driven-development`.

## Task 1: Freeze types, fixtures and pure error gates

**Files:** create example entry point/support mod/fixtures/metrics and `tests/spme_calibration.rs`; update frozen spec with actual materialized fixture hashes.

**Consumes:** spec sections 2, 4 and 7.2. **Produces:** schema `modena-spme-calibration-v1`, fixed Snapshot/Protocol values, strict validation, pure `metrics`, target gate and `prepare` command. No force evaluation yet.

- [ ] Write failing metric and validation tests with the exact vector-RMS denominator:

```rust
#[test]
fn vector_rms_denominator_and_raw_net_force_are_explicit() {
    let r = Evaluation { energy: 0.0, strain_w: 0.0,
        forces: vec![[3.0,4.0,0.0],[0.0,0.0,0.0]] };
    let c = Evaluation { energy: 0.0, strain_w: 0.0,
        forces: vec![[3.0,4.0,0.0],[0.0,0.0,1.0]] };
    let m = metrics(&r,&c,2,1e-8).unwrap();
    assert!((m.reference_force_rms-12.5_f64.sqrt()).abs()<1e-14);
    assert!((m.force_rms_abs-0.5_f64.sqrt()).abs()<1e-14);
    assert!((m.force_rms_rel.unwrap()-0.2).abs()<1e-14);
    assert_eq!(m.max_atom_force_error,1.0);
    assert_eq!(m.net_force,[3.0,4.0,1.0]);
}
```

- [ ] Add tests `relative_zero_denominator_is_undefined`, `force_mean_is_not_removed`, `all_candidate_side_gates_are_required`, `finite_dimensions_charges_and_scales_are_validated`, `unsupported_units_pbc_or_triclinic_rejected`, `allocation_and_term_caps_are_checked_before_work`, `fixture_materialization_is_byte_stable`, `holdout_ids_are_disjoint`.
- [ ] Implement strict serializable records; canonicalize pairs by min/max index; accept exact duplicate identical policy once, reject conflicting duplicates, invalid indices or self pairs; scales must lie in [0,1]. Validate minimum separation >=0.2*a and nonempty/nonzero charges. Snapshot hashing occurs after deterministic wrapping and exact materialization.
- [ ] Implement fixture generator `z=((i+1)*1664525+(axis+1)*1013904223+seed*69069)%65536`, displacement `amplitude*a*(2*z/65535-1)` in u64/f64 as fixed in spec. Neutral selection seeds/amplitudes 7/.02,19/.05,31/.08; neutral confirmation 43/.035,59/.065. Add background +0.5 charge fixture, special rectangular neutral/+0.7/scaled fixtures, anisotropy, full box-wrap regressions, unperturbed near-zero-D crystal.
- [ ] Implement prepare family/replicas/threads/extended flags; exact one-family materialized JSON with `FixtureRole`; invalid special+replicas3 is rejected. Declare units reduced/all-three PBC; source/protocol/coordinate hashes are checked before all later phases.
- [ ] Run `rtk cargo test -p modena-core --test spme_calibration`; first observe missing-interface failures, then passing pure tests. Confirm prepare performs no force call using a test-only counter or stub evaluator.
- [ ] Commit owned files: `:white_check_mark: Add: frozen SPME fixtures and error gates`.

## Task 2: Implement and independently validate the full Ewald oracle

**Files:** `reference.rs`, records/status handling in mod, integration tests. **Consumes:** Snapshot/ReferenceSettings/Limits and exact equations in spec 3. **Produces:** `reference`, streamed full-±k/all-image Evaluation, term counts, independently validated physical invariants.

- [ ] Write failing tests `reference_force_is_negative_energy_gradient`, `reference_w_is_isotropic_strain_derivative`, `full_and_half_k_oracles_agree`, `self_image_energy_is_included`, `charged_background_w_is_three_u`, `bare_pair_removal_inside_and_outside_cutoff`, `reference_rigid_translation_invariant`, `reference_full_box_wrap_invariant`.
- [ ] Code the independent image loop and full ±k direct sum exactly from spec. Stream terms, compensated sum every energy/W/force component, omit only i=j,n=0, include periodic self-image energy/W and uniform background, independently remove bare minimum-image special-pair fractions. Reuse only core `System` representation at the adapter; do not call core Ewald/SPME, CoulombPolicy, exclusion_correction, NeighborList or SimBox::min_image in the oracle.
- [ ] Before each loop, checked-estimate terms, vector capacities, deadlines and accumulated family work. Real bounds `ceil(R/L_axis)+1`; reciprocal bounds `ceil(K*L_axis/(2*pi))`; spherical filter after enumeration. Check deadline every <=1024 terms. Fail with BudgetExceeded instead of partially returning zero.
- [ ] On the 3-ion neutral, charged and scaled rectangular fixtures, compare all 9 force components to central energy FD at h=1e-5*a and h/2; require residual <=`1e-7*F_unit+1e-6*max(abs(F),abs(FD))`. Test W at lambda=1±1e-5 and half step; require <=`1e-7*E_unit+1e-6*max(abs(W),abs(FD))`. Membership in every spherical cutoff must be unchanged during these derivatives. Special-pair removal <=1e-11 scaled units; background algebra <=1e-13*E_unit; full/half ±k <=1e-10 scaled units.
- [ ] Add Madelung refinement check with existing approximate alpha=1.747558775, r0=a/2 and `abs(E_per_pair+alpha/r0)<=2e-5*E_unit` after the strict oracle convergence passes; this is a sign/factor sanity check only. Require oracle integer-wrap/rigid-translation E/F/W agreement <=1e-9 scaled units. Add a secondary core-direct-Ewald comparison only where independent image/tail convergence is demonstrated, to the same reference budgets; never use matching kcuts as the accuracy reference.
- [ ] Run `rtk cargo test -p modena-core --test spme_calibration reference_`; all oracle tests pass before candidate measurements.
- [ ] Commit `:sparkles: Add: independent all-image Ewald oracle`.

## Task 3: Bound reference refinement and handle precision limits

**Files:** reference.rs/mod/tests. **Consumes:** validated reference and frozen protocol. **Produces:** `converged_reference` and ReferenceRecord with every level/crosscheck/delta/status.

- [ ] Write deterministic tests using a small injected reference function returning prechosen Evaluation sequences: `two_consecutive_reference_passes_required`, `cross_kappa_disagreement_blocks_calibration`, `stagnation_above_budget_is_explicit`, `nonfinite_reference_is_precision_limited`, `extended_reference_is_chosen_before_run`, `reference_deadline_preserves_failure_record`. The injected function is a test-only argument to the ladder runner, not a production physics shortcut.
- [ ] Implement exactly four preset `(R/a,K*a)` levels `(4,8),(5,10),(6,12),(7,14)` at kappa*a=1.1; extended protocol appends `(8,16),(9,18)` only if frozen during prepare. Independently check final settings at kappa*a=.9 and1.3. Never expand this ladder automatically.
- [ ] Compute force vector RMS/max, E/N and W/N refinement deltas. Force RMS reference budget `min(1e-8*F_unit,.01*1e-4*D)` above D floor; absolute1e-8 below floor; max1e-7*F_unit; E/N,W/N1e-9*E_unit. Require two consecutive passes plus both cross-kappa passes.
- [ ] Implement normalized score/stagnation/precision rules from spec 7.3. Preserve `ReferenceUnconverged`, `ReferenceStagnated`, `ReferencePrecisionLimited` and BudgetExceeded in JSON; these block every relative-target selected fragment.
- [ ] Run `rtk cargo test -p modena-core --test spme_calibration`; confirm exact failure statuses with synthetic ladders and bounded real three-atom oracle fixture.
- [ ] Commit `:sparkles: Add: bounded Ewald reference convergence checks`.

## Task 4: Evaluate bounded SPME candidates against raw reference metrics

**Files:** reusable evaluator/mod/metrics/integration tests. **Consumes:** frozen per-family snapshots and converged references. **Produces:** complete candidate/error rows, mesh aliases and feasible target sets; no primary timing or selected config yet.

- [ ] Write tests `actual_mesh_deduplication_preserves_aliases`, `one_shift_failure_rejects_candidate`, `relative_undefined_does_not_certify_target`, `different_families_do_not_share_denominators`, `scales_are_applied_once`, `existing_explicit_model_parameters_unchanged`, `mesh_lifecycle_uses_existing_cache_keys`.
- [ ] Implement 126 raw tuples from kappa*a `[1.2,1.6,2,2.4,3,3.6]`, rc/Lmin `[.30,.40,.49]`, kcut*a `[8,12,16,24,32,48,64]`; sort and deduplicate `(kappa,rc,actual_mesh)` while retaining aliases. Do not use an error estimator to prune/approve a candidate. Checked mesh <=1048576 points; preserve current `mesh_for_box` next-5-smooth algorithm.
- [ ] Implement evaluator with persistent production SpmeForceModel and consistent model/list PairScale policies; `NeighborList::with_pair_scales(rc,.3,&scales)`; zero forces/W before each evaluation and return raw results. Charge and background conventions exactly match the reference.
- [ ] Evaluate every selection snapshot and seven fixed sub-grid shifts from spec. Compute R=A/D with per-snapshot reference vector RMS and Dfloor1e-8*F_unit. Require simultaneously R<=t,A<=t*F_unit,M<=10t*F_unit,abs(E/N)<=t*E_unit,abs(W/N)<=10t*E_unit at targets1e-3 and1e-4. Do not subtract net force; retain it and translation variation.
- [ ] Test fixed-mesh energy gradients and isotropic W. Require unchanged production mesh and cutoff memberships at the FD stencil; use fixed-mesh reciprocal helper if at an exact mesh boundary and label that situation, never infer stress from a mesh jump. Existing scalar W is the only strain quantity; do not create six-strain assertions.
- [ ] Run `rtk cargo test -p modena-core --test spme_calibration`; then existing `electrostatics`, `force_hardening`, `pair_scales`, `force_cache`, `final_integration`, `sampling_state` and `--lib spme::tests` suites once. Resolve physics failures before later candidate/timing campaigns.
- [ ] Commit `:sparkles: Add: bounded SPME error sweep`.

## Task 5: Implement ordinary cost ranking and freeze held-out decisions

**Files:** mod/metrics/example binary/tests. **Consumes:** complete feasible sets from Task4; ordinary production compute. **Produces:** raw ordinary cost records, deterministic selected candidate IDs and held-out confirmation results.

- [ ] Write pure tests `selector_requires_all_gates_and_finite_timing`, `five_percent_tie_rule_is_deterministic`, `holdout_failure_does_not_select_another_candidate`, `target_not_reached_emits_no_config`. Use synthetic rows/timing records; no speed-sensitive CI assertions.
- [ ] Implement ordinary persistent model timing; five warmups and five batches*20 calls. Include force zero/W initialization in ordinary batch timing, existing neighbor cache check and normal compute. No coordinate/force result clones in the measured region. Consume energy/W/force slice with black_box.

```rust
let started = std::time::Instant::now();
for _ in 0..calls {
    self.system.zero_forces();
    let mut virial = Virial::default();
    let energy = self.model.compute(&mut self.system,&mut virial,&mut self.neighbors)
        .map_err(|e| CalibrationFailure::InvalidInput(e.to_string()))?;
    std::hint::black_box((energy,virial.w,self.system.forces.as_slice()));
}
let elapsed_ns = u64::try_from(started.elapsed().as_nanos())
    .map_err(|_| CalibrationFailure::BudgetExceeded("timing duration overflow".into()))?;
```

- [ ] Record cold construction/first call and a forced neighbor rebuild separately from warm batches; median/IQR and all raw batches retained. Time every feasible candidate within caps, not only a mesh-size shortlist. If work/time budget prevents this, no fastest-feasible claim.
- [ ] Choose lowest ordinary median; all candidates within5% of best form a tie, then smallest actual mesh product, smallest rc, lexicographic kappa/kcut. Freeze candidate ID for each target before any held-out candidate evaluation.
- [ ] Confirm each frozen ID on every confirmation snapshot and translation using exactly the same R/A/M/E/W gates. Failure retains chosen ID with ConfirmationFailed and emits no fragment for that target; do not fall back. Success may emit reduced-unit explicit fragment labeled empirical helper calibration/profiling pending. Regenerate neither fixture nor tolerance using holdout results.
- [ ] Run pure selector/confirmation integration tests and one three-atom bounded smoke invocation. No timing campaign during any other stream's CPU lease.
- [ ] Commit `:sparkles: Add: SPME measured cost selection and confirmation`.

## Task 6: Add required opt-in actual phase profiling

**Files:** `src/spme.rs`, new `tests/spme_profile.rs`, example profile command, stale comment in `tests/electrostatics.rs`. **Consumes:** unchanged analytic model and Task5 confirmed IDs/baseline. **Produces:** twelve-phase compute_profiled interface and separate diagnostic records. This task is mandatory for full B05 closure; only its invocation is optional in ordinary simulations.

- [ ] Write failing parity, additive and lifecycle tests on 3-ion charged/scaled and64-ion perturbation fixtures. Require each scalar/component U/F/W difference <=1e-12*max(1,abs(value)); preload nonzero force and W values to verify additive semantics. Check required phase fields and repeated cache reuse without timing thresholds.
- [ ] Factor internal compute and reciprocal passes with a compile-time no-op/profile observer. Ordinary ForceModel::compute selects no-op and executes no Instant reads, profile allocations or synchronized profile stats. `compute_profiled` returns `(energy,SpmeProfile)`; retain all existing workspace/plan/scratch logic.
- [ ] Place twelve timings exactly: neighbor, real_space, pair_correction, workspace_setup, influence_setup, spread, fft_forward(includes transpose), influence(mode energy/W/convolution), fft_inverse(includes transpose), gather, self_background, total. Setup/recomputation and reciprocal influence sum must not be conflated. Ordinary production differentiation remains `m4_deriv`; correct only stale `ik` descriptions.
- [ ] Implement profile mode operating only on frozen confirmed IDs and the explicit baseline. Five warmups/five batches20 calls for each; cold/warm records separated. Keep this observer overhead out of ordinary TimingRecord and candidate ranking. Retain individual durations in `phase-times.jsonl` and CSV; total need not exactly equal phase sum.
- [ ] Run `rtk cargo test -p modena-core --test spme_profile`, `rtk cargo test -p modena-core --test force_cache`, `rtk cargo test -p modena-core --test force_hardening`, `rtk cargo test -p modena-core --lib spme::tests`. Review the ordinary instantiated code path for timer/observer allocation calls; do not claim zero overhead without that structural check.
- [ ] Commit `:sparkles: Add: opt-in SPME phase profiling`.

If this task is intentionally staged later, report helper-only status with phase profiling pending and stop short of full B05 closure. Existing reusable FFT plans or total times are not substitute evidence.

## Task 7: Reproducible CLI, retained reports and sequential campaign launcher

**Files:** binary, metrics report code, `bench/spme_calibration_campaign.py`, integration tests, protocol doc. **Consumes:** all implemented modes and data records. **Produces:** concrete artifact tree, exact provenance checks, CPU/memory accounting, tables and cost/error SVG.

- [ ] Implement strict CLI modes from spec: prepare/verify/sweep/confirm/profile/report; unknown/missing flags reject. prepare supports family, replicas, threads and prechosen extended reference; every family is separate. Exit codes0 success,1 I/O/internal,2 invalid/unsupported/budget,3 reference failure,4 unattained/confirmation failure.
- [ ] Invoke `shasum -a 256 ABSOLUTE_PATH` directly with std::process::Command arguments and no shell. Check status and64 hex digest. Verify source, protocol, exact snapshot and prior phase hashes at every boundary; changed physical/source provenance invalidates prior references/costs. No salted DefaultHasher or unchecked user label substitutes for file hashes.
- [ ] Use `create_new`/new directories for generated files; never overwrite production inputs or previous report evidence. Preserve failed/partial candidate rows, reference ladder/error scores, raw ordinary batches, phase records, budget status, diagnostics and actual exit code. Use exact round-trip f64 selected fragments and parse a merged minimal deck with InputConfig to prove compatibility.
- [ ] Render CSV/SVG from retained raw records, escaping all text. Plot ordinary ns/call vs observed relative RMS on log axes with target lines and explicit failed/held-out statuses. Include absolute/max/E/W tables, net force and translation ranges, cold/rebuild/warm distinction and baseline. A failed target still appears; no speedup solely from fewer grid points.
- [ ] Implement sequential launcher using Python stdlib subprocess; one child process at a time, fixed family and thread roster, explicit release Cargo invocation under `/usr/bin/time -l`, environment RAYON_NUM_THREADS matching protocol, stdout/stderr saved in separate files, numeric macOS peak-RSS parsing into `memory.json`. Child deadlines use remaining cumulative family budget; overflow/no RSS record/cap violation prevents FullScopeMeasured. Release binary build occurs before timing lease; no competing benchmark runs.

```python
families = [("cubic",2),("ortho",2),("special",2),("cubic",3),("ortho",3)]
threads = [1,4]
# The launcher executes prepare, verify, sweep, confirm, profile, report
# sequentially for every pair; paths include family/replicas/threads.
# subprocess.run receives argument lists; shell=False remains the default.
```

- [ ] CLI launcher entry: `rtk proxy python3 bench/spme_calibration_campaign.py --out /private/tmp/modena-b05-campaign-v1 --threads 1 4`. It validates a completely new output directory; creates per-family/replica/thread subdirectories; stops computation for a family on reference/physics/budget failure while retaining reports, then proceeds to other frozen independent families. A target exit4 retains the failed decision, completes unaffected target confirmation and measures the physically validated baseline plus confirmed targets in the profile phase; it does not fall back or mutate protocols. Final report receives explicit --campaign, --profiles, --memory and --out paths as in spec7.4.
- [ ] Tests: `selected_fragment_roundtrips_explicit_schema`, `profile_does_not_rewrite_primary_timing`, `source_or_snapshot_change_invalidates_phase`, `missing_hash_utility_fails_before_compute`, `output_directory_is_not_overwritten`, `partial_failure_rows_survive_report`, `rss_over_cap_blocks_full_scope`, `svg_and_csv_are_deterministic`. Test launcher with a stub subprocess executor and recorded tiny fixture output; do not run full timing in CI.
- [ ] Run focused integration tests; manually inspect generated three-atom SVG/table and compare raw values. Commit `:sparkles: Add: reproducible SPME calibration reports`.

## Task 8: Execute the frozen physical/error/cost campaign

**Files:** no algorithm/fixture changes during campaign; add only final small summary to protocol doc after raw records exist. **Consumes:** passed physical/CLI/profile suites, fixed protocol/snapshot hashes, exclusive released measurement lease. **Produces:** actual retained outcomes for each target/family/thread, actual phase records and cost/error data. This is the first full measurement task.

- [ ] Before measurement, recheck source hashes, CPU/toolchain, release profile/features, no competing jobs and exact thread count. Run required regression suites once after the final implementation change:

```sh
rtk cargo test -p modena-core --test spme_calibration
rtk cargo test -p modena-core --test spme_profile
rtk cargo test -p modena-core --test electrostatics
rtk cargo test -p modena-core --test force_hardening
rtk cargo test -p modena-core --test pair_scales
rtk cargo test -p modena-core --test force_cache
rtk cargo test -p modena-core --test final_integration
rtk cargo test -p modena-core --test sampling_state
rtk cargo test -p modena-core --lib spme::tests
rtk cargo build --release -p modena-core --example spme_convergence
```

- [ ] Freeze/save protocol and snapshots without force calls; launcher records immutable hashes before verify/evaluation. Its finite roster is cubic64,ortho64,special3,cubic216,ortho216 at1 then4 threads. Original explicit baseline per family: kappa=1/a, rc=.40*Lmin, kspace_cutoff=12/a; label its measured errors even if it misses target. Baseline need not be a selectable grid candidate and does not change winner selection.
- [ ] Run one sequential launcher command (Task7) under the explicit CPU lease. Reference physical/convergence checks must finish before any candidate timing. Each family follows prepare->verify->error sweep->ordinary timings/frozen decision->held-out confirmation->separate profile->report. No concurrency, grid growth or changed tolerances. Caps exactly as spec7.2/7.3.
- [ ] Stop family computation and retain explicit failure on reference precision/stagnation, physical failure or cap. Retain unattained/holdout target failures while completing unaffected targets and baseline phase diagnostics; do not rerun or fall back until success. A preset extended protocol, if desired later, gets new output/protocol hash and exploratory status. No universal bound or trajectory/general charge-distribution claims.
- [ ] Inspect every target/family row: exact D/R/A/M/E/W, raw net force, translation range, oracle deltas/cross-kappa, raw ordinary and profile timing separation, cold/rebuild/warm values, CPU/work/memory accounting. Reference oracle is all-image/full±k, not current direct Ewald at matching cutoffs.
- [ ] Add a concise measured outcome summary linking/hashing external artifacts. State TargetNotReached/ConfirmationFailed when observed; do not claim those targets accepted. Report full-scope tooling separately from per-target scientific acceptance and measured speedup. Scalar isotropic W only, no B07 tensor claim.
- [ ] Commit only small docs summary and artifact references: `:memo: Document: B05 SPME calibration outcomes`. No large generated measurement files, no source edit after measurements without invalidating/repeating affected evidence.

## Completion review and remaining empirical questions

- [ ] Every spec closure-matrix row has actual evidence or an explicit pending/failure state; no proposal is reported as executed.
- [ ] Both force targets, absolute/max/E/W companions, neutral/background/scaled fixtures, translations, held-out confirmation, bounded work, raw retained cost/error curves and actual phase profiling are present.
- [ ] Existing explicit parameters/defaults and analytic differentiation remain unchanged; no force-mean removal, P3M/ik certificate, full tensor or universal domain claim.
- [ ] The post-B04/B03 interfaces and peak-RSS launcher parser were checked on the actual execution base; source and fixture hashes match the campaign.

No unresolved physics-design choice blocks implementation. Remaining **empirical** questions are whether the frozen grid reaches each target on each family, whether reference/cancellation budgets are feasible in f64, whether any ordinary cost benefit survives equal-error comparison, and how real/spread/FFT/gather costs distribute on the actual host. The controller preflight base/interface recheck and host RSS parser verification are execution prerequisites. This draft establishes no target acceptance, speedup, or completed B05 implementation.

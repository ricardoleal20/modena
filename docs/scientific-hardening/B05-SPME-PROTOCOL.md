# B05 SPME empirical calibration specification

**Status:** adopted implementation contract, 2026-11-06. Controller source/interface preflight is complete; no B05 implementation, scientific tests or measurements have run. This specification supersedes the combined report for design decisions where it explicitly refines them. The companion repository plan is `docs/superpowers/plans/2026-11-06-b05-spme-calibration.md`. Original preparation drafts remain preserved separately.

**Goal:** select existing explicit SPME parameters on a bounded, frozen ionic snapshot domain using an independently converged direct Ewald reference, empirical force targets 1e-3 and 1e-4, companion error gates and measured cost; retain all evidence and report unsupported, unconverged, unattained, or unconfirmed results explicitly.

**Authority:** B05 in `docs/research/2026-11-05-core-improvements/phase3_analysis/ACCEPTANCE-BACKLOG.md`, C2 in `phase2_investigation/CASES.md`, and section 6 of `phase6_revision/REPORT.md`. Those documents propose acceptance criteria and do not assert results. Previous preparation: `/private/tmp/modena-b05-design-research.md`.

**Implementation location:** an offline `modena-core` example and example-only support modules, plus an additive opt-in SPME phase-profile API. Production explicit input schema and defaults remain unchanged. The initial verified B03 base12e4ab1 includes the benchmark library and clap/rayon/serde_json dependencies. Controller preflight confirmed seven SPME/Ewald/input/runner/pair-scale/benchmark source identities match current B04. Early Tasks1–3 add only independent example/test types, fixtures and an oracle; final B04 integration and interface checks precede Task4 and the numerical campaign. Keeping the helper in core examples avoids a new dependency from core to benchmark code and avoids edits to B03-owned files. Reuse stable data conventions after that recheck without broadening this feature into a benchmark-framework refactor.

**Execution workflow:** the controller performs the compatibility/ownership/scaffolding preflight without a separate agent or review task. Execute implementation Tasks 1–8 with `subagent-driven-development`, sequentially: one implementer per task followed by independent review before the next task. This delegation is authorized; it does not authorize work during the current B03 CPU lease. New branches use the user-required `ricardo/` prefix, overriding the host default `ricardoleal20/`.

## Global Constraints

- B03 is delivered. Controller may hand off independent example-only Tasks 1–3 on the verified B03 base while B04 proceeds in a separate checkout. Reconcile the unpublished B05 base with final reviewed B04 before Task 4 candidate evaluation or production profiling. No concurrent work within a B05 task sequence; no competing builds/tests/writes during another case's exclusive measurement lease.
- Preserve explicit `EwaldSpec` parameters, current defaults, unit conversion, runner construction, analytic order-4 B-spline differentiation, existing FFT-plan/workspace reuse, additive ForceModel semantics, and charged-mesh COM behavior.
- Never attribute a P3M/ik estimator or a universal error bound to this analytic SPME implementation. No hidden force-mean subtraction.
- Support only reduced units, finite orthorhombic cells with 3D periodicity, candidate SPME rc < half the shortest cell edge, and the exact frozen fixture families. Reject unsupported geometry/units and unvalidated inputs.
- Force targets are empirical; independent all-image real-space/full-positive-and-negative-k Ewald convergence and held-out confirmation are mandatory before emitting selected configurations.
- Absolute RMS, maximum atom-vector error, energy/atom and scalar isotropic-strain W/atom gates accompany relative RMS. Current Virial is scalar: B05 does not implement or certify the six-component stress tensor; that belongs to B07.
- Freeze fixtures, tolerance, operation/memory budgets and timing protocol before candidate measurements; preserve failed rows and raw data.
- Actual spread/FFT/gather/real-space phase measurements are required for full B05 completion. Their observer is opt-in per call, and primary candidate timing uses the ordinary uninstrumented compute path. Staging the offline-helper alternative may leave profiling explicitly pending; it cannot close full B05.
- No new Rust dependency, production autotuning, interpolation-order change, P3M replacement, tensor extension, automatic input rewrite, or unbounded refinement.
- All future shell commands use rtk; no competing timing campaign or source mutation during another stream's exclusive measurement lease.
- New branches use the user-required `ricardo/` prefix, overriding the host default `ricardoleal20/`.
- The controller owns preflight/scaffolding; Tasks 1–8 use authorized sequential subagent-driven development with one implementer per task and independent review before advancing.

## 1. Observed implementation and constraints

The implementation uses analytic derivatives of the B-spline assignment weights, despite several old comments saying `ik`. `spme.rs:641–681` differentiates the mesh energy with `m4_deriv`; this is the governing code, not the stale label. Correct those labels during the later narrowly scoped documentation change; do not change differentiation.

- `SpmeForceModel::new(kappa, real_cutoff, kspace_cutoff) -> Result<Self>` validates finite positive parameters. `with_pair_scales(Vec<PairScale>) -> Result<Self>` supplies independent Coulomb weights. `ForceModel::compute(&self, &mut System, &mut Virial, &mut NeighborList) -> Result<f64>` accumulates forces and scalar isotropic strain virial; the helper must clone/reset forces and virial before each evaluation.
- `mesh_for_box(&SimBox, kcut) -> [usize; 3]` maps `ceil(L_axis*kcut/pi)` to the next 5-smooth integer, minimum four. `kspace_cutoff` is a mesh-size selector here, not an identical spherical reciprocal truncation to direct Ewald. Deduplicate equal `(kappa, rc, actual_mesh)` candidates.
- `SpmeWorkspace` already owns forward/inverse FFT plans, two complex grids, influence factors, and per-worker scratch; the model's mutex covers mesh validation and use. Influence keys include kappa and orthorhombic edge bits; mesh changes rebuild the workspace. Standalone `spme_reciprocal_forces` and `spme_reciprocal_energy` allocate temporary workspaces and must not be used to claim production steady-state timings.
- SPME rejects triclinic boxes. C2 excludes them. Calibration accepts only 3D periodic finite positive orthorhombic cells, finite coordinates/charges, nonzero charge scale, unique valid pair indices, and no overlapping particles. Minimum particle separation must be at least `0.2*a` in the proposed reduced-unit campaign.
- Current real-space path uses minimum images and a finite cutoff. All emitted candidate cutoffs must be `< L_min/2`; the independent oracle must include omitted periodic real-space images, so real-space truncation error is actually measured. Increasing current direct Ewald cutoff beyond half the box does not create a full image oracle.
- `EwaldSpec` is explicit (`method`, `kappa`, `real_cutoff`, `kspace_cutoff`), with deny-unknown-fields; real-unit conversions are at the I/O boundary. Leave it and `runner::build_electrostatics` unchanged. Initial helper accepts reduced units only and emits reduced-unit fragments. Real units are rejected with a clear message, not guessed or converted twice.
- Both electrostatic models include self energy, conducting/tin-foil periodic reciprocal convention, uniform neutralizing background for net charge, and shared special-pair semantics. The oracle independently implements the same boundary/background convention and removes a bare minimum-image pair fraction; it does not reuse `CoulombPolicy`, `exclusion_correction`, neighbor lists, core Ewald k-vector code, or SPME.
- `conserves_total_momentum` intentionally returns false for charged SPME. Report raw net force and raw force errors. Do not subtract force means in metrics, model, fixtures, integrator, or report.

Existing regression anchors to retain: `tests/electrostatics.rs` (Ewald force FD/Madelung and ionic Ewald/SPME NVE), `tests/force_hardening.rs` (background W, full SPME charged strain, additive contract), `tests/pair_scales.rs` (weighted corrections inside/outside cutoff), `tests/force_cache.rs` (workspace reuse and shared-model cell state), `tests/final_integration.rs` and `tests/sampling_state.rs` (charged-mesh COM/DOF contract), plus the reciprocal energy/force/strain unit tests in `src/spme.rs`. Their existence is not evidence that B05 targets pass; current mesh/direct force comparison uses matching cutoffs and a much looser tolerance.

## 2. Freeze the physical protocol before measuring candidates

Create `docs/scientific-hardening/B05-SPME-PROTOCOL.md` first. It contains the exact fixture generator, candidate grid, oracle settings, tolerances below, timing procedure, and schema version `modena-spme-calibration-v1`. Hash its bytes and record that hash before any candidate error/timing run. Changes after seeing results create a new protocol version and an explicitly exploratory campaign; do not silently tune fixtures/tolerances.

Use `a=2.2`, `q_unit=1`, `F_unit=q_unit^2/a^2`, `E_unit=q_unit^2/a`, and reduced units throughout. Base ionic fixture is the existing `lattice::rocksalt(a,[2,2,2],"Na","Cl",1,1,1,-1)` (64 ions), positions wrapped explicitly before perturbation. Generate each component displacement by a fully specified integer function, not an unstable RNG stream:

```rust
// i: atom index, axis: 0..3, seed: listed below; arithmetic u64.
let z = ((i as u64 + 1) * 1_664_525
    + (axis as u64 + 1) * 1_013_904_223 + seed * 69_069) % 65_536;
let displacement = amplitude * a * (2.0 * z as f64 / 65_535.0 - 1.0);
```

Freeze candidate-selection neutral snapshots `(seed,amplitude)=(7,0.02),(19,0.05),(31,0.08)`; confirmation-only neutral snapshots `(43,0.035),(59,0.065)`. Save exact materialized coordinate/charge/box/pair-scale arrays and hashes before computing candidates; confirmation snapshots do not rank candidates. Add an unperturbed crystal for the near-zero-denominator test, never as the only accuracy fixture.

Other mandatory fixtures:

1. **Net charge/background:** selection seed 19 perturbed 64-ion snapshot, change atom 1 charge from -1 to -0.5 (`Q=+0.5`); confirmation seed 59 with the same change. All background/self terms stay enabled and Q, sum(q²), boundary convention are reported.
2. **Special pairs, rectangular cell:** `L=[9,10,11]`, positions `[[1.1,2.3,3.4],[2.7,3.2,3.8],[6.2,5.8,7.0]]`, charges `[1,-0.7,-0.3]`; pair scales `(0,1,c=0)` and `(0,2,c=0.4)`, LJ weights one. First pair is inside every proposed rc and second outside every rc. Also evaluate c=1 and c=0.7 pair-policy algebra separately. Confirmation positions add `[0.031,-0.017,0.023]` to atom 1 before wrapping. A charged variant changes final charge to +0.4 (`Q=+0.7`). Model and neighbor-list scales are supplied consistently, exactly once.
3. **Orthorhombic anisotropy:** uniformly remap fractional coordinates of the full cubic roster into the orthorhombic cell: neutral selection seeds/amplitudes `(7,.02),(19,.05),(31,.08)`, neutral confirmation `(43,.035),(59,.065)`, background selection `(19,.05)` and confirmation `(59,.065)` with atom 1 charge `-0.5`, and the unperturbed regression crystal. For replicas 2 the edges are `[4.4,5.06,5.72]`; replicas 3 uses `(a*3)*[1,1.15,1.3]`. The originally named seed 31 and seed 43 remain included with their original roles. This controller clarification was frozen before any force evaluation; it makes the section 5 “their anisotropic variants” grouping explicit. This is a separate declared fixture family, not permission to generalize to arbitrary aspect ratios.
4. **Integer periodic wrap:** shift any atom by each full box vector; E/F/W must agree to regression tolerances. This is distinct from sub-grid translation.
5. **Unsupported/degenerate fixtures:** triclinic box; empty/zero-charge input; NaN/Inf; overlapping atoms; invalid/duplicate conflicting pair policies; unsupported units; particle count or mesh over caps. Ensure explicit rejection/status with no configuration emitted.

For each candidate and every nondegenerate selection/confirmation snapshot, use seven rigid sub-grid translations. Set `h_axis=L_axis/mesh_axis`, add `h * t` componentwise to all positions and wrap, with
`t=[[0,0,0],[0.125,0.375,0.625],[0.25,0.5,0.75],[0.5,0.5,0.5],[0.75,0.25,0.5],[0.37,0.61,0.83],[0.93,0.17,0.41]]`.
The oracle reference is invariant under a rigid translation up to its convergence tolerances; compare each translated SPME result to the same converged oracle forces/energy/W after independently verifying that oracle invariance. Report force variation and E/W ranges across translations, plus raw net force for every row. Candidate-specific offsets are exact deterministic probing of mesh phase, not a claim that seven shifts bound all phases.

For force and isotropic-strain finite differences, record/check that no real-space pair, reference image, or reference reciprocal vector crosses a spherical cutoff during either perturbation. Hold kappa and real cutoff fixed. SPME W differentiates the functional at a **fixed mesh**, not a discontinuous re-selection of mesh size. Require `mesh_for_box` at both strained cells to equal the central actual mesh before comparing production full-model W with FD. If a candidate kcut lies exactly at a mesh threshold, choose and record a kcut strictly inside the interval producing the identical mesh (it defines exactly the same finite functional); otherwise mark that FD case `MeshBoundary` and use the existing public fixed-mesh reciprocal energy/force helper plus independently checked real/correction/background pieces for the derivative check. Never treat a difference across a mesh or neighbor cutoff jump as a force/strain error certificate. Report changes of mesh/caches under a finite box change as a separate lifecycle regression.

## 3. Independent direct Ewald oracle and convergence proof

Implement the reference in example-only code. It consumes flat positions/charges/orthorhombic cell and a canonical pair-scale list; it must not call an engine force evaluator or engine neighbor/min-image helper. Use independent componentwise `d -= L*(d/L).round()` for bare special-pair corrections. For the full sum, enumerate image vectors explicitly.

For `d_ij,n = r_i-r_j + n*L`, full ordered i,j,image enumeration (exclude i=j,n=0) and full positive/negative reciprocal enumeration, use:

```text
U_real = (1/2) sum' q_i*q_j*erfc(kappa*r)/r, r=|d_ij,n| < R
F_real,i = sum_(j,n)' q_i*q_j*[erfc(kappa*r)/r^3
                       +(2*kappa/sqrt(pi))*exp(-kappa²*r²)/r²]*d_ij,n
W_real = (1/2) sum' q_i*q_j*[erfc(kappa*r)/r
                              +(2*kappa/sqrt(pi))*exp(-kappa²*r²)]
k = 2*pi*(m_x/L_x,m_y/L_y,m_z/L_z), 0 < |k| < K
C=2*pi/V; A_k=exp(-k²/(4*kappa²))/k²
rho = sum_j q_j*(cos(k.r_j)+i*sin(k.r_j))
U_rec = C*sum_k A_k*|rho|²
F_rec,i = 2*C*q_i*sum_k A_k*k*(rho.re*sin(k.r_i)-rho.im*cos(k.r_i))
W_rec = C*sum_k A_k*|rho|²*(1-k²/(2*kappa²))
U_self = -kappa/sqrt(pi)*sum_i q_i²; W_self=0
U_bg = -pi*Q²/(2*V*kappa²); W_bg=3*U_bg; F_bg=0
For each canonical (i,j,c): U -= (1-c)*q_i*q_j/r_min
                         W -= (1-c)*q_i*q_j/r_min
                         F_i -= (1-c)*q_i*q_j*d_min/r_min³; F_j += same
```

Use compensated accumulation for energies/virials and force components. `W=-dU(lambda*r,lambda*L)/d lambda at lambda=1` with kappa fixed (not a full tensor). Self-image real energy/virial is included, and self-image forces cancel by image symmetry. For each axis, real-image enumeration bounds `|n_axis| <= ceil((R + max_ij |r_i,axis-r_j,axis|)/L_axis)` safely before spherical filtering; wrapped coordinates allow `ceil(R/L_axis)+1`. Reciprocal bounds are `ceil(K*L_axis/(2*pi))`. Record counts, caps, and truncation settings.

Default reference ladder fixes `kappa_ref=1.1/a` and evaluates `(R/a,K*a)=(4,8),(5,10),(6,12),(7,14)`. Require **two consecutive** refinement differences satisfying all reference budgets; otherwise emit `ReferenceUnconverged`, not a selected config. Evaluate final settings independently at `kappa_ref=0.9/a` and `1.3/a` and require agreement with the final 1.1/a result to the same budgets; if not, no calibration. These are empirical convergence checks, not a rigorous infinite-series error bound.

Reference budgets on each snapshot: force vector RMS delta <= `min(1e-8*F_unit,0.01*1e-4*D)` where D is the final reference force vector RMS when D > `1e-8*F_unit`; for a smaller D use absolute budget `1e-8*F_unit`. Max atom-vector force delta <= `1e-7*F_unit`. Energy/atom and W/atom delta <= `1e-9*E_unit`. A failed ladder is reported with its raw deltas and cannot be repaired by increasing thresholds after candidate inspection. Explicit later `--reference-extended` protocol may allow two additional predetermined levels `(8,16),(9,18)`, retaining the same budgets and operation caps; it is a new recorded protocol option, not unbounded iteration.

Validate this oracle before using it: its full analytic forces match central energy differences on the 3-atom rectangular neutral/charged/scaled fixtures; analytic W matches isotropic strain FD; full ±k versus an independently coded half-space sum agrees; translation invariance; known background coefficient; special-pair bare-fraction identities. Cross-check current `EwaldForceModel` at high kappa and sub-half-box rc against the oracle only where the independent ladder shows omitted real-image/tail effects under the same budget. Current direct Ewald is a useful secondary implementation cross-check, not the independent oracle itself. Keep a Madelung sequence check (increasing R/K and finite-cell lattice convergence) to catch sign/factor mistakes; do not call the existing loose 0.02 Madelung tolerance a 1e-4 force certification.

## 4. Error metrics, gates, and domain language

Per snapshot s with N atoms, define raw `delta_i=F_spme,i-F_ref,i`:

```text
D_s = sqrt(sum_i |F_ref,i|²/N)                 # vector RMS, NOT component RMS
A_s = sqrt(sum_i |delta_i|²/N)                 # same vector RMS convention
R_s = A_s/D_s                                 # only when D_s > 1e-8 F_unit
M_s = max_i |delta_i|                          # maximum atom-vector norm
E_s = (U_spme-U_ref)/N                         # signed + absolute reported
W_s = (W_spme-W_ref)/N                         # signed + absolute reported
net_force = sum_i F_spme,i                     # three components and norm
```

Never normalize by candidate forces, maximum reference force, signed mean force, energy zero, or a reference force averaged across distinct fixtures. Preserve per-snapshot denominators and worst-case reductions. Report absolute/component maxima too, with names distinct from atom-vector maximum.

For target `t in [1e-3,1e-4]`, a candidate passes only if **every selection snapshot and translation** satisfies `R_s <= t`, `A_s <= t*F_unit`, `M_s <= 10*t*F_unit`, `|E_s| <= t*E_unit`, and `|W_s| <= 10*t*E_unit`, with oracle converged. These numerical side gates are proposed protocol tolerances, frozen before measurement; they are not already demonstrated attainable. Near-zero D snapshots report `RelativeUndefined` and are tested against absolute/energy/strain regression budgets; they do not rank or certify relative-target calibration. Relative-undefined arbitrary input cannot receive a relative-target configuration.

After selection, the exact selected candidate must pass the same gates on every held-out confirmation snapshot/translation. Failure produces `ConfirmationFailed`; the helper does not silently select the next candidate based on holdout results. A new selection/confirmation campaign requires a new manifest. If no feasible grid candidate exists, emit `TargetNotReached` with the full Pareto curve and no config. Failed or capped oracle -> `ReferenceUnconverged`/`BudgetExceeded`, also no config.

Calibration wording: **"Observed target met on these exact snapshots, shifts, and boundary convention."** It is not a universal error bound. A report records charge distribution, N, density, minimum separation, aspect ratio, Q, sum(q²), pair scales, box/position/charge hashes, seed/amplitude, method `spme-order4-analytic`, and raw reference convergence data. Only exact fixture-hash matches are `CalibratedSnapshots`. Another snapshot in a similarly summarized domain is `NotCalibrated`, because scalar summaries cannot certify arbitrary charge arrangements. Changed box, charges, pair scales, interpolation algorithm, units, or cutoffs invalidate the report. Dynamic trajectories and NPT box changes require separate sampling and validation; no runtime guarantee is made by this helper. Partial PBC and triclinic input are rejected, not labeled approximately calibrated.

## 5. Bounded sweep, cost selection, and output

Initial candidate grid, per fixture/cell family:

```text
kappa*a in [1.2,1.6,2.0,2.4,3.0,3.6]
rc/L_min in [0.30,0.40,0.49]
kspace_cutoff*a in [8,12,16,24,32,48,64]
```

126 raw candidates, sorted lexicographically; deduplicate actual mesh triples as above. No expected-error formula prunes or certifies candidates. Hard limits: 256 atoms/reference input, 1,048,576 mesh points/candidate, 126 candidate tuples/family, 25,000 candidate evaluations/campaign, 4 reference refinement levels (+2 only explicit extended protocol), 100,000,000 real ordered-pair-image evaluations and 100,000,000 atom-k terms per individual reference call. Count estimated work safely before allocation; use checked mesh products. If a whole campaign exceeds a declared total evaluation cap, reject before partial selection. All caps, elapsed time, and skipped reasons appear in JSON. No automatic grid expansion.

Family grouping is explicit: the 64-ion cubic neutral/background snapshots share one triple and must all pass it; their anisotropic variants form a second cell family; the rectangular 3-ion neutral/charged/special-pair fixtures form a third family. A winner is calibrated only for its own family; the three-atom special-pair regression is additionally mandatory before any ionic-family campaign, but its separate rc grid does not secretly alter the cubic winner. Repeat the ionic families independently for 216 ions. The 25,000 cap is per family campaign, including translation and confirmation calls, and excludes separately counted bounded oracle work and timing repetitions.

Pipeline: materialize/hash protocol and fixtures -> oracle and physical checks -> reference convergence -> all selection force/error rows -> feasible set per target -> timed cost for feasible candidates -> fixed winning candidate -> held-out confirmation -> output. Unsuccessful physical/error gates precede timing; timing does not rescue a physics failure. For uninstrumented total-cost ranking, measure every feasible candidate (within 126 cap), using its persistent model and neighbor list. Do not use an analytic mesh-size proxy as measured cost. If bounded timing budget cannot cover all feasible candidates, report shortlist/proxy status and do not claim fastest feasible.

The explicit frozen baseline is `kappa=1/a`, `real_cutoff=0.40*L_min`, `kspace_cutoff=12/a` for each family, evaluated and retained even when its observed error misses a target. It is not silently inserted into the selectable 126-tuple grid. Profile its actual phases after its derivative/additive/strain physics checks, alongside any confirmed target candidates; baseline target failure is reported rather than hidden.

For each feasible candidate, use the 64-ion selection fixture and a separately frozen 216-ion (`replicas=[3,3,3]`) family only after its own oracle/selection protocol completes. Report separate winners per N/family; do not combine to manufacture a universal winner. The 216-ion family uses the same perturbation seeds and amplitudes, background rule, and translations; its explicit confirmation fixtures are equally required. Larger N can be exploratory profiling only and is uncalibrated unless a new reference-supported campaign is registered.

Timing: release build; 1 thread primary and a separately labeled 4-thread campaign; fixed CPU/toolchain/config; no concurrent B03 benchmark. Keep model/list alive across calls; reset force and W each call, no cloning in timed region. Five warmups, then five batches of 20 force evaluations with `black_box` of U/W/forces; emit every batch, median ns/evaluation, IQR, observed min/max. List build cache check remains in timed compute. Report cold construction/first-call timing separately and one forced-neighbor-rebuild evaluation separately, never mixing these with warm timings. Choose minimum median among candidates passing all gates; declare candidates within 5% of fastest median a tie, select lower actual grid-point count, then lower cutoff, then lexicographic kappa/kcut for deterministic output. This is a noise-aware practical rule, not a statistical significance claim.

Output files in an explicit user-supplied directory: `protocol.json`, exact `snapshots.json`, `reference.json`, `candidates.jsonl` (all accepted/failed/capped rows), `summary.json`, `cost-error.csv`, `cost-error.svg` (log force-error vs ns/evaluation, target guides, failures and confirmation status labeled), and `selected-1e-3.toml`/`selected-1e-4.toml` only after confirmation succeeds. For multiple families, place each in its own N/cell-named subdirectory. SVG can be generated with escaped text and deterministic simple polyline/scatter markup; no plotting dependency needed. Include energy/W/max-error tables, not only the force RMS chart. Selected fragment is exactly:

```toml
# Example of valid output syntax; these numbers are illustrative, not selected.
# Real output comments name the report path and protocol/snapshot/source hashes.
[electrostatics]
method = "spme"
kappa = 1.0
real_cutoff = 2.0
kspace_cutoff = 12.0
```

The implementation writes actual selected finite numbers with round-trip f64 precision. The user merges the fragment into an existing `units="reduced"` deck. The helper never overwrites a production input, restart, weights, or prior report. Default existing output directory is an error unless explicitly writing a new unique run subdirectory.

## 6. Files and interfaces for later implementation

Create:

- `docs/scientific-hardening/B05-SPME-PROTOCOL.md`: frozen science/timing protocol, status language, usage, rejection rules, URLs, and manifest schema.
- `crates/modena-core/examples/spme_convergence.rs`: thin example entry point, manual std::env argument parsing, file handling, JSON/TOML output and process exit status.
- `crates/modena-core/examples/spme_calibration/mod.rs`: serializable types, input validation, bounded sweep orchestration, deterministic selection; expose its submodules to the integration test by path.
- `crates/modena-core/examples/spme_calibration/reference.rs`: independent all-image, full-±k Ewald and convergence bookkeeping, no core evaluator calls.
- `crates/modena-core/examples/spme_calibration/fixtures.rs`: deterministic fixture materialization and exact selection/confirmation splits.
- `crates/modena-core/examples/spme_calibration/metrics.rs`: pure error metrics, budget checks, selector and SVG/CSV rendering.
- `crates/modena-core/tests/spme_calibration.rs`: integration tests importing the example support module with `#[path="../examples/spme_calibration/mod.rs"] mod calibration;`. No public core API for the offline selector.
- `crates/modena-core/tests/spme_profile.rs`: additive/physics parity and phase/cache regression tests for the opt-in profiling interface.
- `bench/spme_calibration_campaign.py`: sequential standard-library host launcher for immutable family/thread campaigns, stderr/stdout retention and `/usr/bin/time -l` peak-RSS records. No B03-owned source is modified.

Modify later, separately: `crates/modena-core/src/spme.rs` only for the additive phase-profile interface and corrected analytic-force comments; `crates/modena-core/tests/electrostatics.rs` only to correct its stale `ik` explanatory comment. No planned edits to io/input.rs, runner.rs, core Ewald, `Cargo.toml`, or B03 benchmark files.

Concrete example support interfaces (all `pub(crate)` within example/test compilation):

```rust
#[derive(Clone, serde::Serialize, serde::Deserialize)]
pub struct Snapshot {
    pub id: String,
    pub cell: [f64;3],
    pub positions: Vec<[f64;3]>,
    pub charges: Vec<f64>,
    pub scales: Vec<Scale>,
}
#[derive(Clone, serde::Serialize, serde::Deserialize)]
pub struct Scale { pub i: usize, pub j: usize, pub coulomb: f64 }
#[derive(Clone, Copy, serde::Serialize, serde::Deserialize)]
pub struct Parameters { pub kappa: f64, pub real_cutoff: f64, pub kspace_cutoff: f64 }
#[derive(Clone, serde::Serialize, serde::Deserialize)]
pub struct Evaluation { pub energy: f64, pub strain_w: f64, pub forces: Vec<[f64;3]> }
#[derive(Clone, Copy, serde::Serialize, serde::Deserialize)]
pub struct ReferenceSettings { pub kappa: f64, pub real_cutoff: f64, pub kspace_cutoff: f64 }
pub fn reference(snapshot: &Snapshot, settings: ReferenceSettings,
                 limits: &Limits) -> Result<Evaluation, CalibrationFailure>;
pub fn converged_reference(snapshot: &Snapshot, protocol: &Protocol)
    -> Result<ReferenceRecord, CalibrationFailure>;
pub fn evaluate_spme(snapshot: &Snapshot, p: Parameters) -> Result<Evaluation, CalibrationFailure>;
pub fn metrics(reference: &Evaluation, candidate: &Evaluation,
               n: usize, force_floor: f64) -> Result<ErrorMetrics, CalibrationFailure>;
pub fn campaign(protocol: &Protocol, snapshots: &[Snapshot])
    -> Result<CampaignReport, CalibrationFailure>;
pub fn select(rows: &[CandidateReport], target: f64)
    -> Result<Parameters, CalibrationFailure>;
```

Define remaining named structs as exact named-field serialized records in `mod.rs`: `Limits` carries caps above; `Protocol` carries `schema`, `a`, `q_unit`, `targets`, candidate arrays, reference ladder, limits, fixture-set ID, translation array, tolerance fields, thread count, warmup/batch counts; `ReferenceRecord` carries settings/evaluations/deltas/work counts and convergence flags; `ErrorMetrics` carries `force_rms_abs`, optional `force_rms_rel`, `reference_force_rms`, `max_atom_force_error`, `max_component_force_error`, signed/absolute E_per_atom and W_per_atom, `net_force`; `CandidateReport` carries parameters/actual mesh/per-fixture metrics/status/raw batch ns and summary; `CampaignReport` carries protocol/snapshot/provenance records, all candidate reports, per-target decisions/confirmation results; `CalibrationFailure` is an enum with `InvalidInput(String)`, `UnsupportedDomain(String)`, `BudgetExceeded(String)`, `ReferenceUnconverged(String)`, `RelativeUndefined(String)`, `TargetNotReached(f64)`, `ConfirmationFailed(f64)`, `Io(String)`. No scientific error is serialized as a successful zero result. Validation returns explicit context and fixture ID.

The binary modes:

```text
spme_convergence prepare --out DIR --family cubic|ortho|special --replicas 2|3 --threads 1|4 [--reference-extended]
spme_convergence verify --protocol DIR/protocol.json --snapshots DIR/snapshots.json
spme_convergence sweep --protocol DIR/protocol.json --snapshots DIR/snapshots.json --out NEW_DIR
spme_convergence confirm --campaign DIR/summary.json --out NEW_DIR
spme_convergence report --campaign DIR/summary.json --out NEW_DIR
```

`prepare` performs no force evaluation and emits byte-stable JSON excluding timestamp from hashable protocol/fixture content. `verify` performs oracle/regression gates and records results in `verification.json` beside the frozen protocol (it does not rewrite protocol/snapshots). `sweep` requires successful verify records with matching hashes, recalculates references if provenance differs, and writes no selected configuration before confirmation. `confirm` reads frozen selected candidate(s), evaluates held-out fixtures, and emits configuration only on pass. `report` re-renders immutable result records, does not rerun/tune. Exit codes: 0 successful requested phase, 2 invalid/unsupported input or cap, 3 reference failure, 4 target not reached/confirmation failed, 1 I/O/internal error. The user may inspect a failure report even when command exits nonzero.

For deterministic SHA-256 without a Rust dependency, the example uses `std::process::Command` (no shell interpolation) to invoke the existing host `shasum -a 256 ABSOLUTE_PATH` utility on explicit paths. Validate exactly 64 hex digits and the expected single-file result, check process status, and fail before campaign work when unavailable. Source and materialized protocol/snapshot bytes are hashed in prepare and rechecked by verify/sweep/confirm, not trusted from user labels. The supported first campaign host is the current macOS environment with shasum (RTK is agent-shell tooling only); portability to a different hashing utility is an explicit later adapter. The protocol calls this host prerequisite out; no new Rust package or hidden network call is introduced.


## 7. Normative execution details and closure

This section resolves the combined research report's implementation details. Its precise API, budget, status, and profiling requirements govern the implementation.

### 7.1 Required and optional work

Required for the complete B05 deliverable: protocol/fixtures, independent oracle validation and convergence, error gates, bounded candidate selection, ordinary-path measured costs, fixed held-out confirmation, retained candidate/error/cost/provenance records, and **actual opt-in phase profiling** of real-space, spread, forward FFT, reciprocal influence, inverse FFT, gather, and setup/correction/background overhead. Phase profiling is optional for each normal force call, so it introduces no clock reads into the ordinary path. Implementing it and running its separate diagnostic campaign are required before reporting full B05 scope delivered. A helper-only intermediate result must say `OfflineHelperComplete_PhaseProfilingPending` and cannot represent profiling as inferred from total time or existing cached FFT plans. Full scope delivered can still have `TargetNotReached` or `ConfirmationFailed` at a target: those are measured scientific outcomes, not tool implementation failures or positive acceptance of the accuracy target.

Optional extensions excluded from this iteration: production input RMS target field, runtime retuning, six-component stress, arbitrary snapshots accepted as a calibrated domain, larger N, other interpolation orders, real-unit support, partial PBC, triclinic SPME, alternate estimators. A separate immutable extended-reference protocol with the two preset extra levels is allowed; it must be chosen at prepare time, not enabled after a candidate outcome.

### 7.2 Exact serializable records

All records below live in example-only `spme_calibration/mod.rs`, derive `Clone, Debug, Serialize, Deserialize`, and use `#[serde(deny_unknown_fields)]` on structs. Float fields must be finite unless explicitly represented by `Option`; reject nonfinite values rather than emit JSON null or infinity. Use the Snapshot/Scale/Parameters/Evaluation/ReferenceSettings definitions in section 6, with the following exact extensions: `Snapshot` also has `role: FixtureRole`, `family: String`, `units: String`, `pbc: [bool;3]`, `seed: Option<u64>`, `amplitude: Option<f64>`; `units` is always `reduced` and PBC always all true for prepared fixtures. Roles are `Selection`, `Confirmation`, and `RegressionOnly`. The binary may read only the fixed prepared fixture schema; custom input is rejected as uncalibrated in v1. Each prepare command writes exactly one family protocol and that family's fixtures. `--family special` requires replicas=2 and prepares the fixed three-atom fixture, not 64 ions. For ortho fixtures use edges `(a*replicas)*[1,1.15,1.3]`, i.e. `[4.4,5.06,5.72]` for replicas=2 and `[6.6,7.59,8.58]` for replicas=3; remap the cubic fractional coordinates. No candidate from another family is inherited.

```rust
pub enum FixtureRole { Selection, Confirmation, RegressionOnly }
pub struct Limits {
    pub max_atoms: usize,                    // 256
    pub max_mesh_points: usize,              // 1_048_576
    pub max_candidates: usize,               // 126 per family
    pub max_accuracy_calls: u64,             // 25_000 per family, selection+confirmation
    pub max_timing_calls: u64,               // 20_000 per family, includes cold/rebuild/warmup
    pub max_profile_calls: u64,              // 1_000 per family
    pub max_real_terms_per_call: u64,        // 100_000_000
    pub max_atom_k_terms_per_call: u64,      // 100_000_000
    pub max_real_terms_per_family: u64,      // 5_000_000_000
    pub max_atom_k_terms_per_family: u64,    // 500_000_000
    pub max_owned_heap_bytes: u64,           // 268_435_456 (256 MiB)
    pub max_peak_rss_bytes: u64,             // 536_870_912 (512 MiB), acceptance cap
    pub max_family_wall_seconds: u64,        // 1_800, cumulative verify/sweep/confirm/profile
    pub max_threads: usize,                  // 4
}
pub struct Tolerances {
    pub reference_force_rms_abs_units: f64,  // 1e-8 F_unit
    pub reference_force_target_fraction: f64,// 0.01 times tightest target times D
    pub reference_force_max_units: f64,      // 1e-7 F_unit
    pub reference_energy_per_atom_units: f64,// 1e-9 E_unit
    pub reference_w_per_atom_units: f64,     // 1e-9 E_unit
    pub relative_denominator_floor: f64,     // 1e-8 F_unit
    pub candidate_force_rms_abs_factor: f64, // 1 times target F_unit
    pub candidate_force_max_factor: f64,     // 10 times target F_unit
    pub candidate_energy_factor: f64,        // 1 times target E_unit per atom
    pub candidate_w_factor: f64,             // 10 times target E_unit per atom
    pub cost_tie_fraction: f64,              // 0.05
}
pub struct Protocol {
    pub schema: String,                     // modena-spme-calibration-v1
    pub fixture_set: String,                // c2-spme-v1
    pub family: String,                     // cubic-64, ortho-64, special-3, cubic-216, ortho-216
    pub replicas: u32,                      // 2 or 3; special-3 prepared alongside replicas 2
    pub a: f64, pub q_unit: f64,             // 2.2, 1.0
    pub targets: Vec<f64>,                  // [1e-3,1e-4]
    pub kappa_a: Vec<f64>,                  // [1.2,1.6,2.0,2.4,3.0,3.6]
    pub cutoff_lmin: Vec<f64>,              // [.30,.40,.49]
    pub kcut_a: Vec<f64>,                   // [8,12,16,24,32,48,64]
    pub reference_ladder: Vec<[f64;2]>,      // [R/a,K*a], four or prechosen six levels
    pub reference_kappa_a: f64,             // 1.1
    pub crosscheck_kappa_a: [f64;2],         // [.9,1.3]
    pub translations: Vec<[f64;3]>,         // exact seven rows in section 2
    pub limits: Limits, pub tolerances: Tolerances,
    pub threads: usize,                    // 1 or 4
    pub warmups: usize, pub batches: usize, pub batch_calls: usize, // 5,5,20
}
pub struct WorkCounts { pub real_terms: u64, pub atom_k_terms: u64 }
pub struct ErrorMetrics {
    pub force_rms_abs: f64, pub force_rms_rel: Option<f64>,
    pub reference_force_rms: f64, pub max_atom_force_error: f64,
    pub max_component_force_error: f64,
    pub energy_per_atom_signed: f64, pub energy_per_atom_abs: f64,
    pub w_per_atom_signed: f64, pub w_per_atom_abs: f64,
    pub net_force: [f64;3], pub net_force_norm: f64,
}
pub struct ReferenceLevel {
    pub settings: ReferenceSettings, pub evaluation: Evaluation,
    pub delta_from_previous: Option<ErrorMetrics>, pub work: WorkCounts,
}
pub struct ReferenceRecord {
    pub snapshot_id: String, pub snapshot_sha256: String,
    pub levels: Vec<ReferenceLevel>, pub crosschecks: Vec<ReferenceLevel>,
    pub chosen_level: Option<usize>, pub status: Status,
}
pub struct MetricRow {
    pub snapshot_id: String, pub shift_index: usize,
    pub metrics: ErrorMetrics, pub passing_targets: Vec<f64>, pub status: Status,
}
pub struct TimingRecord {
    pub reset_included: bool,               // true for ordinary primary batches
    pub cold_ns: u64, pub rebuild_ns: u64,
    pub batch_total_ns: Vec<u64>, pub calls_per_batch: usize,
    pub median_ns_per_call: f64, pub q1_ns_per_call: f64, pub q3_ns_per_call: f64,
}
pub struct CandidateReport {
    pub id: usize, pub family: String, pub parameters: Parameters,
    pub actual_mesh: [usize;3], pub aliases: Vec<f64>,
    pub rows: Vec<MetricRow>, pub timing: Option<TimingRecord>, pub status: Status,
}
pub struct TargetDecision {
    pub target: f64, pub selected_id: Option<usize>, pub selection_status: Status,
    pub confirmation_rows: Vec<MetricRow>, pub confirmation_status: Status,
}
pub struct PhaseRecord {
    pub candidate_id: usize, pub snapshot_id: String, pub threads: usize,
    pub cold: bool, pub calls: usize, pub phases_ns: Vec<[u64;12]>,
}
pub struct Provenance {
    pub git_head: String, pub source_sha256: std::collections::BTreeMap<String,String>,
    pub protocol_sha256: String, pub snapshot_sha256: String,
    pub reference_algorithm: String,        // all-image-full-pm-k-v1
    pub rustc_verbose: String, pub cargo_version: String,
    pub os: String, pub cpu: String, pub threads: usize,
    pub profile: String, pub features: Vec<String>, pub command: Vec<String>,
    pub utc_started: String, pub family_elapsed_ns: u64,
    pub measured_peak_rss_bytes: Option<u64>,
}
pub struct CampaignReport {
    pub protocol: Protocol, pub provenance: Provenance,
    pub snapshot_ids: Vec<String>, pub references: Vec<ReferenceRecord>,
    pub candidates: Vec<CandidateReport>, pub decisions: Vec<TargetDecision>,
    pub phase_records: Vec<PhaseRecord>, pub status: Status,
}
pub enum Status {
    Prepared, Verified, Passed, Failed, Pending, CalibratedSnapshots,
    NotCalibrated, RelativeUndefined, ReferenceUnconverged, ReferenceStagnated,
    ReferencePrecisionLimited, BudgetExceeded, UnsupportedDomain,
    TargetNotReached, ConfirmationFailed, MeshBoundary,
    OfflineHelperCompletePhaseProfilingPending, FullScopeMeasured,
}
pub enum CalibrationFailure {
    InvalidInput(String), UnsupportedDomain(String), BudgetExceeded(String),
    ReferenceUnconverged(String), ReferenceStagnated(String),
    ReferencePrecisionLimited(String), RelativeUndefined(String),
    TargetNotReached(f64), ConfirmationFailed(f64), Io(String),
}
```

For section 6 function signatures, use `CalibrationFailure` as the error type consistently in every helper (`reference`, `evaluate_spme`, `metrics` included), replacing the earlier generic String error shorthand. `select` considers only rows whose physical gates and ordinary timing are complete and finite, and returns a Parameters value from the same candidate ID retained in TargetDecision. `campaign` orchestrates only one family/protocol; never pool families for force denominators or cost ranking.

The Madelung sanity fixture uses the existing repository constant `alpha=1.747558775`, r0=a/2, and verifies `abs(E_per_pair + alpha/r0) <= 2e-5*E_unit` after the independent oracle ladder passes. This approximate historical constant is only a sign/normalization sanity anchor, not a 1e-4 force calibration or the reference-convergence budget. Full ±k and half-space values must additionally agree to 1e-10 scaled units, integer-wrap and rigid-translation oracle E/F/W to 1e-9 scaled units, and stable-image core Ewald cross-checks to the already stricter reference budgets. An unperturbed crystal never substitutes for perturbed force fixtures.

The regular timed evaluator is a reusable example-only `SpmeEvaluator { model: SpmeForceModel, system: System, neighbors: NeighborList }` with `new(snapshot:&Snapshot, parameters:Parameters)->Result<Self,CalibrationFailure>`, `evaluate(&mut self)->Result<Evaluation,CalibrationFailure>` for untimed accuracy, and `ordinary_batch(&mut self, calls:usize)->Result<u64,CalibrationFailure>` for timed normal compute. `ordinary_batch` includes force zeroing and `Virial::default()` once per call, calls `ForceModel::compute`, black-boxes energy/W/force contents, and returns elapsed nanoseconds for the batch; it does not clone coordinates or forces in the timed region. Report this reset policy explicitly. `evaluate` may copy returned forces after compute; that copy is not primary timing.

### 7.3 CPU, memory, precision and failure rules

Validate/estimate every image and k loop's work with checked arithmetic before allocating or executing. Stream reference image/k contributions; do not store an O(N² images) table. Check deadlines after at most 1,024 reference terms, between accuracy calls, and between timing batches; timing must not call budget instrumentation inside the measured region. The fixed operation counters are the primary CPU work budget. Cumulative family wall limit is 1,800 seconds over verify/sweep/confirm/profile with at most four Rayon workers: a bound of 7,200 worker-seconds plus bounded serial orchestration. Store elapsed work in phase records to prevent resetting budgets by splitting commands. On deadline, count, or allocation failure, emit BudgetExceeded and retain partial raw records with no selected configuration newly emitted. Memory budget covers simultaneous owned vector capacities, checked before reserve/resize; report opaque FFT-plan allocation separately. Record actual process peak RSS with the host `/usr/bin/time -l` wrapper; exceeding 512 MiB invalidates the memory budget in the final report and blocks FullScopeMeasured. A requested cap cannot be silently increased to make a target pass.

Reference refinement always runs its finite preset ladder. Define normalized convergence score S for a delta as max(force RMS/budget, max atom error/budget, energy/atom/budget, W/atom/budget). Two consecutive finite scores <=1 plus final cross-kappa checks are required. If final scores stay >1 and the last two refinements each improve by less than 10% (`S_new >= 0.9*S_old`), use ReferenceStagnated; otherwise exhausted ladder without two passes is ReferenceUnconverged. Bit-identical or tiny deltas below budgets are not alone called stagnation, because they may be convergence; cross-kappa agreement and the validated derivative/oracle checks remain mandatory. Nonfinite intermediate values or analytic/finite-difference cancellation that cannot meet the registered residuals after h and h/2 checks produce ReferencePrecisionLimited, never a relaxed reference threshold. Near-zero reference force denominators produce RelativeUndefined, not numerical target attainment.

When a tool phase fails, write its exact failure and exit nonzero; never reinterpret partial data as a measured fastest feasible winner. After a confirmation fails, do not choose another candidate using holdout data. The frozen candidate remains recorded as failed. Candidate grid, oracle ladder or tolerance changes require a new protocol hash and separate exploratory/confirmatory label.

### 7.4 Phase-profile API and artifact contract

Add in `spme.rs`:

```rust
#[derive(Clone, Debug, Default)]
pub struct SpmeProfile {
    pub neighbor: std::time::Duration,
    pub real_space: std::time::Duration,
    pub pair_correction: std::time::Duration,
    pub workspace_setup: std::time::Duration,
    pub influence_setup: std::time::Duration,
    pub spread: std::time::Duration,
    pub fft_forward: std::time::Duration,
    pub influence: std::time::Duration,
    pub fft_inverse: std::time::Duration,
    pub gather: std::time::Duration,
    pub self_background: std::time::Duration,
    pub total: std::time::Duration,
}
// In impl SpmeForceModel; core Result alias and existing types.
pub fn compute_profiled(&self, system: &mut System, virial: &mut Virial,
                        neighbors: &mut NeighborList) -> Result<(f64, SpmeProfile)>;
```

This twelve-field layout supersedes the earlier eleven-field sketch: `influence_setup` separately measures cached-key check/recomputation inside reciprocal setup; `workspace_setup` covers mutex acquisition, mesh validation and workspace construction/reuse. Forward/inverse FFT fields each include their own transposes. `influence` is the reciprocal mode energy/W/convolution pass. Total surrounds complete profiled compute and is not claimed to be exactly the sum of fields because timer/observer overhead exists. Ordinary ForceModel compute instantiates a compile-time no-op observer; no clocks, profile vector, allocation or shared statistics are introduced in that path. Twelve raw field durations are retained in the stated struct order in PhaseRecord. Profile and normal U/F/W must agree to `1e-12*max(1,|value|)` per component and scalar on the three-atom and 64-ion fixtures.

Add command `spme_convergence profile --campaign CONFIRMED_DIR/summary.json --out NEW_DIR`. It accepts only confirmed candidate IDs and the declared baseline, runs each on the same frozen fixtures with five warmups and five batches of 20 **profiled** calls, and records cold first-call and warmed phase distributions separately. It never rewrites primary TimingRecord values or reselects candidates. Include this diagnostic overhead in labels. `report` merges hash-matched confirmation and phase records and sets FullScopeMeasured only if all required records and budgets are complete; positive target status remains separately Passed/TargetNotReached/ConfirmationFailed.

The sequential launcher writes immutable `memory.json` for every child, then calls report with explicit inputs: `report --campaign CONFIRM_DIR/summary.json --profiles PROFILE_DIR/phase-times.jsonl --memory FAMILY_DIR/memory.json --out REPORT_DIR`. Profile and report accept the retained decision manifest even when one or both targets fail; profile then measures the physically validated baseline and any independently confirmed target, while failed target IDs remain failed. Nonzero exit4 does not erase reports or skip baseline phase diagnostics. Reference/physics/budget failures stop computation for that family and render the retained failure report without requiring profile data. Missing or over-cap process RSS blocks FullScopeMeasured. The prepared protocol/hash and typed decision IDs are always rechecked before merging records.

Artifact extension: `phase-times.jsonl`, `phase-times.csv` and `memory.json` are mandatory for full scope. Primary `candidates.jsonl`, ordinary raw batch totals, all error rows (including failures), reference ladder deltas and convergence scores, frozen holdout decisions, environment, source hashes, and timing/profile separation must remain independently inspectable. Selected fragments may be emitted after successful confirmation as the helper alternative; their report remains marked profiling pending until the separate profile phase completes. Do not claim measured speedup before matched-accuracy ordinary path results exist.

### 7.5 Closure matrix

| Requirement | Evidence needed | Failure leaves |
|---|---|---|
| Independent reference | all-image/full ±k validated equations, ladder deltas, cross-kappa agreement, no shared force code | uncalibrated; no config |
| Target 1e-3 / 1e-4 | per-snapshot/translation relative+absolute+maximum+E+W gates and held-out confirmation | target explicitly unaccepted |
| Physics invariants | analytic F=-grad U, fixed-mesh isotropic W, background, pair policy, translation/wrap, additive and cache regressions | implementation incomplete |
| Cost/error selection | retained all-candidate errors and ordinary primary cold/rebuild/warm timing at matched errors | no fastest-feasible or speedup claim |
| Phase profiling | actual twelve-phase records, same U/F/W, observer absent from ordinary timing | full B05 profiling pending |
| Bounded/reproducible execution | immutable hashes, elapsed/work/memory counts, process peak RSS, command/toolchain/thread/hardware details | budget/provenance not accepted |

## 8. Provenance and sources checked during this preparation

Read checkout HEAD was `aec5bfb4a3066e1df217052d19af6a13472ed6de`; working files are live and may differ from HEAD because another stream exists. SHA-256 of files read:

```text
121ed7db2606935411f7bcbb9229bd5063cdc7c60d5ffd448b1b547777e1c2a8 crates/modena-core/src/spme.rs
11d294c94eb5cefa77ad16f586f564cf3f26b14d835b15a6354971028c90bf2b crates/modena-core/src/ewald.rs
e55b9be9d3d07af674ec7928bbe4ecc192324a5e99549836f377d2e3f1a62af8 crates/modena-core/src/io/input.rs
1e46b318945079878cdcab7b344faabb72c194692c0acbb41b84afa8335d3540 crates/modena-core/src/runner.rs
c3e73916b8ef4c171515957c2944c9d7655ee35c3e742116822790acceb278eb crates/modena-core/tests/electrostatics.rs
842d365e0a216b0486c49fc976cb7ff3e353b957d549ff0b44f3cf026c17ace0 docs/research/2026-11-05-core-improvements/phase2_investigation/CASES.md
91eea5ee745963bc01bedf56c08726b9eaf6517b189277003a6b716e58965442 docs/research/2026-11-05-core-improvements/phase3_analysis/ACCEPTANCE-BACKLOG.md
```

Primary sources opened 2026-11-06:

- Essmann et al., *A smooth particle mesh Ewald method* (1995), DOI https://doi.org/10.1063/1.470117 ; full paper https://biomolmd.org/mw/images/e/e0/Spme.pdf . Its analytic gradients and Coulombic lattice/exclusion conventions support defining the oracle and preserving the energy-gradient formulation; it does not certify this checkout's target errors.
- Deserno & Holm, *How to mesh up Ewald sums I* (1998), https://arxiv.org/abs/cond-mat/9807099 and https://doi.org/10.1063/1.477414 . Its comparison of assignment/influence/differentiation supports identifying actual method ingredients rather than relying on a PME label.
- Deserno & Holm, *How to mesh up Ewald sums II* (1998), https://arxiv.org/abs/cond-mat/9807100 and https://doi.org/10.1063/1.477415 . The published estimate concerns P3M; none is imported as an SPME-analytic error certificate.

The fixture grid, numerical budgets, bounded helper design, timing tie rule, operation caps, and oracle-convergence policy above are proposed engineering/scientific protocol decisions, not claims copied from these papers or results already observed.


### Controller adoption clarifications

The offline convergence helper is the B05 backlog's stated alternative; full
B05 still requires the actual opt-in production phase profiler and its separate
campaign. Selected configurations remain empirical for the frozen domain.
Canonical example APIs return CalibrationFailure, including the section6
interfaces previously sketched as Result<_,String>; detailed section7 statuses
remain authoritative. Runtime hash subprocesses invoke shasum directly without
a shell or RTK dependency; agent-shell commands retain the RTK prefix.
Early Tasks1–3 remain example-only and require source identity plus separate
checkout. Final B04 base synchronization/interface validation precedes Task4.
Physical equations, fixtures, targets, reference ladder, hold-out rules, budgets
and timing procedure are unchanged from the prepared specification.

### Execution checkout and resource reservation

The controller verified that the completed B01 temporary checkout was clean and
its agents terminal, backed up its historical HEAD as
`ricardo/b01-temp-pre-b05-20261106`, and reused it on
`ricardo/b05-spme-calibration` at B03 `12e4ab1`. Its path is
`/private/tmp/modena-b01`; the directory name is historical. Existing ignored
artifacts and shared target remain intact. The separate native managed-worktree
creation remains pending, without restart or a claim of terminal failure.
No performance lease is active. Early Tasks1–3 are independent example-only
work; final B04 synchronization and source/interface verification precede
Task4. All team compute/builds/writes pause during any exclusive campaign.

### Task 1 concrete preparation interface and source policy

Frozen before any oracle/SPME force calls: `prepare` writes one family's
`protocol.json`, `snapshots.json`, per-ID `fixtures/*.json`, and `manifest.json`
into a new directory. Existing directories/files are never overwritten.
Hashable physical JSON has no timestamps. Individual snapshot hashes cover
pretty JSON plus one final LF, after wrapping and canonical pair policies;
collection hashes cover the same serialization convention in fixture order.
Wrap regression IDs retain the axis/sign of the canonicalized full-box probe;
the later physical regression must also apply that full-box shift to the raw
evaluator input, rather than interpreting wrapped arrays as the raw probe.

The sorted source bundle is root `Cargo.toml`, `Cargo.lock`, core `Cargo.toml`,
all actual `crates/modena-core/src/**/*.rs`, the actual example entry point,
and all actual `examples/spme_calibration/**/*.rs`. It intentionally invalidates
old certification after source/base changes, including later B04 integration.
The protocol Markdown is hashed separately from that physical source bundle;
protocol and snapshot JSON are hashed separately as well. The historical seven
controller interface checks and historical section 8 research reads are not an
exhaustive runtime source bundle. Tests and logs are report evidence, not physical
source dependencies. Task 1 preparation is not a build/binary/source binding
certificate; later launch-freeze work must establish that binding.

Controller ruling replaces the planned synthetic force-call counter with actual
prepare contract tests and review of the pure module boundary: the preparation
code imports the lattice builder and has no evaluator, force-model construction,
or compute path. It defines future serialized records without fake physics
implementations. Numerical phases will recheck this boundary separately.

### Task 1 materialized SHA-256 roster (before force evaluation)

Default reference ladder, one thread; pretty JSON plus final LF.
The snapshots hash is independent of thread count and extended-reference choice.

| Family | Protocol JSON SHA-256 | Snapshots JSON SHA-256 |
| --- | --- | --- |
| cubic-64 | `31df92c0d443c2b8e22a6c171bd5ded0a7f882b5fa9c39e9fbb5c7b691003817` | `b3c8833c77185a93d799b0f2eee83df95695e9bf007cbd1b9033b3f943ece514` |
| cubic-216 | `01968a5a4a7e670e472c13a8bc6b0ba6319cb0cd0dea9e1f3af366f12d282491` | `13de8ab7275ef2c317cbe73ebfac405c9030f5dacd6d3b28bf29b81cb622f3b9` |
| ortho-64 | `6d08dbda363728d10359a0a50847fc3166d36eed1d97ac5600831c6d1be2697b` | `9a1f626df4aabe9338d6f31f5a41f938de24fee65826c10ea633fce9c347eefe` |
| ortho-216 | `0a9d3931b80cfbcd36775605494659211bd7ce5aba17e56b53e0b04a3ad93789` | `9c43eab24ee2e13a035772d30f487923c66c7d8673bfaba51a1fc9ce2b5e3dc6` |
| special-3 | `446238deeda0dfe273273f5e1c4b4800afada80a00884f5319a1385a694f0fd8` | `43c697afef02cf27c9f22f2acb69fa3208bbbda564c6d345485928488593cc62` |

Individual canonical snapshot JSON hashes:

```text
54626050ea3af6c6b6a62ce0e39776e781fc18e96020fa9fbd870594af6cf034 cubic-64-background-seed-19
e018779cc665f7a5a3113fcaffd38519f70670302c6dc61ff3b49406f94588d0 cubic-64-background-seed-59
af2f2b53f25da4783a08c31189f9d04c3a1a4467e6fe8db0a2261fbb6bf2f547 cubic-64-crystal
4c5487f8f7017d62bd6243edde2fd4c9f097df9ab9d8d73bfc5757ba9683f306 cubic-64-neutral-seed-19
3af638f9c9008d9a87a618d39af6bc53dc56871c323a680f2e8d098df444d9ed cubic-64-neutral-seed-31
f182d82943d76856235b124747bc03ab4764c431a19eb2a09359493953db93fd cubic-64-neutral-seed-43
eedd6d22a2ac3649c8a5034c411cbec2194ff3e533ff2ea533afcdc596987167 cubic-64-neutral-seed-59
3bbca1d7cd35bd3bb7ba0ae3495fba1170d75e08de91fe5b81f78294316074e8 cubic-64-neutral-seed-7
840c65a7a832910de1a72ade8c4819b580b6745748a2a524d9444a6f8ba567e9 cubic-64-neutral-seed-7-wrap-0-minus
ac8db954cb266c031c0ce82ac5ea86ac888a7e931661d176bad6a0116f0c5587 cubic-64-neutral-seed-7-wrap-0-plus
d539d53da2ec18ed8257983259dc9bc21128753188debd67e3ca9a49be0c63f0 cubic-64-neutral-seed-7-wrap-1-minus
1423e7f9d8d415782f8ad9ed8764625742f40ef64bad29e9791cebe8077b6c5c cubic-64-neutral-seed-7-wrap-1-plus
8e62cba1034df350887751993b6b2b3198192c54afc443382ad9a825f3acca34 cubic-64-neutral-seed-7-wrap-2-minus
a70846a2e5a549ca9e6715c501972cf149611f2315ff31ee6db0c8d6096bb3e2 cubic-64-neutral-seed-7-wrap-2-plus
1b78fe459b665c01c8aea7cc7660bcbf813bdbec3d2440fe1eec35cb59a21d55 cubic-216-background-seed-19
2bba3b29277fd7920eea10fed4799b4eb727c344cbea8e5023782291d4da01fa cubic-216-background-seed-59
304d1d9ef5da4a55dbc1d1214022ef1ce62b3e18af2168cd66cf232391cf6e61 cubic-216-crystal
8b7ac8a44dde676698f5ee7abe5d0b711fed91a743d45a9a4bbdb89615d52bff cubic-216-neutral-seed-19
f4f8dae04af9ab4e660b2cc65e67ff6ee3037903f3702b5e1d3e53d98702bc32 cubic-216-neutral-seed-31
4b3fbb65c8023f5b323984c1de38ee0c71c89a0733673e6b68c0b4cb6243bac1 cubic-216-neutral-seed-43
d8d7d51b398c99a243c95fd78c649cf011ba3d65d576c860b9961b60a9bd15a8 cubic-216-neutral-seed-59
644ba481a9553913bf9c7ab6f513d1fa771c0b07bf6e101030eea0b4b3d94edd cubic-216-neutral-seed-7
e681eded1ef7bcd961dc6928dd17bdc6a18d3e707db1cfa7dda6b3a6e6a4e6b8 cubic-216-neutral-seed-7-wrap-0-minus
39d302992d9c0455b15a9223e30b16411c92459526e2146d7e7f3453b5dbfd15 cubic-216-neutral-seed-7-wrap-0-plus
f1ecab4c6c11aaa55ea0352730626b200e002510757e16cad8fcb15f77a82326 cubic-216-neutral-seed-7-wrap-1-minus
8f2264e97099315dc42ea1fc4c9d1d17f57d920f5cef218268c829d4bc481bac cubic-216-neutral-seed-7-wrap-1-plus
ce5648e8a4ce05d84d35f911a998c2c4136767f037499c9d72f9f68d61a5994d cubic-216-neutral-seed-7-wrap-2-minus
8136536a03e0930a4735abd8e9094aa8a401fda4dd810886d36c7b5405e226b9 cubic-216-neutral-seed-7-wrap-2-plus
247a927795d724bfa022af406894c0f659ebdfa6c92fe9c9ed41102a7b909a67 ortho-64-background-seed-19
095b5412594f3808baeb9cf45fbcb7e2bebb97fc2b40f202c820cf933229d44a ortho-64-background-seed-59
7a059b1d96bc7518cd74db6590dd5ebd65d2b1d01396ee31fd4acc9700ac78a8 ortho-64-crystal
e5945a66d5362ea2f2615b8cb7f7c6fcff7de2f19b430a65b2a0922b75b74e0a ortho-64-neutral-seed-19
ec4439fe43b0bedf07137dbf9ffdc6a34476a60abbd040a29262303593c35485 ortho-64-neutral-seed-31
bed41ef42fb7fd762d742cbf991865b1a46772c7142c3d726727a8674f489d7b ortho-64-neutral-seed-43
2e3c519b4d2c5d4abcf90fa1c01c5374d26b2860532fe9a860f0d2b77a563e97 ortho-64-neutral-seed-59
9fc5661dba87135a66172525720b8f3dffbbc64d546080d070d017067c1f6f44 ortho-64-neutral-seed-7
f15ae42ca850dca611eea28188e0ef8377ff7370938ea1477dfb1b36a8db1f14 ortho-64-neutral-seed-7-wrap-0-minus
ff49cffad918f56a4a0fa9cfe64f368216a7dbb2664cfcda33a88703a6ebe8a9 ortho-64-neutral-seed-7-wrap-0-plus
822e132109700006eb098dd5a7a30ab4226eed3371b2b1ff0aa75b77f28713fd ortho-64-neutral-seed-7-wrap-1-minus
d786b7b9eb96d0a933bcb828d1cca014541b59712ac2ff696147e66b0a853a19 ortho-64-neutral-seed-7-wrap-1-plus
b39fab9076f96880c6e60af22d44ed2bfc7ac64da07c45d8ad72bd06dd0fc82f ortho-64-neutral-seed-7-wrap-2-minus
25ee0f76cc4a5480ffbb3ac9acde7c0d85496f829fc2cc1030135ea309234362 ortho-64-neutral-seed-7-wrap-2-plus
e3fcc43ab90ba6c65c6bb6182b9efa2d4eaaa64e631ba196b00feb7b9ca9259b ortho-216-background-seed-19
76cf5966c78ee71de268ffeba45bade1a152cd00fc38a979b9b2150bafb6b1bd ortho-216-background-seed-59
07068df58a2178ae2b96bbdaad3838d22f50cbba38695b0cb5c61e44c112c232 ortho-216-crystal
88a58d9af101e014619b7e44802e96a5edfe3b2cdd9cddc643e7c57d95f33fd5 ortho-216-neutral-seed-19
7d0dfc93b9fab73882dc4e093fecc0e4a9750e22006b90b3180fb2901dd8fda5 ortho-216-neutral-seed-31
872776dbc98b5bfc4c41c4a0c7c911de2d6443ae0962eb703753b09f8352f9cf ortho-216-neutral-seed-43
a0916936eb518ad44fc86c48ea620ba7802d70f234dfa0d81539974d831dfcc4 ortho-216-neutral-seed-59
f1258a8588d61f25eb7e9d4564247f6f6a58aa97ef411625d4943a05e67364b6 ortho-216-neutral-seed-7
c471bed189e4a5335250701d269803cdb356a824aa946643dda72d45c6cd9f3c ortho-216-neutral-seed-7-wrap-0-minus
cd4eb339ddb2a16b0a161b5162e6e123c0f215193a20aaad1c216f9dfc94f4c9 ortho-216-neutral-seed-7-wrap-0-plus
5d2e53a58282214a3974c5a91bb4b934aa928dfa4f92a561c8088b33fe07f9ec ortho-216-neutral-seed-7-wrap-1-minus
e0cfa7ee2b8f9a92033972c116fa09c31c129be94296094ed525322ceab178bb ortho-216-neutral-seed-7-wrap-1-plus
59969f053394688725ff96a444c4fb536fc986cb648fec439aaf2e4fe1acb463 ortho-216-neutral-seed-7-wrap-2-minus
dd4f0459e8ca49ad4f0bf2fe2859eb70991f8e6f505fe44f8e9fa5774e9018bd ortho-216-neutral-seed-7-wrap-2-plus
0051f7f566befa725ed93e1e1066762d64b908e8504210db9386d0929abd8acf special-3-confirmation-charged-bare
b0e9ec135f70f1d45f443c6d796870964166dfd969840149cda433ab10a8f023 special-3-confirmation-charged-c07
9cd003f6ce98a08f53de0e32a03b3599e54d7009a6f57fd37ec5338d1cab95cb special-3-confirmation-charged-c1
d6127856805f1e25c830c69c11bf0a644d606ea2c3b8f723e90d534c083ba5f0 special-3-confirmation-charged-scaled
5005c4edf026dc4048cd25bd895a7fe44072ffd3285529ac379f112dc17adfc9 special-3-confirmation-neutral-bare
f6fdcaa7737c65ffc052b31de7b0c9d7e0793b057d28f6eb63d983bb82780a3b special-3-confirmation-neutral-c07
3f6e021fe5ed70e603d61a29c9918414e242cd4ee8d6cae06e709e290b638e13 special-3-confirmation-neutral-c1
67950f0cf4f4cd323057d6362fa741ebc57c90be6cd6f0ab953c114ccde12ba6 special-3-confirmation-neutral-scaled
56961c3ff8d6aa74207653342296452583434643abb1d34633ddc36cbf3f6927 special-3-selection-charged-bare
9ad3e5dd63647beb58acc0eb23470feb5ddd6d6ff1f102d1dd6fd7aa4da78816 special-3-selection-charged-c07
5f4c6b54dbf80fb2c692c343f76078bf99221c4192d72e07e1149faf6a4e70ea special-3-selection-charged-c1
53af84ef194811f077f77ef5836efe7eeac03625ad6de70d6a834435beb3fba1 special-3-selection-charged-scaled
695770d2ec23535ad6eb1e77fde3bc9a0f1effc345602ad980ed25310699795a special-3-selection-neutral-bare
459bcf435f18b9b77c8f4158ee3f479f6bc0336271fd8041158fecefa1922543 special-3-selection-neutral-bare-wrap-0-minus
9228307d4619c0a606712682eda32c2bb8d1d081afd37d87568faa057c6ea0d2 special-3-selection-neutral-bare-wrap-0-plus
1e24f7e7d92fa7a9af71f6670a7853b0a2256757ca44ee451d0bdd7abefd8296 special-3-selection-neutral-bare-wrap-1-minus
0de08bb53658bd48b0ed9fa2142dd63c2e4966aacd6b86db2da3e49c64720db8 special-3-selection-neutral-bare-wrap-1-plus
11a78c4eff024756ac194a1f4e189f564fa14bb864a3faed5eb930844e831532 special-3-selection-neutral-bare-wrap-2-minus
be39d5214fe77b09ee361167c7bdb540942d1e678caf537e7264c2fe24caa42f special-3-selection-neutral-bare-wrap-2-plus
edc66f4a206e2d79d94cd7af707a7fb34348b4ab7fb95d7074dc4e4081c5d233 special-3-selection-neutral-c07
3d93e082d213eb586a09852cc38f61fbb81f86a2ed0e1a57eb2473153a8fc1b0 special-3-selection-neutral-c1
d6352e2752f50a0684ab61d8f1d3694e3c69585149355ade22715c8b94a31795 special-3-selection-neutral-scaled
```

No reference, candidate accuracy, timing, or profile outcome was observed when
this materialization roster was recorded. These are input identity records only.


### Task 2 derivative geometry registration (2026-11-06)

The independent oracle derivative regression uses the fixed three-ion rectangular
selection fixtures (neutral/charged, bare/scaled), kappa=1.1/2.2, R=8.73 and
K=3.75 in absolute reduced units. Force steps remain h=1e-5*a and h/2; isotropic
strain steps remain 1e-5 and 0.5e-5. All original residual thresholds apply.
Exact real-image and reciprocal integer membership sets are compared before
perturbed energies; a changed membership fails the derivative test.

The first trial registered K=3.71 but its strain membership precondition failed
before either strained energy was evaluated: the four modes (±5,±2,0) have
|k|=3.7099640827968456 and cross K under lambda=.99999. This failed raw trace is
retained. The controller separately registered K=3.75 using geometry only,
without energy/force/W evaluations: a covering [-12,12]^3 reciprocal cube has
872 members at all five lambda values, with minimum cutoff margin
0.0036247514070448794. A covering [-3,3]^3 real-image cube has 14 ordered members
at all five lambda values for selection and confirmation geometry; minimum
real margins are 0.05911002364431539 and 0.024110251081122414. The controller
receipt is `controller-fd-geometry-preflight.json` in the Task 2 evidence directory.
This repairs an invalid regression stencil; it changes no convergence ladder,
candidate grid, fixture, step, cap, or physical error threshold. The first raw
failure log is retained, but its original source epoch was not saved and is
unavailable. Later registered source epochs and their raw results
remain separate; no historical source was reconstructed.


### Task 2 separately registered extended Madelung regression (2026-11-06)

The default four-level cubic-64 crystal sequence failed the two-consecutive-pass
prerequisite: the (5,10) to (6,12) W/atom change was
1.746796041146581e-9, above 1e-9*E_unit=4.545454545454545e-10; only its final
refinement passed. This remains a `ReferenceUnconverged` negative regression,
with the original raw run retained as `madelung-first.log`. It is not accepted
using the approximate Madelung constant or its loose sanity tolerance.

Before any extended-reference or candidate evaluation, the controller authorized
separate immutable cubic-64 and cubic-216 prepared protocols with the existing
`--reference-extended` option: exactly (4,8),(5,10),(6,12),(7,14),(8,16),(9,18).
The two consecutive strict differences and both final cross-kappa checks
(0.9/a and 1.3/a), original caps and tolerances remain required. Their prepare
artifacts are separate `extended-cubic-64` and `extended-cubic-216` directories
under the Task 2 evidence root. No additional refinement level is allowed.
Any positive Madelung sign/factor sanity result belongs only to this extended
regression epoch after those strict prerequisites actually pass. The default
campaign protocol and the later candidate campaign remain unchanged.


### Task 2 observed oracle evidence (2026-11-06)

The independent streamed all-image/full-±k implementation passed 72 force FD
component checks (four fixed selection variants, all nine components, two
steps), maximum absolute residual 1.3239451202018415e-11, and eight isotropic
strain checks, maximum absolute residual 4.4044601299475516e-11. At the registered
R=8.73/K=3.75 stencil, all exact geometry sets and actual streamed membership
fingerprints agreed: 14 included real tuples and 872 reciprocal vectors. The
independent pair-difference half-space reciprocal sum agreed with the full sum;
maximum energy/W differences were 5.551115123125783e-17 and
1.1102230246251565e-16. Periodic self-image, background, bare special-pair,
rigid-translation and raw signed full-box-wrap regressions passed their original
gates. Removing self images in an explicitly induced test mutation caused the
self-energy gate to fail by 0.1633567233; restoring the image loop passed.

The separately prepared six-level crystal regressions passed the final two
strict refinements and both final cross-kappa checks for 64 and 216 ions.
Energy per ion pair was respectively -1.58869508603016474 and
-1.58869508603016496, passing only the registered Madelung sign/normalization
sanity check. The default four-level regression remains `ReferenceUnconverged`.
The secondary core-direct-Ewald comparison passed strict reference budgets on
the four three-ion variants only after the separately registered high-kappa
full-image refinements, cross-kappa comparisons and omitted-image/real-tail
comparisons passed. Matching core/reference cutoffs were not used as an accuracy
certificate.

Stable-toolchain reference-filtered tests: 14 passed; remaining original pure
tests: 15 passed. Rust 1.90: all 29 passed. Relevant clippy (`-D warnings`) and
format checks passed. Raw failed and successful logs, exact derivative membership
sets, work/heap/deadline records, controller geometry receipt, separate extended
prepare artifacts and source epochs are retained under
`/private/tmp/modena-b05-target/b05-evidence/task2-reference-20261106`.

Work counters count actual examined ordered real tuples (including spherical
rejects, excluding only zero self) and actual atom visits in the two included-k
passes; conservative rectangular preflight bounds are separate. The 216-ion
extended regression used 218068416 real tuples and 5338656 atom-k visits across
eight calls. Per-call accounted vector/string capacities peaked at 27682 bytes
for that fixture, excluding caller-owned retained data unless supplied via the
explicit baseline. These are reference allocation records, not a process-RSS or
complete-family-memory acceptance claim. Family orchestration must retain the
budget context, caller-owned capacities and elapsed time across phases. No
candidate SPME evaluation, ranking, target acceptance, phase profiling, production
convergence workflow or full-campaign claim is made by these oracle regressions.

### Post-B04 integration and main-campaign preset, 2026-11-09

B05 Tasks1–3 were rebased onto integrated main `4fcca11527d62402d23a810a573f9bf45dc6540e` before Task4. All eight B05-owned files and the entire B05 patch remained byte-identical during synchronization; the prior approved implementation is preserved. The production SPME/Ewald/input/runner/pair-policy interfaces are unchanged. Earlier prepared/reference artifacts retain their original source epochs and are not reused as current-source certification.

The forthcoming main families at both thread counts use the already defined `--reference-extended` option, chosen before any SPME candidate outcome: `(R/a,K*a)=(4,8),(5,10),(6,12),(7,14),(8,16),(9,18)`, kappa*a1.1, and both cross-kappa*a.9/1.3 checks. The failed four-level reference remains a negative result. All work/memory/time caps, tolerances, fixture seeds, candidate tuples and held-out roles remain unchanged. New current-source prepare artifacts must precede candidate evaluations. This preset does not assert convergence or target attainment; failures still block calibration.


### Task 4 preregistration: restored B04 source and process heap guard (2026-11-10)

Before the first Task4 production-model/candidate evaluation, the controller fixed the following implementation choice. The example/test compilation installs a capped global Rust allocator. It checks/reserves requested live Rust heap before alloc/alloc_zeroed/realloc growth, rolls back failed allocation reservations, and releases on shrink/deallocation. The unchanged256MiB owned-heap cap applies to the entire example process, a conservative superset of family-owned Rust containers including private neighbor/SPME/FFT-plan storage. Its counters report requested Rust allocation bytes, not per-model capacities, allocator overhead, or OS RSS. The512MiB host peak-RSS gate remains separate. ReferenceBudget's explicit external retained vector/string-capacity baseline remains independent; global totals are never added to that baseline.

An allocation-free Unix stderr marker `MODENA_B05_RUST_HEAP_CAP` and exit86 distinguish heap-cap termination. `MODENA_B05_TIMING_ALLOCATOR_CALLBACK` and exit87 distinguish forbidden allocator traffic. A null System allocation rolls back its reservation and emits `MODENA_B05_SYSTEM_ALLOCATION_FAILED` with exit88; this is a resource failure, not a successful zero evaluation. The parent CLI must retain already committed partial records and classify either termination; Task4 supplies the guard and isolated-child proofs, not the later CLI failure-persistence implementation. A no-allocation guard entered/exited outside primary clocks rejects any alloc/dealloc/realloc callback in the guarded region. Accepted warmed ordinary paths therefore have no allocator callback/budget instrumentation inside the measured region. Cold/setup/profile observations are separate. No Task4 primary timing or selection is performed.

Oracle translation certificates may be reused only for an exact materialized snapshot hash, actual candidate mesh, chosen reference settings, and prepared manifest source/protocol identity. The candidate kappa/rc values sharing that mesh have identical probing translations. Every certificate retains all six nonzero shifted raw oracle E/F/W evaluations, deltas, audits and actual cumulative work/time. It is reusable only after all six pass both the registered translation residuals and the stricter applicable reference budgets. Failed/partial/mismatched certificates are retained and cannot authorize reuse. Project the deduplicated checks and whole selection/confirmation call roster before a full-grid campaign; exceeding any unchanged family cap yields BudgetExceeded with no selected configuration.

Task4's first bounded actual-model tests are fixed before outcomes: the six-level prepared special-3 selection roster (eight neutral/charged and bare/scaled/c=1/c=.7 snapshots), all seven shifts of one predefined grid candidate kappa*a1.6,rc/Lmin.40,kcut*a24, raw metrics at the unchanged targets1e-3/1e-4, and fixed-mesh force/isotropic-W finite differences on its charged/scaled fixture at h1e-5 and5e-6. Pair-scale algebra, integer wraps, additive/reset behavior and cache changes are regressions. Separate warm no-callback children use kcut*a32 at one/four threads to exercise the production parallel-FFT branch, without clocks or performance claims. Every actual run follows fresh current-source `prepare --reference-extended`; source/settings receipts and all failures remain under the new Task4 evidence epoch. Historical default four-level failures remain unchanged. No full-grid, primary-cost, holdout candidate or selected-configuration campaign is authorized by this functional test registration.


Task4 functional evidence update (same physical protocol): the first warm-path guard probe accepted the one-thread child but rejected the four-thread child with exit87 on the production parallel-FFT branch after five warmups. This is a truthful instrumentation rejection, not accepted primary cost data. Four-thread ordinary batches that trigger any allocator callback must remain invalid and cannot rank/select a configuration. No core optimization, weaker guard, or instrumentation-inclusive primary timing is substituted. The retained negative and later callback-kind diagnostic are part of Task4 evidence; cold/profile observations and future independently valid cases remain separate.


### Task 5 ordinary-cost functional preregistration (2026-11-10)

Before any accepted Task5 release cost observation, the bounded functional smoke
is fixed to a fresh current-source six-level `special-3`, one-thread prepared
protocol and the existing grid candidate kappa*a1.6, rc/Lmin.40, kcut*a24.
It evaluates exactly the eight selection snapshots and seven shifts using the
unchanged reference/error/certificate gates, then records a separate ordinary
cost diagnostic on `special-3-selection-charged-scaled`: construction/first call,
one forced neighbor rebuild, five warmups, and five batches of twenty calls.
This single-candidate diagnostic cannot establish a complete feasible set or a
fastest candidate, and does not produce a selected configuration or use holdout
candidate results. All original targets, fixtures, grid tuples, six reference
levels, work/memory/time caps and thread policies remain unchanged.

Pure synthetic selector/confirmation/phase-hash controls exercise the complete
exact grid and all target decisions without using synthetic costs as measured
scientific evidence. An isolated release-only timing-observer failure control
uses the same three-ion candidate and stops after the second completed warm
batch when its observer returns an induced I/O error; completed cold/rebuild and
two batch totals must survive, with no accepted TimingRecord. All observers and
allocator guard entry/exit execute outside primary clocks. This is a lifecycle
control, separate from the smoke and from the final Task8 campaign.

Accepted ordinary batches retain the no-allocator-callback guard, force/W reset,
production cache checks and black-box consumption of raw U/W/force slices.
Incomplete/nonfinite/rejected costs cannot rank. Quartiles of the five retained
batch costs are sorted order statistics at indices1 and3, with median index2;
raw totals also retain observed minimum/maximum without another estimator.
The 5% tie rule remains relative to the global smallest median, followed by
actual mesh product, rc, kappa/kcut and finally stable ID for identical keys.
All target IDs freeze before any heldout candidate call; a confirmation failure
retains its frozen ID, blocks that fragment and never triggers a fallback.
Source-bound reference reload requires the current prepared manifest identity
and the authoritative phase-file hash from the orchestrator's trusted journal,
plus restored cumulative work/calls/elapsed. Revalidating a prepared directory
alone never gives imported reference records a new origin.


## Task6 opt-in profile preregistration (2026-11-10)

Before Task6 physical/profile observations, fix unchanged analytic-force parity tolerance per scalar/component at 1e-12*max(1,abs(value)). New focused tests use the exact special-3 selection charged/scaled coordinates and cubic-64 background seed19 amplitude.05 recipe; cold/repeated calls, nonzero force/W preload, additive comparison, coordinate/box/kappa/mesh cache transitions are diagnostic controls, without duration thresholds. Required scoped cache/hardening/SPME suites follow final implementation once. Actual example smoke, if run, must first prepare fresh current-source extended-six fixtures and retain source hashes under a unique task6 evidence directory. Ordinary ranking, allocator guard, fixed caps and physical formulas remain unchanged. No full campaign is authorized here.

Task6 controller clarification, before profile observations: phase diagnostics use exactly the first Selection snapshot in validated prepared order, unshifted, matching ordinary primary timing (cubic/ortho-{64,216}-neutral-seed-7 or special-3-selection-neutral-bare). Baseline and each distinct confirmed target use the same representative. Each model consumes106profile calls: cold1, warmups5, and five20-call batches. At most318calls fit the unchanged1000cap. These are representative-only phase observations; full selection/confirmation snapshot and translation accuracy gates remain unchanged. The authoritative --campaign profile CLI adapter is Task7-owned alongside its summary/journal schema; Task6 provides the bounded runner, individual records/writer and completed-call observer. No temporary provenance schema substitutes for that chain.

Task6 baseline control roster, fixed before outcomes: explicit baseline parameter tuple; first Selection representative; central evaluation plus nonzero full force/W preload additive evaluation, then first3 atoms x3 axes at central-difference h=1e-5 and5e-6 (36calls), plus isotropic strain at the same two steps (4calls), total42accuracy calls. Force/W analytic-FD absolute residual is the existing Task4 2e-8 bound; additive tolerance is per-component1e-12*max(1,abs(value)). All real pair memberships and actual meshes must remain central; a boundary stops validation, never relaxes a threshold. Derivative coverage is all special3 atoms but only first3 ionic atoms; full-component additive/profile parity still covers all atoms. Retain raw energies/forces/W, signed label/step, mesh/pair identities and derivative residual inputs. Continuing accuracy/pair/deadline caps apply before and during work. Baseline diagnostic ID126 is reserved outside selectable0..125 and never enters ranking/freeze. A bounded special3 release smoke validates those42controls then records106profile calls, followed by at most one cold profiled call whose completed observer deliberately returns Io to test partial retention; ordinary costs/ranking are not run. Missing/source-mismatched/duplicate tokens and upfront profile budget rejection are no-compute controls. Actual smoke requires fresh current-source prepare with six-level preset; no oracle rerun or full campaign.


### Task7a runtime physical-verification extraction preregistration (2026-11-10)

Before runtime observations, the existing registered test controls are extracted
into bounded example-only APIs. Special3 uses the exact neutral/charged,
bare/scaled Selection roster and the already registered derivative geometry
kappa1.1/2.2,R8.73,K3.75; steps, full/half-k, background, self-image, bare-fraction,
translation/wrap, secondary core-Ewald ladders and strict tolerances remain
those of the existing integration controls. No candidate timing or target
acceptance occurs. Cubic64 and216 crystal sanity run separately in their own
family verify phases with the prechosen six-level preset; strict absolute
convergence may have RelativeUndefined status and cannot certify a force target.
The optional comparison of both current-source crystal records retains1e-9/a.

The runtime observer journals source/fixture-bound admission reservations before
risky calls and raw completions before later work. Admission reservations replace
actual usage only if a child exits in flight; they are not added again after a
completed event. Oracle and auxiliary work use the same continuing ReferenceBudget;
independent membership/half-space loops charge conservative rectangular bounds
and checkpoint within1024terms. Secondary core Ewald uses bounded3atom work with
before/after deadline checks because production has no inner callback. An expired
post-call checkpoint prevents a passing proof; this is not an uninterrupted
core-call deadline guarantee. External launcher deadlines remain required.

Special physical prerequisite usage belongs to special3 at its own thread count
and must be restored when that family reaches its frozen campaign roster slot.
Cubic64/216 controls belong to their corresponding continuing family budgets.
No uncharged global pool, source relabeling or status-only physics proof is allowed.

Task7a smoke roster is fixed before its observations: four special Selection
fixtures each have42derivative oracle evaluations,42independent membership
enumerations, full/no-mode/half-space comparisons, one background evaluation,
eight translation evaluations,19raw-wrap evaluations, six secondary reference
evaluations and one bounded core call. Each of the two bare fixtures adds
nine bare-fraction evaluations; one self-image fixture adds one oracle call.
Totals:331oracle evaluations,172auxiliary enumerations/sums,4core evaluations;
reference-call accounting includes auxiliaries (503), model accuracy calls4.
Each cubic crystal has six ladder levels and two cross-kappa evaluations,
charged to that separate family's existing budget. Whole planned work and
retained capacity bounds are checked first; cap rejection is retained honestly.


## Task7b immutable integration contract and preregistration (2026-11-10)

This append-only integration clarification precedes Task7b current-source observations and the Task8 final campaign. Factual final results will live in [B05-SPME-RESULTS.md](B05-SPME-RESULTS.md), with the delivered domain recorded in VALIDATED_SCOPE. Those result documents are outside execution source identity. Do not amend this hashed protocol after the final measurement merely to add result links.

The exact frozen campaign roster remains cubic2, ortho2, special2, cubic3, ortho3, first at one thread then at four threads. Every launcher prepare explicitly chooses `--reference-extended`; generic prepare still supports the default four-level ladder. Special3 physical verification runs once before ionic work at each thread count. Its verification receipt, actual counters, and active duration return to the same special3 slot later in this roster; idle time spent on other families is excluded. Cubic64 and cubic216 Madelung work occurs in their own verify phases, and cubic216 consumes the authoritative cubic64 proof for replication comparison. Orthorhombic verification never applies the cubic Madelung formula.

Baseline126 is nonselectable. It is evaluated by the same role/certificate/row engine on every selection and confirmation snapshot times seven shifts, with held-out baseline work only after both target IDs are frozen. Its ordinary107-call cost uses the same first Selection central snapshot. Whole-campaign admission includes these baseline calls, a potential third confirmation mesh,42 baseline functional controls, and up to318 diagnostic profile calls. The baseline kcut12 selection certificate may reuse an exact grid mesh. No selectable grid entry, cap, tolerance, formula, or roster is changed. Profiles use baseline plus each distinct confirmed ID only,106calls/model on the first Selection central snapshot; target exit4 does not block the eligible baseline. Reference/physics/resource failures stop all further computation for the affected family.

The binary accepts exactly prepare, verify, sweep, confirm, profile, report. `cost-smoke` is no longer a production mode. Six-phase usage, with absolute paths recommended and hashes supplied from the trusted parent journal, is:

```text
spme_convergence prepare --out PREPARED --family special --replicas 2 --threads 1 --reference-extended
spme_convergence verify --protocol PREPARED/protocol.json --snapshots PREPARED/snapshots.json --elapsed-ns ACTIVE_NS
spme_convergence verify --protocol IONIC/protocol.json --snapshots IONIC/snapshots.json --special SPECIAL/verification.json --special-sha256 TRUSTED_SPECIAL_SHA --elapsed-ns ACTIVE_NS
# cubic216 additionally requires:
# --madelung CUBIC64/verification.json --madelung-sha256 TRUSTED_CUBIC64_SHA
spme_convergence sweep --protocol PREPARED/protocol.json --snapshots PREPARED/snapshots.json --verification-sha256 TRUSTED_VERIFY_SHA --out SWEEP --elapsed-ns ACTIVE_NS
spme_convergence confirm --campaign SWEEP/summary.json --expected-sha256 TRUSTED_SWEEP_SHA --out CONFIRM --elapsed-ns ACTIVE_NS
spme_convergence profile --campaign CONFIRM/summary.json --expected-sha256 TRUSTED_CONFIRM_SHA --out PROFILE --elapsed-ns ACTIVE_NS
spme_convergence report --campaign CONFIRM/summary.json --expected-sha256 TRUSTED_CONFIRM_SHA --profiles PROFILE/phase-times.jsonl --profiles-sha256 TRUSTED_PROFILE_DATA_SHA --profile-receipt-sha256 TRUSTED_PROFILE_RECEIPT_SHA --memory FAMILY/memory.json --memory-sha256 TRUSTED_MEMORY_SHA --out REPORT
```

For a family failure, report consumes the last retained failure receipt with the same campaign/memory/hash/out flags and omits all three profile flags. Compute requires a release executable and exact matching `RAYON_NUM_THREADS=1` or4. Unknown, duplicate, missing, or malformed options are errors. Invalid JSON is exit2. Source-path/hash I/O failures are exit1. Failure exit3/2 stops a family; target exit4 retains decisions and permits unaffected confirmation and eligible profiling.

Supplemental `modena-spme-phase-v1` receipts are written as `verification.json` beside prepared files and `summary.json` in later phase directories. They contain the current prepared manifest digest, prior receipt path/SHA, cross-family proof inputs, individual artifact hashes, cumulative reference/candidate/timing/profile usage, and primary/finalization failures. These are supplemental decision manifests, not replacements for canonical scientific records. The report directory's `summary.json` is canonical CampaignReport. The external parent expected hash anchors each imported chain; an imported file's self-computed hash alone is never authority. Checked special/Madelung loaders enforce successful complete physical invariants in addition to this chain. Failed receipts remain raw-reportable and cannot mint proofs. All source inventory includes the Python launcher bytes, Rust example/support/core and Cargo bytes, plus the separately hashed protocol specification.

Every reference level, translation reference, candidate call, ordinary stage, baseline functional control, and profiled call journals admission and completion before later risky work. Events carry the current prepared identity and sequence. Files are create-new and synced. Admission is a conservative cumulative reservation; completion replaces it with actual cumulative usage, never adds both. An allocator/process exit preserves an unfinished reservation and raw completed rows without inventing an evaluation. Parent recovery retains a torn tail raw and creates a separate complete-line prefix for reports. Primary, nested reference finalization, and outer finalization diagnostics remain inspectable. Exit86=requested Rust heap cap,87=forbidden primary allocator callback,88=null System allocation; these are runtime failures, never scientific target unreachability.

The standard-library launcher invocation remains:

```text
python3 bench/spme_calibration_campaign.py --out /private/tmp/modena-b05-campaign-v1 --threads 1 4
```

It builds with `cargo build --release -p modena-core --example spme_convergence` using absolute `CARGO_TARGET_DIR=REPOSITORY/target/b05-build` before acquiring its measurement lease. Each measured command runs the resulting exact hashed Rust release executable under `/usr/bin/time -l` and the launcher's direct-child wait4 monitor. This separates actual Rust exit from Cargo build, monitor completion and time-wrapper status. The monitor records actual wait status and child RSS; numeric macOS time-wrapper peak RSS remains mandatory for512MiB acceptance. The requested Rust heap cap remains256MiB. Parent timeout receives the remaining1800s cumulative active-family budget. Every child's command, environment thread count, stdout, stderr, status, RSS and active duration is retained; missing/over-cap RSS prevents FullScopeMeasured. Release source/binary/launcher/protocol identities and exact flags are registered before compute and checked throughout. No shell execution, overlapping benchmark child, or automatic cap/grid reduction is permitted.

Task7b's single actual diagnostic is preregistered as special3, replicas2, thread1, six reference levels. After successful source-bound physical verification, a test-only diagnostic resumes the same family counters/time, generates all eight selection reference ladders, evaluates only explicit baseline126 on all56 selection rows using48 shifted reference calls, measures its107ordinary calls, runs42functional controls and106profile calls on first Selection central. No heldout, selectable grid sweep, frozen winner, selected fragment, or FullScopeMeasured claim occurs. Expected completed counters are615reference calls (503physical+64selection+48shifts),102accuracy calls (4physical+56baseline+42functional),107timing and106profile calls. Source registration records exact binary/test executable hashes, prepared/receipt expected hashes, environment and absolute paths before this observation. This diagnostic is explicitly NotCalibrated; generated tables/SVG display actual raw baseline errors/costs only. Existing negative four-thread allocator evidence is retained; no guard relaxation or production optimization is authorized.

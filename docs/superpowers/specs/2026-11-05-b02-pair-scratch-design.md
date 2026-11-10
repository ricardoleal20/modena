# B02: measured native CPU pair scratch reuse

The owner requested B09 followed immediately by B02. B09 is delivered in open
PR #55 at `ce23e8e78bc6b0745b113a146609c84e46332902`; B02 is a separate branch
based on that reviewed code, so both features remain usable in this checkout.

## Scope and invariants

Profile the native CPU pair path before changing its allocation/reduction.
Evaluate reusable scratch for `PairForceModel`, preserving the existing
100000-half-pair serial threshold, LJ specialization, generic potential path,
periodic geometry, scales/exclusions and additive force/energy/virial contract.
One-thread and `compute_deterministic` arithmetic/order remain bit-identical.
Multi-thread summation can change within declared f64 tolerances. B09's bounded
snapshot list must still pass directly into the shared force pass, including
zero skin; do not reintroduce an unrestricted rebuild.

No System-layout, neighbor-algorithm, GPU/MPI or new physical model changes.
Public constructor signatures remain. If private scratch makes external struct
literals incompatible, document `new`/`from_boxed` migration and use a
BREAKING CHANGE footer; do not hide that compatibility cost.

## Owned scratch and contention

Use model-owned reusable force buffers with a bounded active partition count
not exceeding useful worker/pair work. Clear force and scalar state correctly,
handle N/pool changes, empty systems and repeated use, and release ownership on
drop. Dense scratch retains O(TN) memory; high-water capacities must be reported.
An initial dense design may be compared with touched/block accumulation if the
profile or measurements require it. No global model-ID cache or production
test-only instrumentation.

ForceModel remains Send+Sync. Avoid blocking a Rayon worker on a model lock while
another caller forks/joins in the same pool. One simple bounded design is a
nonblocking `try_lock` workspace with the existing serial path on contention;
short checkout/return ownership is another option if justified. The selected
policy must be documented and tested for concurrent/nested pool calls.

## Profiling and acceptance

Before production edits, add an isolated allocator profile/regression and
preserve ordinary release CLI/profile binaries built in this worktree from a
clean baseline commit. Record compiler/flags/source/binary hashes. Profile
actual pair count, force-sized alloc/realloc counts/bytes and memory across N
and threads after pool/list warmup. A stable, unreconstructed, unresized force
pass must allocate zero N-sized force buffers after warmup. Tests must show RED
on the old parallel path and GREEN on the candidate. Use test/example-only
allocator hooks. Existing `tests/force_cache.rs` provides the isolation pattern.

Check real half-pair counts: N2048 can still be below the serial threshold.
Use N864/2048 controls and N5000/10000 parallel cases, threads1/4/8, with recorded
geometry and pair counts. Include composition with nonzero accumulators,
scaled/excluded pairs, generic dispatch, resizing, repeated calls and shared
models. Compare E/F/W against the unchanged serial/scalar path and independent
potential consistency tests; declare tolerances before comparing timings.
For large static comparisons use E/W absolute tolerance1e-8 plus relative1e-12,
and every force component absolute1e-10 plus relative1e-12; one-thread results
require exact equality. Do not relax these bounds merely to accept a candidate.

Measure ordinary release *full CLI* before/after binaries, with matched LJ
sigma=epsilon=mass=1, cutoff2.5, skin0.3, FCCa1.6796, dt0.005, T1, seed87287,
unitsreduced and output policy. Alternate AB/BA per paired repetition. Start
with1000steps,5pairs per case and separate discarded warmups. Record actual
steps, wall/user/system CPU, RSS, load, raw logs/inputs and every sample. Use
longer/repeated selected cases if observed noise prevents a conclusion; never
select only favorable samples. Instrumented allocator timing is profiling,
not the primary speed claim. Preserve one-thread physical results exactly;
multi-thread static E/F/W and scientific C1/timestep gates must pass.

Retain a variant only with reproducible runtime or measured memory benefit in
its declared case, with full-step cost and neutral/regressing cases explicit.
Eliminating malloc alone is not evidence of a useful speedup. If dense reuse
does not help, evaluate a focused touched/block alternative; reject ineffective
production changes rather than manufacture a performance claim.

## Verification and delivery

After source stabilizes, run one full actual-ASE workspace/Python/worker/tooling
campaign, formatting, default/GPU Clippy, Rust1.90 default/GPU and the required
standalone-harness check. Run67-observation scientific smoke. GPU checks are
compile/lint only; no hardware/MPI/other-engine claims. Preserve B09 evidence.

Read-only task and broad final reviews precede PR publication. B02 PR targets
`ricardo/b09-ase-sessions` while #55 is open; clearly state dependency and
eventual base synchronization after B09 squash, without force-pushing a
published branch. No new merge has been authorized.

## Global constraints

- Branch prefix `ricardo/`; preserve unrelated files and prior evidence.
- Rust floor1.90; no new mandatory dependency or downloaded model.
- Shell commands start `rtk`; literal gitmoji plus imperative Action commits.
- No force-push, protected-branch mutation or merge without authorization.
- Finish B02 in this task without stopping to ask whether to continue.

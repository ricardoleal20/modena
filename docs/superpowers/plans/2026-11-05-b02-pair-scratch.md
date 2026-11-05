# B02 CPU Pair Scratch Implementation Plan

> **For agentic workers:** use subagent-driven-development with one implementer and separate task review at a time. Preserve the B09 implementation/evidence and unrelated edits.

**Goal:** remove repeated N-sized native parallel force allocations only when measured full-step or memory results justify the change.

**Architecture:** model-owned bounded reusable partitions, unchanged serial path, explicit nonblocking contention policy; independent profile and matched before/after release binaries.

**Tech Stack:** Rust1.90, Rayon, std synchronization/allocator instrumentation in tests/examples, stdlib Python benchmark tooling.

**Spec:** `docs/superpowers/specs/2026-11-05-b02-pair-scratch-design.md`.

## Global Constraints

- Branch prefix `ricardo/`; preserve unrelated files and prior evidence.
- Rust floor1.90; no new mandatory dependency or downloaded model.
- Shell commands start `rtk`; literal gitmoji plus imperative Action commits.
- No force-push, protected-branch mutation or merge without authorization.
- Finish B02 in this task without stopping to ask whether to continue.

## Task 1: profile, reusable scratch and measured acceptance

Files: `crates/modena-core/src/force_model.rs`, a focused private scratch module
if needed, `crates/modena-core/tests/pair_scratch.rs`, isolated profile example
`crates/modena-core/examples/pair_scratch_profile.rs`,
`bench/pair_scratch.py`, `bench/tests/test_pair_scratch.py`, migration/limits
guide `docs/B02_PAIR_SCRATCH.md` and scoped changelog.

Interfaces: retain `PairForceModel::new`, `from_boxed`, `ForceModel::compute`,
`compute_deterministic` and B09's `compute_with_built_neighbors` contract.
The profiler is standalone; no public test hooks. The benchmark accepts:

```bash
rtk proxy "$PY" bench/pair_scratch.py --before "$BASELINE" --after "$CANDIDATE" \
  --sizes 864 2048 5000 10000 --threads 1 4 8 --steps 1000 \
  --repetitions 5 --output "$EVIDENCE"
```

- [x] Add allocator profile/isolated warmed regression first; measure actual pair counts and observe the zero-force-buffer-allocation assertion fail on old parallel code. Keep serial controls separate.
- [x] Commit profile/test/driver support before changing production, build and retain clean-source baseline CLI/profile binaries and raw profiling data. The meaningful regression is:

```rust
// After pool + neighbor warmup, while N and pool stay unchanged:
system.zero_forces();
virial.w = 0.0;
let before = force_buffer_allocations.load(Ordering::SeqCst);
pool.install(|| model.compute(&mut system, &mut virial, &mut neighbors).unwrap());
let observed = force_buffer_allocations.load(Ordering::SeqCst) - before;
assert_eq!(observed, 0, "warm pair pass allocated an N-sized force buffer");
```

- [x] Implement bounded owned scratch; reset forces/PE/W, retain serial arithmetic, handle resize/contention safely and preserve the prebuilt bounded snapshot path. Add exact/tolerance E/F/W, composition, scales/exclusions, generic, empty/repeat/resize and nested-concurrency coverage.
- [x] Run focused tests to GREEN, format/lint/default compile; self-review and commit candidate source/tests before timing.
- [x] Build ordinary release candidate in the same worktree/compiler settings. Run failure-sensitive driver tests (actual subprocess boundary, incomplete steps/physics failure/fresh outputs), then alternating paired full-CLI measurements and isolated allocation/memory profiles with provenance and fixed numerical gates.
- [x] Evaluate all cases, not just favorable ones. If dense scratch lacks benefit, measure a focused touched/block alternative; retain only a justified variant, documenting controls/noise and any regressions.
- [x] Write full Task1 report with RED/GREEN, commands, source/binary hashes, all raw rows, memory-vs-N/threads, medians and numerical policy; commit guide/results and get independent read-only task review, with covering fixes/rereview if needed.

## Task 2: complete validation, evidence and closure

Files: `docs/B02_VERIFICATION.md`, guide/architecture/validated-scope/changelog
as required, existing tooling/smoke outputs retained under ignored target.
Consumes Task1 accepted source and before/after artifacts; introduces no
unmeasured force-kernel changes.

- [ ] Run fresh actual-ASE full workspace, Python/worker/tooling tests once at stabilized source; preserve failures and correct concrete issues with scoped tests.
- [ ] Run formatting, default/GPU all-target Clippy with `-D warnings`, Rust1.90 default/GPU workspace checks and the locked standalone harness check; compilation only for GPU.
- [ ] Build scientific release example and run existing67-observation smoke. Retain versions, code/binary hashes, raw series and exact results; no tolerance relaxation or new ensemble/engine certification.
- [ ] Record reproducible commands and accurate source epochs in verification/guide, preserving raw profile/timing rows and B09 evidence. Run diff checks and obtain read-only task review.
- [ ] Controller: whole-branch review, grouped fixes/scoped rereview if needed, then assigned/labeled PR against B09 with Summary, Changes, Test plan, Refs. No merge; report B09/B02 results together.

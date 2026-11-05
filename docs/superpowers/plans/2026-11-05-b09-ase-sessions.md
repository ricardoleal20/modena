# B09 ASE Sessions Implementation Plan

> **For agentic workers:** use subagent-driven-development task by task with test-first implementation and separate read-only reviews. Preserve unrelated edits.

**Goal:** implement both B09 adapter directions and measure the effect of persistent EMT transport.

**Architecture:** one locked, terminal-on-error process session per persistent provider; separate strict native LJ snapshot CLI consumed by an optional ASE calculator. Preserve one-shot compatibility.

**Tech Stack:** Rust 1.90, std process/IO/synchronization, serde JSON, Python optional ASE/NumPy.

**Spec:** `docs/superpowers/specs/2026-11-05-b09-ase-sessions-design.md`.

## Global Constraints

- Branches begin `ricardo/`; preserve unrelated files and the old B01 evidence.
- Rust floor 1.90; no new mandatory Python/ASE dependency or downloaded model.
- Shell commands begin `rtk`; commits use literal gitmoji plus imperative Action.
- No merge of B09/B02 is authorized by the earlier B06/B01-specific exception.
- Finish and review B09 before implementation of B02; publish reviewable PRs.

## Task 1: persistent pinned EMT provider

Files: `crates/modena-core/src/external.rs`, new focused session module under
`external/`, bundled `ase_emt_worker.py`, external-provider tests; Python config
and MCP schema; affected interface tests and external-provider documentation.

Consumes `AseEmtSpec`/`AseEmtForceModel` and B01 transactional evaluation.
Produces opt-in `transport="persistent"` with the same physical identity.

- [x] Write regressions that repeated evaluations use one process, return changed E/F, and reject stale IDs/identity while preserving accumulators; run to observe missing support.
- [x] Implement bounded request/response session ownership, deadline, cleanup and terminal state; preserve B01 one-shot default. Extend worker with a session loop that constructs EMT once and validates/refreshes snapshots.
- [x] Cover EOF/provider failure, oversized/malformed output, timeout, drop/reaping, restart/physical signatures and safe cancellation; compare actual installed ASE for changed snapshots/cells.
- [x] Run focused core/CLI/MCP/Python tests and format/lint/compiler-floor checks. Record exact commands/results/source in a task report and commit.
- [x] Read-only spec/quality review; correct findings with covering tests and scoped rereview.

## Task 2: MODENA native LJ calculator for ASE

Files: new focused core snapshot module and tests, CLI snapshot command/tests,
`python/modena/ase.py`, Python calculator tests, Python/CLI usage docs.

Consumes existing core `PairForceModel` and `LennardJones`; produces a strict
Angstrom/eV JSON snapshot command and `ModenaCalculator` with energy/forces only.

- [x] Add independent two-particle shifted-force LJ energy/forces assertions and finite differences; run failing CLI/calculator tests before implementation.
- [x] Implement strict snapshot validation and existing core force dispatch; add optional ASE calculator without importing it from base `modena`.
- [x] Verify position/number/cell/parameter/PBC invalidation, failure-cleared results, absent stress, invalid input and actual CLI E/F units. Reject constraints, nonzero initial charge/spin, mixed species and unsupported cells.
- [x] Run focused tests plus full workspace/Python/default+GPU lint and Rust 1.90 checks. Record skips and limits; commit.
- [x] Read-only spec/quality review and scoped fixes if required.

## Task 3: measured evidence and closure

Files: reproducible B09 benchmark driver under `bench/`, `docs/B09_VERIFICATION.md`,
provider/calculator scope and changelog documentation.

- [ ] Build release from this branch; use pinned local ASE and alternate one-shot/persistent order over repeated equivalent Cu full runs. Retain timing rows, identity, command, hardware/source/binary hashes and E/F/thermo equality; calculate median ratios from raw rows.
- [ ] Record cold/full-run versus steady transport limits; no unmeasured solver speed claim.
- [ ] Check final source diff against design and run appropriate regressions for any changes. Whole-branch read-only review, one grouped fix wave if needed, covering tests and scoped rereview.
- [ ] Push `ricardo/b09-ase-sessions` and create assigned/labeled PR with Summary, Changes, Test plan, Refs; no merge. Continue B02 from a separate `ricardo/` branch.

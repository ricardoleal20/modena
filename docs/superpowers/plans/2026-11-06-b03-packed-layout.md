# B03 Packed Pair Layout Implementation Plan

> **For agentic workers:** Use subagent-driven-development, one implementation task and separate read-only review at a time. Preserve earlier evidence.

**Goal:** decide whether packed/cluster views improve full-step cost with correct logical identity.

**Architecture:** experimental ForceModel implementations in the nonpublished benchmark crate, exact logical mapping and independent scalar oracle; retain production until acceptance.

**Tech Stack:** Rust1.90, existing core/Rayon where available, stdlib Python driver.

**Spec:** docs/superpowers/specs/2026-11-06-b03-packed-layout-design.md.

## Global Constraints

- Rust floor1.90; no mandatory dependency; branches ricardo/; literal gitmoji + imperative Action commits.
- Agent shell commands start rtk, product subprocesses do not depend on RTK.
- Preserve earlier evidence, no force-push or extra merges.
- E/W absolute1e-8+relative1e-12; F component absolute1e-10+relative1e-12; identity/mask sets exact.
- Public System and production one-thread/deterministic path stay unchanged.

## Task 1: independent layout experiment and correctness gates

Files: create crates/modena-bench/src/lib.rs, packed_pair.rs and tests/packed_pair.rs;
add a standalone example layout_step.rs; touch its Cargo.toml only for existing
workspace dependencies actually needed. No production-core edits in this task.

Interfaces: consume public System, SimBox, NeighborList, ForceModel, PairPotential,
LennardJones and existing velocity Verlet APIs. Produce benchmark-only layout
selector Interleaved/Soa/Aosoa4/Cluster4, a ForceModel for each packed view and
a full-step example with JSON results suitable for a future Python driver.
Read the spec first for domains/gates; inspect existing interfaces before coding.

- [x] Establish focused baseline with force_hardening and pair_scales tests, saving command/status evidence under a new target/b03-evidence directory, preserving all B02 files.
- [x] Write independent scalar oracle and exact mapping/mask tests first. The key assertions are `assert_eq!(unpack(pack(logical_ids)), logical_ids)` and equality of BTreeSet `(logical_i,logical_j,scale_bits)` from masks vs neighbor half pairs. Perturb positions/box and force rebuild to catch stale masks; padded indices never appear.
- [x] Record meaningful RED before implementation, then implement layouts, safe masking/permutation, current-coordinate packing and additive force scatter. Test E/F/W numerical gates and edge cases in the spec; fail unsupported input explicitly.
- [x] Build a full-step experiment example with CLI mode/N/threads/steps/rebuild-policy/seed and JSON counters. Integrate using existing velocity-Verlet with same initial state for each mode; include packing, neighbor rebuild, cluster reconstruction/scatter inside elapsed region. Report setup independently, checksum final state and include actual completed steps/rebuild/pair counts and finite E/F/W.
- [x] Run focused tests GREEN, fmt/clippy and Rust1.90 check for changed targets; self-review, commit with required format and report exact commands/results/source hashes. No timing acceptance claim in Task1.

## Task 2: frozen matched campaign and decision

Files: bench/layout_step.py and bench/tests/test_layout_step.py, docs/B03_PACKED_LAYOUT.md.
Consume Task1 example binary JSON/counters; preserve original raw subprocess outputs.

- [ ] Write failure-sensitive real-subprocess driver tests for wrong mode/N/steps, missing/nonfinite physics, invalid output reuse and process failure. Product argv uses direct tools, never RTK. Validate fixed inputs and output identity independently of a declared checksum.
- [ ] Build and freeze ordinary release executable, source/compiler/input hashes. Pilot to fix duration, then execute the spec matrix with five AB/BA pairs and separate warmups including rebuild/no-rebuild scenarios. Run variants sequentially; retain every failed/unfavorable case and CPU/RSS/load rows.
- [ ] Apply the spec acceptance bar; repeat five pairs in each candidate favorable case. Record all cases, logical/mask correctness, packing/sort/rebuild costs, memory high water and scoped conclusion. If rejected, leave production unchanged. If accepted, controller amends a concrete promotion spec/plan before kernel integration; never silently enable new default.
- [ ] Commit evidence guide, tooling and report; obtain independent task review.

## Task 3: closure validation and delivery

Files: docs/B03_VERIFICATION.md, docs/CORE_RESEARCH_EXECUTION.md and validated-scope/changelog only as applicable.

- [ ] Run tests appropriate to the final diff, Rust floor/fmt/clippy and existing scientific smoke when production is touched; for experimental-only changes, verify unchanged production hashes and new full-step numerical/convergence checks without relabeling old runs as fresh.
- [ ] Record source epochs, all commands, limits, acceptance/rejection, then whole-branch review and its grouped fix wave/scoped rereview if necessary.
- [ ] Publish assigned/labeled B03 PR against B02 if still open (otherwise synchronize base safely), with four normative sections; no merge. Update overall execution record and continue B04 immediately.

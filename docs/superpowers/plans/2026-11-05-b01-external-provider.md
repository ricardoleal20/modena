# B01 External Provider Implementation Plan

> **For agentic workers:** execute inline with test-driven development, one verified unit at a time. Preserve unrelated working-tree files. B06 integration remains dependent on the user's administrator-exception answer.

**Goal:** Run pinned ASE EMT as MODENA's first concrete external total force provider and demonstrate snapshot parity through the actual runner.

**Architecture:** A bundled one-shot Python bridge returns explicit E/F/stress plus installed model identity. Rust validates complete snapshots/results before accumulation, keeps native force construction intact and incorporates provider identity into checkpoints. Python/MCP expose the same bounded domain.

**Tech Stack:** Rust 1.90+, existing serde/serde_json, tempfile for bounded process transport; Python 3.11+, optional ASE >=3.26,<4. Reference environment installed locally for this task: ASE 3.29.0.

**Spec:** [B01 design](../specs/2026-11-05-b01-external-provider-design.md).

## Global constraints

- Core arithmetic stays reduced; external wire uses Angstrom/eV, with explicit conversion and tensile-positive stress.
- Only standard EMT metals, neutral unconstrained fixed-cell real-unit NVE/NVT and PBC=[true,true,true]; total composition only.
- Pin installed ASE version and EMT source SHA-256; no silent result/virial fallback.
- Native default build needs no Python/ASE; real-reference tests explicitly require the local optional environment.
- Run every shell command through rtk. Keep model metadata separate from transport paths and force-field types.

## Task 1: strict snapshot bridge and provider

Files: create `crates/modena-core/src/external.rs`, `crates/modena-core/src/external/ase_emt_worker.py`, `crates/modena-core/tests/external_provider.rs`; export from lib.rs; add the already-locked tempfile runtime dependency. Add the System atomic-number field/init/resize here because the provider consumes it; restart persistence and runner assignment remain Task 2.

- [ ] Write Rust tests that import `AseEmtSpec`, `AseEmtForceModel` and `describe_ase_emt`, and assert explicit unsupported pins/PBC/elements and non-mutating failed evaluation. Run to observe the missing provider failure.
- [ ] Implement serializable spec/identity/evaluation types, strict bridge envelopes, bounded timeout subprocess and SI unit conversion. `describe_ase_emt(python: &str, timeout_ms: u64) -> Result<ProviderIdentity>` returns discovery; `AseEmtForceModel::new(spec: AseEmtSpec, factors: UnitFactors) -> Result<Self>` pins it.
- [ ] Add bridge Python tests against installed EMT for exact displaced-Cu snapshot outputs and malformed requests; run them before bridge implementation, then implement the bridge. Use direct `Atoms(...); atoms.calc=EMT(asap_cutoff=False)` as the reference, not the bridge's calculation helper.
- [ ] Validate fresh evaluate requests, all identities, property shapes and finiteness before accumulation. Test timeouts/process errors and property omissions through actual subprocess fixtures; do not expose test-only production APIs.
- [ ] Run focused Rust/Python boundary checks and commit the verified provider unit.

## Task 2: real runner, identity and restart

Files: modify `io/input.rs`, `config/validation.rs`, `runner.rs`, `io/restart.rs`, `force_model.rs`; add runner/restart cases to `tests/external_provider.rs`. Consume the System atomic-number field from Task 1.

- [ ] Add a TOML acceptance test using `[model] type="ase_emt"`, explicit pins/Z/total composition, a 4-Cu FCC structure and real units. Assert that unsupported coupling/composition is rejected before process launch. Observe failure first.
- [ ] Add `System.atomic_numbers`, backwards-compatible checkpoint defaults, and the external ModelSpec variant. Keep native physical signatures unchanged; derive the external signature from actual runtime provider identity plus explicit Z/PBC/parameters.
- [ ] Branch force construction around the pair-specific path. Native pair tables and composite policy remain in their existing branch; an external total model receives no automatic pair/bonded/electrostatic terms. Declare provider-managed neighbors explicitly, with no invented cutoff.
- [ ] Add provider-stress pressure metadata. For observed E/F/stress, compute W from `-V trace(stress)` and use existing thermodynamic unit rendering. Propagate failed-step error events for an evaluation error.
- [ ] Verify skew-cell row/column mapping, finite-difference force and strain W with spec tolerances, translation invariance, selected-Z mutation and split restart identity/mechanics through real ASE. Run native regressions, then commit.

## Task 3: CLI, Python and MCP surface

Files: modify CLI main/tests, Python config/exports/tests/optional dependencies, MCP schema/tools/capabilities/tests.

- [ ] Write CLI `provider-info --python PATH` and Python typed rendering tests; observe missing-command/type failures. Add discovery and typed `AseEmt` fields, preserving native rendering and path resolution.
- [ ] Write an MCP JSON request test for external scope/pins and provider metadata. Describe optional runtime execution, cancellation latency and lack of pair-kernel cost estimates; reject unsupported cost estimates instead of inventing zero or classical estimates.
- [ ] Expose identity/capabilities in final run metadata and live status consistently; do not duplicate native defaults into provider pins.
- [ ] Execute real CLI/Python/MCP integration against the installed optional reference environment and commit.

## Task 4: verification and reproducible usage

Files: create `docs/B01_EXTERNAL_PROVIDER.md`, `docs/B01_VERIFICATION.md`; update validated scope, README/architecture/changelog and optional ASE test instructions.

- [ ] Record exact installed versions/model hash, snapshot/finite-difference tolerances and executable commands. Distinguish direct-EMT parity from physical material validation and one-shot transport from B09 persistence.
- [ ] Run `cargo fmt --all -- --check`, default/GPU all-targets Clippy with `-D warnings`, default/GPU Rust 1.90 checks and complete workspace tests.
- [ ] Run complete Python/benchmark-tool tests, existing 67-observation smoke and optional real-ASE tests. Record skips as skips, not reference passes.
- [ ] Request independent code review under the requesting-code-review skill, fix actual findings test-first, and commit reviewed B01 changes on its separate branch. Do not claim B06 merged until GitHub confirms the merge.

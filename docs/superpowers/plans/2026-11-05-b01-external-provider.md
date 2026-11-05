# B01 External Provider Implementation Plan

> **For agentic workers:** use subagent-driven development with test-driven implementation, one reviewed unit at a time. Preserve unrelated working-tree files. The owner authorized administrator merges for B06, then B01, on 2026-11-05; new branches use the `ricardo/` prefix.

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

- [x] Write Rust tests that import `AseEmtSpec`, `AseEmtForceModel` and `describe_ase_emt`, and assert explicit unsupported pins/PBC/elements and non-mutating failed evaluation. Run to observe the missing provider failure.
- [x] Implement serializable spec/identity/evaluation types, strict bridge envelopes, bounded timeout subprocess and SI unit conversion. `describe_ase_emt(python: &str, timeout_ms: u64) -> Result<ProviderIdentity>` returns discovery; `AseEmtForceModel::new(spec: AseEmtSpec, factors: UnitFactors) -> Result<Self>` pins it.
- [x] Add bridge Python tests against installed EMT for exact displaced-Cu snapshot outputs and malformed requests; run them before bridge implementation, then implement the bridge. Use direct `Atoms(...); atoms.calc=EMT(asap_cutoff=False)` as the reference, not the bridge's calculation helper.
- [x] Validate fresh evaluate requests, all identities, property shapes and finiteness before accumulation. Test timeouts/process errors and property omissions through actual subprocess fixtures; do not expose test-only production APIs.
- [x] Run focused Rust/Python boundary checks and commit the verified provider unit.

## Task 2: real runner, identity and restart

Files: modify `io/input.rs`, `config/validation.rs`, `runner.rs`, `io/restart.rs`, `force_model.rs`; add runner/restart cases to `tests/external_provider.rs`. Consume the System atomic-number field from Task 1.

- [x] Add a TOML acceptance test using `[model] type="ase_emt"`, explicit pins/Z/total composition, a 4-Cu FCC structure and real units. Assert that unsupported coupling/composition is rejected before process launch. Observe failure first.
- [x] Add `System.atomic_numbers`, backwards-compatible checkpoint defaults, and the external ModelSpec variant. Keep native physical signatures unchanged; derive the external signature from actual runtime provider identity plus explicit Z/PBC/parameters.
- [x] Branch force construction around the pair-specific path. Native pair tables and composite policy remain in their existing branch; an external total model receives no automatic pair/bonded/electrostatic terms. Declare provider-managed neighbors explicitly, with no invented cutoff.
- [x] Add provider-stress pressure metadata. For observed E/F/stress, compute W from `-V trace(stress)` and use existing thermodynamic unit rendering. Propagate failed-step error events for an evaluation error.
- [x] Verify skew-cell row/column mapping, finite-difference force and strain W with spec tolerances, translation invariance, selected-Z mutation and split restart identity/mechanics through real ASE. Run native regressions, then commit.

## Task 3: CLI, Python and MCP surface

Files: modify CLI main/tests, Python config/exports/tests/optional dependencies, MCP schema/tools/capabilities/tests.

- [x] Write CLI `provider-info --python PATH` and Python typed rendering tests; observe missing-command/type failures. Add discovery and typed `AseEmt` fields, preserving native rendering and path resolution.
- [x] Write an MCP JSON request test for external scope/pins and provider metadata. Describe optional runtime execution, cancellation latency and lack of pair-kernel cost estimates; reject unsupported cost estimates instead of inventing zero or classical estimates.
- [x] Expose identity/capabilities in final run metadata and live status consistently; do not duplicate native defaults into provider pins.
- [x] Execute real CLI/Python/MCP integration against the installed optional reference environment and commit.

## Task 4: verification and reproducible usage

Files: create `docs/B01_EXTERNAL_PROVIDER.md`, `docs/B01_VERIFICATION.md`; update validated scope, README/architecture/changelog and optional ASE test instructions.

- [x] Record exact installed versions/model hash, snapshot/finite-difference tolerances and executable commands. Distinguish direct-EMT parity from physical material validation and one-shot transport from B09 persistence.
- [x] Run `cargo fmt --all -- --check`, default/GPU all-targets Clippy with `-D warnings`, default/GPU Rust 1.90 checks and complete workspace tests.
- [x] Run complete Python/benchmark-tool tests, existing 67-observation smoke and optional real-ASE tests. Record skips as skips, not reference passes.
- [x] Controller: independent whole-branch review, test-first corrections and scoped closure review completed. All four final findings and related UTF-8 issue are closed at `05b5642`. B01 remains on its separate branch; B06 merge is not implied.

Final B01 status: implementation and local verification/review complete within
the documented ASE EMT scope. The discovered pressure factor, failed-instance
safety, cancellation and interpreter-resolution defects were corrected and
reviewed. [B01 verification](../../B01_VERIFICATION.md) records source-specific
evidence; [review closure](../../B01_REVIEW.md) records decisions and limitations.
B06 PR #53 was subsequently squash-merged with the owner's administrator
exception. B01 is prepared on `ricardo/b01-ase-emt-provider`; its rebase onto
merged main preserved the reviewed tree exactly. See the verification record's
shipping provenance. GitHub's account billing lock prevents remote CI from
starting, so local verification is not described as remote CI success.

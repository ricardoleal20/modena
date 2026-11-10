# B09 Task 2 report: native LJ calculator for ASE

## Scope and source epoch

Implemented the optional inverse adapter, `from modena.ase import
ModenaCalculator`, backed by the strict `modena snapshot` CLI. Base `import
modena` remains free of ASE/NumPy imports and mandatory dependencies.

Review base: `d5e6f0d8d47b1898039c5364e93af3b0b55d4fe2` (Task 1 production
through `2586f43784c82aa7cfbaba4c4f031ec0b0526037`, plus controller's plan-only
commit). No Task 1 implementation was changed.

Initial native/adapter source manifest recorded at 2026-11-05T16:16:36Z under
`target/b09-evidence/task2/provenance-initial.json`:

- Source epoch SHA-256 (sorted JSON mapping of nine task source/docs paths to
  their SHA-256): `fb7d4d983885480e0d04cbc1fa4f0fbae2eed8b0271ac1e0ad6f85ab9a29556f`.
- Actual debug CLI SHA-256:
  `5da2355cfcea23797d5cdf51f93dec79a8a62133df28f918b6fb40eeaa3ab6a9`.
- macOS 26.6.2 arm64; actual ASE environment Python 3.14.5, ASE 3.29.0,
  NumPy 2.5.4. Existing ASE/NumPy `a.shape` deprecation warnings were visible.

After the test-fixture and inherited harness-lock integration corrections below,
the final eleven-file manifest at 2026-11-05T16:23:05Z is in `provenance.json`:
`777835691054292c23c81f457f0713f8cccdd2535e1b76abb204886aa771067a`.
All native/adapter production source and the debug executable are identical
between these epochs. Only one existing MCP test fixture and the standalone
harness lockfile were added to the final source manifest.

The report itself is excluded from the source epoch. Evidence is ignored under
the shared target directory, not committed as generated artifacts. Commands
ran in `/Users/ricardoleal20/.codex/worktrees/modena-b09/modena` with
`CARGO_TARGET_DIR=/Users/ricardoleal20/Documents/Code/Misc/modena/target`.

## Files and behavior

- `crates/modena-core/src/snapshot.rs` / public `snapshot` module in `lib.rs`:
  versioned typed input/output, explicit Angstrom/eV units, homogeneous atomic
  numbers 1..118, finite coordinates/parameters, full PBC, positive diagonal
  orthorhombic cell, safe cutoff-plus-skin radius, resource/allocation guards.
  Uses existing `LennardJones`, `NeighborList`, `PairForceModel::compute` with
  one species index. No duplicate production LJ pair arithmetic. Singular
  minimum-image neighbors are rejected before the shared pair kernel can skip
  them. Returned energy/forces must be finite.
- `crates/modena-cli/src/main.rs`: strict JSON `snapshot` from stdin, bounded
  to 8 MiB, rejects unknown/missing/ill-shaped fields and trailing input,
  writes only one result JSON on success. No integration/output file creation.
- `python/modena/ase.py`: explicit optional import, energy/forces only,
  subprocess argument list without shell, one executable path/PATH name,
  positive finite timeout, exit/timeout/invalid-output diagnostics. Strict
  numeric protocol/unit/shape/finiteness checks before committing results.
  ASE atom changes and `set()` changes to model/binary/timeout invalidate
  properties. All accesses validate constraints/charge/spin, including cached
  accesses. Failures and unsupported property requests reset results.
- `crates/modena-core/tests/snapshot.rs`,
  `crates/modena-cli/tests/snapshot.rs`, `python/modena/tests/test_ase.py`:
  scientific oracle, finite differences, actual executable boundary and cache,
  malformed/unsupported input/output and failure behavior.
- `crates/modena-cli/README.md`, `python/README.md`: complete request and
  calculator examples, shifted-force convention, units, supported subset,
  unsupported stress, parameter changes, executable/deadline and resource caps.
- `crates/modena-mcp/tests/final_integration.rs`: controller-authorized scope
  exception correcting the expected schema key list to include Task 1's
  `transport`, after the full workspace run exposed the stale expectation.
- `bench/validation/harness/Cargo.lock`: controller-authorized inherited B01
  dependency integration repair for the required Rust 1.90 CI harness check.
  Adds the missing `tempfile` dependency and its absent transitive packages;
  every previously locked package/version/checksum is preserved.

## TDD and focused evidence

1. Wrote CLI/calculator boundary tests before implementation, then ran:

   `rtk proxy env CARGO_TARGET_DIR=/Users/ricardoleal20/Documents/Code/Misc/modena/target cargo test -p modena-cli --test snapshot --locked --offline -- --nocapture`

   RED: exit 101, 0 passed / 2 failed, both exposing the missing `snapshot`
   subcommand (`cli-red.log`).

   `rtk proxy /Users/ricardoleal20/Documents/Code/Misc/modena/target/b01-ase-venv/bin/python -m unittest python/modena/tests/test_ase.py -v`

   RED: exit 1, six explicit `missing optional ModenaCalculator` assertion
   failures, one base optional-import isolation pass (`python-red.log`).

2. Implemented snapshot and adapter. CLI GREEN: exit 0, 2/2; calculator
   initial GREEN: exit 0, 7/7 (actual dependencies, no skips), in
   `cli-green-attempt.log` / `python-green-attempt.log`.

3. Strengthened malformed-response types (numeric strings/booleans must not
   become accepted floats). Ran `... -m unittest python/modena/tests/test_ase.py
   -k transport -v`: RED exit 1, `CalculationFailed not raised` in
   `response-types-red.log`. Added strict JSON numeric type checks, then full
   focused calculator GREEN below.

4. Final focused commands:

   `rtk proxy env CARGO_TARGET_DIR=/Users/ricardoleal20/Documents/Code/Misc/modena/target cargo test -p modena-core -p modena-cli --test snapshot --locked --offline -- --nocapture`

   GREEN: core 2/2, CLI 2/2, no ignored tests (`native-focused-green.log`).

   `rtk proxy /Users/ricardoleal20/Documents/Code/Misc/modena/target/b01-ase-venv/bin/python -m unittest python/modena/tests/test_ase.py -v`

   GREEN: 7/7, no skips (`python-focused-green.log`). Two test fixture mistakes
   during refinement were corrected from independent arithmetic/input intent:
   a mistyped energy literal and a mixed-species mutation that temporarily
   changed a Cu number to Cu. These are recorded in `core-focused.log` and
   `python-green.log`; neither was used as evidence of production RED/GREEN.

For sigma=3.4 Angstrom, epsilon=0.0103 eV, r=4 Angstrom and rc=8.5
Angstrom, independent direct 12-6/derivative arithmetic yields
U_sf = -0.008978481029880884 eV, F0x = +0.005608613407781233 eV/Angstrom.
Tests use direct analytic formulas/literals (not production potential helpers),
Newton pair symmetry, zero transverse force, cutoff zero, periodic translation,
and centered energy finite differences. CLI tests run in empty temporary
directories and assert no files are created, on both successful/failed calls.

## Full validation

Commands below use the shared `CARGO_TARGET_DIR` above. Actual ASE checks set
`MODENA_ASE_PYTHON=/Users/ricardoleal20/Documents/Code/Misc/modena/target/b01-ase-venv/bin/python`
and `MODENA_BINARY=/Users/ricardoleal20/Documents/Code/Misc/modena/target/debug/modena`.
Exact argument arrays, original statuses and elapsed seconds are retained in
`checks.json`; scoped follow-ups in `scoped-checks.json` and
`python-worktree-checks.json`, with matching log files.

| Command | Result |
| --- | --- |
| `rtk proxy cargo test --workspace --locked --offline -- --nocapture` | Original complete progression: 349 passed / 1 failed / 0 ignored in 239.38 s; stopped at the stale MCP schema expectation, leaving the last MCP suite and docs unexecuted. All core/CLI suites passed and actual provider branches executed (no `SKIP` lines). `workspace.log` retained unchanged. |
| `rtk proxy cargo test -p modena-mcp --locked --offline -- --nocapture` | After the test-only fixture correction: 47 passed / 0 failed / 0 ignored, including the previously unexecuted final sweep suite; actual ASE branch executed. `mcp-after-fixture.log`. |
| `rtk proxy cargo test --workspace --doc --locked --offline -- --nocapture` | Exit 0; one preexisting intentionally ignored core builder example, zero failed docs. `workspace-docs.log`. |
| `rtk proxy /Users/ricardoleal20/Documents/Code/Misc/modena/target/b01-ase-venv/bin/python -m unittest discover -s python/modena/tests -v` | Exit 0; all 33 pass, zero skips, actual optional ASE/NumPy. `python-ase-worktree.log`. |
| `rtk proxy /Users/ricardoleal20/.pyenv/versions/3.14.5/bin/python -m unittest discover -s python/modena/tests -v` | With `MODENA_ASE_PYTHON` unset: exit 0, 33 collected, 26 pass and 7 explicit skips (6 calculator tests without ASE/NumPy, 1 actual-provider opt-in test). `python-base-worktree.log`. |
| `rtk proxy cargo fmt --all --check` | Exit 0 before and after fixture correction. |
| `rtk proxy cargo clippy --workspace --all-targets --locked --offline -- -D warnings` | Exit 0. |
| `rtk proxy cargo clippy --workspace --all-targets --locked --offline --features modena-core/gpu -- -D warnings` | Exit 0. |
| `rtk proxy cargo +1.90.0 check --workspace --all-targets --locked --offline` | Exit 0, actual compiler floor. |
| `rtk proxy cargo +1.90.0 check --workspace --all-targets --locked --offline --features modena-core/gpu` | Exit 0, actual compiler floor with GPU feature. |
| `rtk proxy cargo +1.90.0 check --locked --manifest-path bench/validation/harness/Cargo.toml --all-targets --offline` | Exit 0 after minimal nested lock synchronization. `msrv-harness-synchronized.log`. |

Combined non-overlapping workspace coverage is **351 passing tests plus one
intentionally ignored doctest**, using unchanged successful core/CLI coverage
and corrected complete MCP coverage. This is combined coverage, not a claim
that the original full workspace command exited zero. No unchanged heavy core
suite was repeated. Default GPU/MPI runtime suites are feature-gated to zero
tests; GPU compilation/lint passed, no MPI or GPU hardware run is claimed.

The first Python attempts are retained as failed environment checks:
`python-ase.log` had two existing workflow tests denied temporary writes in the
managed worktree; `python-base.log` used system Python 3.9 (unsupported, missing
`tomllib`). Corrected Python 3.14.5 runs from writable `/private/tmp` passed
first (`python-ase-writable.log`, `python-base-supported.log`), and the exact
worktree commands above were then rerun with permission to create those
preexisting temporary cwd fixtures. No source changes were needed.

The harness check initially failed because its nested lock lacked `tempfile`,
added to core in B01 (`59ef2c6`). `.github/workflows/ci.yml:34` requires this
locked Rust 1.90 check. The identical failure was reproduced from an archived,
unchanged `d5e6f0d` baseline in a temporary directory
(`msrv-harness-baseline.log`), so it predates Task 2. After controller
authorization, `rtk proxy cargo +1.90.0 check --manifest-path
bench/validation/harness/Cargo.toml --all-targets --offline` synchronized the
lock. Diff inspection and `harness-lock-delta.json` confirm no removals or old
pin/version/checksum changes: only `modena-core` gains `tempfile`, plus nine
absent packages (`tempfile`, `fastrand`, `once_cell`, `rustix`, `bitflags`,
`errno`, `linux-raw-sys`, `windows-sys`, `windows-link`). The final locked check
then passed. No manifest or unrelated dependency upgrade was made.

## Self-review and limits

Reviewed the implementation against the Task 2 brief/spec: no ASE imports from
base package; no mass/velocity/integration or thermodynamic unit conversion;
no native model generalization; no scalar-virial-to-tensor inference. Stress
and any other unimplemented property explicitly fail. Unsupported constraints,
charges/spin, species/cells/PBC are rejected; constraint changes cannot bypass
validation via ASE cache hits. Subprocess failures and malformed responses
cannot retain old properties. The strict CLI returns no partial stdout on
error and bounds input before JSON parsing.

This is one process per snapshot and native single-species shifted-force LJ
only. No stress, mixed-species model, parameter database, rotated/skew cell,
partial PBC, constrained or charged/spin calculation is provided. Resource caps
are 8 MiB JSON, 100000 atoms, and 1000000 dense neighbor cells; they do not
constitute a large/dense system throughput claim. Near-singular neighbor pairs
with r² < 1e-14 Angstrom² are rejected. MPI/runtime GPU scientific validation
is outside this task; default/GPU compile/lint coverage is listed above.
Task 3 owns release measurement, full-run persistent transport evidence and
whole-branch review. No push or merge was performed in this task.

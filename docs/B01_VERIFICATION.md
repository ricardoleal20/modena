# B01 local verification

Executed on 2026-11-05 from `/private/tmp/modena-b01`, branch
`feat/b01-ase-emt-provider`, source head
`4214832a010281c551a31b888012a04dd0191111`. The worktree additionally contained
the documentation/plan changes and a comment-only MCP test clarification in
this documentation commit; production source was unchanged. Local checks below
passed. The four findings from subsequent whole-branch review are addressed in
one final fix wave, including independent SI pressure evidence below. Scoped
rereview of the combined corrections remains pending.

B06 remains open as PR #53: GitHub Actions cannot start because of the account
billing lock, and the required administrator-exception answer is pending. These
local results do not establish remote CI success or a merge. [B06's record](B06_VERIFICATION.md)
retains its historical revision and measurements.

## Environment and installed identity

macOS 26.6.2 (25G83), arm64; rustc 1.97.1
(`8bab26f4f`, 2026-07-14), cargo 1.97.1 (`c980f4866`, 2026-06-30).
Compiler-floor checks used rustc 1.90.0 (`1159e78c4`, 2025-09-14).
The isolated interpreter was Python 3.14.5 at
`/Users/ricardoleal20/Documents/Code/Misc/modena/target/b01-ase-venv/bin/python`.
Release `provider-info` actually discovered:

| Field | Value |
| --- | --- |
| Provider / model / protocol | `ase_emt` / `EMT` / `modena-ase-emt-v1` |
| Task / composition | `energy_forces_stress` / `total` |
| Installed ASE / NumPy | `3.29.0` / `2.5.4` |
| EMT source/model-data SHA-256 | `8ce8bb2f90e8afad8f091b8c3f600ea260534e5c270a6b2c82cae1817cdf3881` |
| ASAP cutoff | `false` |
| Wire units | Å, eV, eV/Å, eV/Å³ |

This is an installed EMT source hash, not learned weights. See the
[usage/migration guide](B01_EXTERNAL_PROVIDER.md) for fields and limits.

## Fresh automated checks

All commands ran through `rtk proxy`, from this worktree, with
`CARGO_TARGET_DIR=/Users/ricardoleal20/Documents/Code/Misc/modena/target`.
Cargo dependency/build commands used `--locked --offline`. Log names below are
relative to the local evidence directory described at the end.

| Command after `rtk proxy` | Result / log |
| --- | --- |
| `cargo fmt --all -- --check` | Exit 0; `fmt.log` |
| `cargo clippy --locked --offline --workspace --all-targets -- -D warnings` | Exit 0; `clippy-default.log` |
| Same Clippy command with `--features modena-core/gpu` before `--` | Exit 0; `clippy-gpu.log` |
| `cargo +1.90.0 check --locked --offline --workspace --all-targets` | Exit 0; `msrv-default.log` |
| Same check with `--features modena-core/gpu` | Exit 0; `msrv-gpu.log` |
| `cargo test --locked --offline --workspace -- --nocapture` | 320 reported passes, 1 ignored, 39 suites, exit 0 in 210.034 s; `workspace-default.log` |
| `cargo package --list --locked --offline -p modena-core --allow-dirty` | Exit 0; includes `src/external/ase_emt_worker.py`; `package-list.log` |
| Python wrapper suite, no ASE opt-in | 25 tests: 24 executed successfully, 1 explicitly skipped; `python-default.log` |
| Python wrapper suite, ASE opt-in | 25 passed, no skips; `python-ase.log` |
| `cargo test --locked --offline -p modena-core --test external_provider -- --nocapture`, ASE opt-in | 18 passed, no skips; `core-ase.log` |
| `cargo test --locked --offline -p modena-cli -p modena-mcp -- --nocapture`, ASE opt-in | 50 passed, no skips; `interfaces-ase.log` |
| Direct worker Python reference suite | 3 passed; `worker-ase.log` |
| Benchmark-tool Python suite | 22 passed outside the macOS timer sandbox restriction; `bench-tools-unsandboxed.log` |
| Final `cargo fmt --all -- --check` / `git diff --check` after documentation edits | Both exit 0; `final-fmt.log`, `final-diff-check.log` |

The default workspace's 320 successful Rust functions include five optional
checks that printed `SKIP` and returned: three core, one CLI and one MCP. They
did not execute ASE. The single ignored item is the existing illustrative
builder doc block. Optional-reference commands above separately executed ASE.
GPU feature checks compile/lint only; no GPU hardware or native MPI runtime
campaign was performed.

The benchmark suite initially failed four unchanged tests under the sandbox:
`test_failure_retains_evidence_and_fails`,
`test_real_release_cli_uses_valid_no_file_deck_and_reports_drift`,
`test_success_has_actual_steps_numeric_renderer_and_provenance` and
`test_missing_executable_retains_nullable_hash_and_failure`.
`/usr/bin/time -l /usr/bin/true` reproduced
`sysctl kern.clockrate: Operation not permitted`. The authorized local rerun
outside that restriction passed all 22, without weakening scripts/assertions.
`bench-tools.log` retains the initial failures. This tests benchmark tooling,
not performance. The direct independent ASE reference retained its upstream
`ase/atoms.py:432` NumPy-2.5 shape-assignment DeprecationWarning; no broad
warning filter or dependency change was introduced.

To reproduce the Python and optional checks, let `ASE_PYTHON` be the isolated
interpreter and `TARGET` the Cargo target directory, then use:

```bash
rtk proxy env -u MODENA_ASE_PYTHON MODENA_BINARY="$TARGET/release/modena" \
  "$ASE_PYTHON" -m unittest discover -s python/modena/tests -v
rtk proxy env MODENA_ASE_PYTHON="$ASE_PYTHON" MODENA_BINARY="$TARGET/release/modena" \
  "$ASE_PYTHON" -m unittest discover -s python/modena/tests -v
rtk proxy env CARGO_TARGET_DIR="$TARGET" MODENA_ASE_PYTHON="$ASE_PYTHON" \
  cargo test --locked --offline -p modena-core --test external_provider -- --nocapture
rtk proxy env CARGO_TARGET_DIR="$TARGET" MODENA_ASE_PYTHON="$ASE_PYTHON" \
  cargo test --locked --offline -p modena-cli -p modena-mcp -- --nocapture
rtk proxy "$ASE_PYTHON" crates/modena-core/tests/test_ase_emt_worker.py
rtk proxy "$ASE_PYTHON" -m unittest discover -s bench/tests -v
```

The benchmark-tool real-CLI test expects `target/release/modena` in the checkout.
For this worktree, ignored `target` was linked to the shared target directory;
`CARGO_TARGET_DIR` alone does not change that test's hard-coded path.

## What the reference checks establish

Independent direct `Atoms(...); atoms.calc=EMT(asap_cutoff=False)` calculations
are compared with the public provider on displaced four-Cu periodic snapshots,
both orthorhombic and skew. Cell vectors are rows. Fixed absolute bounds:

| Quantity | Tolerance / perturbation |
| --- | --- |
| Snapshot energy | `1e-10 eV` |
| All force components | `1e-10 eV/Å`, relative tolerance zero |
| Python stress components | `1e-12 eV/Å³`, relative tolerance zero |
| Rust scalar virial `-V trace(stress)` | `1e-9 eV` |
| Selected force central difference | `h=1e-5 Å`, error at most `2e-6 eV/Å` |
| Fixed-fractional isotropic-strain virial | `h=1e-6`, error at most `2e-6 eV` |
| Common translation / net force | Energy/force invariance and force-sum bounds `1e-10` in respective wire units |

Actual runner tests also cover selected Z mutation, split continuation versus
uninterrupted mechanics, changed pins/runtime identity, native signature
compatibility, known stress-derived virial/pressure plumbing using the current
unit factor, and failed-step events. The final-fix independent SI/ASE pressure
regression below additionally validates the physical atm conversion.
Subprocess fixtures exercise malformed/missing/nonfinite properties, conversion
overflow, timeout, bounded diagnostics and accumulator atomicity without ASE.
This is direct-model parity and bounded integration evidence; it is not
material validation, a converged ensemble study or independent-engine parity.

The actual Python test executes discovery → rendering → check → run and exact
metadata equality. The actual MCP test executes newline-delimited stdio
JSON-RPC `run_probe`, in-process JSON probe and async lifecycle/final status.
The status test accepts populated metadata even if already completed; it does
not prove a guaranteed nonterminal live observation. Its comment now states
that scope. External estimates explicitly error without launching Python.

## Rebuilt binaries, CLI runs and scientific smoke

The shared cache initially contained the old B06 scientific binary with hash
`d6ef67f94c0c166251ed7d3bfb8eb286af3680041ed5477cd8398dbdfbaa4373`
and old CLI hash `cbc18e44b6e14ec3191fd9e9efe9b246e9ebeb1f91e8e0c22aabad9cadc5050d`.
Before attributing new results, both were rebuilt from this B01 worktree:

```bash
rtk proxy env CARGO_TARGET_DIR="$TARGET" cargo build --locked --offline --release \
  -p modena-cli -p modena-core --example scientific_validation --bin modena
rtk proxy "$TARGET/release/modena" provider-info --python "$ASE_PYTHON"
rtk proxy "$ASE_PYTHON" bench/validate.py \
  --binary "$TARGET/release/examples/scientific_validation" \
  --output "$TARGET/b01-scientific-smoke"
```

Build exit 0; `release-build.log` records compilation from `/private/tmp/modena-b01`.
Smoke at `2026-11-05T08:07:42.009860+00:00`, one Rayon thread: **67/67
existing scalar observations passed**, schema 2, six analyzed series and two
groups of three replicas. Expectations/tolerances were unchanged. Its raw
stdout hash equals the earlier deterministic B06 output, but the executable
hash and captured worktree provenance are new. Those existing observations
cover native numerical fixtures, not external EMT material validation.

The ordinary CLI example used discovered exact pins, Z=[29,29,29,29], explicit
Cu mass 63.546, full PBC, real units, two Velocity-Verlet steps at 0.0001 ps,
10 K initial temperature and seed 42. Discovery/check/run all exited 0;
metadata's provider descriptor exactly equals discovery, pressure is `unconstrained_provider_stress`,
and outputs include three CSV rows, three diagnostic samples, trajectory and
final checkpoint. The short Langevin NVT variant used the same pins/mass/cell,
`dt=0.0001`, `gamma=1.0`, target 10 K and seed 42. It also exited 0 and verified
`ensemble=nvt`, `kinetic_dof=12`, `output_units=real`, the same provider,
`verdict=not_applicable` and three full-rate samples. These are integration
checks only. The [usage guide](B01_EXTERNAL_PROVIDER.md#reproduce-a-short-cu-workflow)
reproduces the NVE input; substitute `langevin_baoab`, gamma/target/seed to
exercise NVT.

| Artifact | SHA-256 |
| --- | --- |
| Rebuilt release CLI | `8639f2beae459c8ca5ef7225438d7ac7feacc6a5bb65afc7632a37aae7bb5d1d` |
| Rebuilt scientific executable | `e35c176e655d93a93467d812d601c6df5ccf499046ceb3406749d6a87f4d921c` |
| Smoke raw stdout (1,534,589 bytes) | `394797f8789becbbcad3fa556f29154851875bf8be25951e4d779dd87ad8dac8` |
| Smoke validation JSON | `97fe8b6e0874159a83f40098001d22e4d46c9bf494145196e328edae5c4dc270` |
| NVE `cu.toml` | `6d58b962fa20b19a22916a5cd6bacc6ab100ecf081938ccdb2b15eb302ed19cf` |
| NVE `run.stdout` | `4fc5c93a83ebc5f17d0fc8bd1db8e3730c9dfc8a1f7c9b2c83ebf0382e539660` |
| NVE `cu.checkpoint` | `6f65ad58c396337dce04c19e6f1379be6f1a479fa7ad93fb8ea158d4bde29787` |
| NVE `metadata.json` | `f1bc648bc460e0a11b4b6195b27f359bca43bb93fa0dc091a662a7028124dd19` |
| NVT `cu.toml` | `56c03628a3e07aa00d9acabb71af0948d6d5106fa7f3525a8a4f7fcadef2d588` |
| NVT `run.stdout` | `6f1897dde96d90df6c92ea63ff3c1edcbeecd493b6fa09af118eec7d8fe50742` |
| NVT `cu.checkpoint` | `189ee1631acee2f215ac2446e04c18fa62de1b1c75c1b711877a66c3509edb74` |
| NVT `metadata.json` | `53f8e0e6aaace4246f1c379cb596cf6b3f7689cced35f4afe635d40abb2c15eb` |

## Corrected pressure units and final-fix evidence

The final fix corrects the former 69477.8 scale once in shared units. From
thermochemical kcal (4184 J), Angstrom, exact Avogadro constant and standard
atmosphere (101325 Pa), the derived factor is
`4184/(6.02214076e23*1e-30*101325) = 68568.42296625092` atm per kcal/(mol Å³).
Definitions: [NIST SP811 Appendix B.8](https://www.nist.gov/pml/special-publication-811/nist-guide-si-appendix-b-conversion-factors/nist-guide-si-appendix-b8)
and [BIPM mole](https://www.bipm.org/en/si-base-units/mole), checked 2026-11-05.
This factor is calculated from those definitions; it changes real output,
pressure-target normalization and inverse-pressure compressibility. Reduced
configurations and real fixed-cell mechanics remain unchanged.

The immutable pre-fix cold-Cu evidence stays at
`b01-evidence/controller-pressure-probe/`: old CLI -11215.5992549632 atm versus
direct ASE/SI -11068.801164315795 atm. Corrected CLI regression gives
-11068.8011643184 atm (CSV rounding), within 1e-7 atm of the independent
reference, including an interpreter relative to the input directory and a
preserved venv symlink. Unit/non-unit conversion tests use independent literal
SI factors in both directions; MC/Berendsen target and compressibility tests do
not reuse the production conversion helper for expected values.

Fresh final-fix logs, command records, provenance, corrected cold-Cu files and
hashes live separately under `b01-evidence/final-fix/`. New scientific smoke
artifacts live under `b01-scientific-smoke-final-fix/`. Earlier Task 4 logs,
manifests and binary hashes above describe their historical source epoch and
are preserved; rebuilding does not relabel those old artifacts.

Regression coverage additionally checks exact fixed-cell pressure-only v2
migration, stochastic/extended-state continuation, unrelated-factor/provider/Z
rejections, old pressure-coupled rejection with explicit state-only acceptance,
terminal failed instances that preserve good checkpoints, synchronized
setup/first-step cancellation and CLI/MCP relative interpreter handling.
The [guide](B01_EXTERNAL_PROVIDER.md#checkpoints-and-public-api-migration)
documents the narrow compatibility exception and public Cancelled/constructor
API. A failed instance does not promise buffer rollback. Cancellation uses safe
boundaries and one in-flight provider deadline plus bookkeeping, without a hard
global wall-time guarantee. Known ASE/NumPy warnings remain visible.

Final scoped rereview and B06 merge remain controller-owned and pending.

Final fix checks executed from base `81414606560881ae978738e330a3ba1ee4114e4b`
plus the recorded working patch (commit range is recorded in the local fix report):

| Fresh final check | Result |
| --- | --- |
| Native workspace | 330 passed, 0 failed, 1 existing ignored doc example; 6 explicit optional ASE skips |
| Actual ASE provider | 24 passed, no optional skips |
| Actual CLI/MCP interfaces | 53 passed, no optional skips |
| Python wrapper | 25 actual-ASE tests passed; native-only run has 24 executed + 1 explicit skip |
| Python worker | 3 passed; known ASE/NumPy deprecation warning visible |
| Benchmark tooling | 22 passed against rebuilt release CLI, authorized macOS timer access |
| Scientific smoke | 67/67 observations, schema 2; no tolerances relaxed |
| Format, default/GPU all-target clippy, Rust 1.90 default/GPU checks | All passed |

Rebuilt release CLI SHA-256:
`f1e374f87c8e6f264ad6525d9344612ddf7660d055f41c904b159d4ee1a6dea1`.
Rebuilt scientific executable SHA-256:
`26879a5e33913f20b88ce3a06608b33988b2d118d0b2a876e546193cfb2dbca5`.
The retained corrected cold-Cu comparison has absolute error
`2.6047928258776665e-9 atm`. Full logs and per-artifact hashes are under the
separate final-fix evidence directory; earlier manifests above are unchanged.

## Local evidence locations

These are **untracked local build outputs**, not files supplied by a fresh
checkout. Shared target root is
`/Users/ricardoleal20/Documents/Code/Misc/modena/target`:

- `b01-evidence/`: full named logs, `commands.jsonl`, pre/post-build provenance,
  verification driver and summary. [Local workspace log](../target/b01-evidence/workspace-default.log).
- `b01-evidence/cu-cli/` and `cu-cli-nvt/`: complete input, stdout/stderr,
  provider/metadata/diagnostic JSON, CSV, trajectory, checkpoint and per-file
  `hashes.json`. [Local NVE manifest](../target/b01-evidence/cu-cli/hashes.json),
  [local NVT manifest](../target/b01-evidence/cu-cli-nvt/hashes.json).
- `b01-scientific-smoke/`: raw `validation.stdout`, `validation.stderr` and
  structured `validation.json` with command, source head, dirty status,
  compiler/interpreter/binary hash and statistical-analysis provenance.
  [Local smoke JSON](../target/b01-scientific-smoke/validation.json).

The core package listing confirms the Python worker is bundled inside
`modena-core`; installation does not depend on a worker outside that package.
No publish, push, merge, billing or protection change was performed.

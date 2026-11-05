# B02 native CPU pair scratch verification

The fresh local campaign passed **357 Rust tests**, **35 Python tests with
actual ASE**, **five worker tests**, **34 benchmark-tool tests**, and all
**67 scientific observations**. One existing illustrative builder doctest was
ignored. No optional ASE Rust checks were skipped. GPU verification was
compilation/lint only; MPI runtime and other-engine performance were not tested.

## Source epoch and environment

Executed from the clean managed worktree
`/Users/ricardoleal20/.codex/worktrees/modena-b02/modena`, source
`360defbdb2e64b61d37aabcd36c76fb0d19ea100`, on 2026-11-05,
**18:32:58–18:38:25 UTC**. This includes the accepted Task 1 review fixes at
`bf2ef75`; subsequent changes in this verification commit are documentation
only. No force kernel, assertion, tolerance or retained timing row changed.

Local evidence lives in the ignored shared target directory:
`target/b02-evidence/task2-360defb/`. `provenance.json` records the clean source
epoch, hashes of 728 tracked files, compiler/interpreter versions and exact
environment. `commands.json` records each command, working directory, source,
start/end time, status, duration, environment and full-log SHA-256. The retained
`validation-driver.py` reproduces this particular campaign; use a fresh output
directory on another run. These artifacts and executable copies are local
outputs, unavailable in a fresh clone.

Host: Apple M3 Pro, 11 logical CPUs, macOS 26.6.2 arm64. Compiler: rustc 1.97.1
(`8bab26f4f`, LLVM 22.1.6); compiler-floor checks actually used rustc 1.90.0
(`1159e78c4`). Actual-reference Python was 3.14.5 at
`/Users/ricardoleal20/Documents/Code/Misc/modena/target/b01-ase-venv/bin/python`,
with ASE 3.29.0 and NumPy 2.5.4. The supported no-ASE interpreter was
`/Users/ricardoleal20/.pyenv/versions/3.14.5/bin/python`; discovery confirmed ASE
absent. The venv executable symlink was preserved.

`CARGO_TARGET_DIR` was
`/Users/ricardoleal20/Documents/Code/Misc/modena/target`; the worktree's ignored
`target` link points there. `MODENA_ASE_PYTHON` selected the exact venv and
`MODENA_BINARY` selected a debug CLI freshly built from this worktree before
tests. `OMP_NUM_THREADS=OPENBLAS_NUM_THREADS=1` limited optional-library resource
pressure. `RAYON_NUM_THREADS`, `RUSTFLAGS` and `CARGO_ENCODED_RUSTFLAGS` were
unset. Native tests retained normal Rayon pool selection; only the scientific
driver explicitly set one Rayon thread. macOS timer queries and temporary
fixture working directories ran with host access. The initial load averages
were 3.8213/3.7607/4.0991 on this shared host; this was no new timing campaign.

## Fresh checks

Commands below are shown after `rtk proxy`; Cargo commands used the shared
target and `--locked --offline`. The full workspace test ran once on the
stabilized source. An initial debug build also succeeded before the driver's
logged build; it did not execute tests.

| Command / check | Result | Log |
| --- | --- | --- |
| `cargo build --locked --offline -p modena-cli --bin modena` | Exit 0; fresh debug CLI preserved | `debug-build.log` |
| `cargo test --workspace --locked --offline -- --nocapture` | 357 passed, 0 failed, 1 ignored across 43 reported suites; exit 0, 268.696 s | `workspace-actual-ase.log` |
| Actual Python: `python -m unittest discover -s python/modena/tests -v` | 35 passed, no skips | `python-actual-ase.log` |
| Base Python, `MODENA_ASE_PYTHON` unset, same discovery command | 26 executed successfully; nine explicit skips out of 35 | `python-base-no-ase.log` |
| `python crates/modena-core/tests/test_ase_emt_worker.py` | Five passed | `worker-actual-ase.log` |
| `cargo fmt --all -- --check` | Exit 0 | `fmt.log` |
| `cargo clippy --workspace --all-targets --locked --offline -- -D warnings` | Exit 0 | `clippy-default.log` |
| Same Clippy command, `--features modena-core/gpu` before `--` | Exit 0; GPU compile/lint only | `clippy-gpu.log` |
| `cargo +1.90.0 check --workspace --all-targets --locked --offline` | Exit 0 | `msrv-default.log` |
| Same Rust 1.90 check with `--features modena-core/gpu` | Exit 0 | `msrv-gpu.log` |
| `cargo +1.90.0 check --locked --offline --manifest-path bench/validation/harness/Cargo.toml --all-targets` | Exit 0 | `msrv-harness.log` |
| `cargo build --locked --offline --release -p modena-cli -p modena-core --bin modena --example scientific_validation` | Exit 0; executable copies preserved | `release-build.log` |
| `python -m unittest discover -s bench/tests -v` | 34 passed, no skips, against fresh release CLI | `bench-tools.log` |
| `python bench/validate.py --binary "$EVIDENCE/scientific_validation" --output "$EVIDENCE/scientific-smoke"` | 67/67 observations passed | `scientific-smoke.log` |
| `git diff --check` | Exit 0; repeated after documentation edits | `diff-check.log`, `final-diff-check.log` |

The no-ASE run skips eight calculator tests because optional ASE/NumPy are
absent, plus the explicit actual-provider opt-in check. All nine execute in
the actual-ASE run. The only ignored Rust item is
`config::builder`'s illustrative doctest. The upstream ASE/NumPy shape-assignment
deprecation warning remains visible in worker output; existing proptest
source-persistence notices are retained. No broad warning filter or assertion
relaxation was added.

The full workspace includes cutoff/C1 consistency, timestep-dependent NVE,
bounded zero-skin snapshots, actual ASE runner/direct reference/restart,
MCP/CLI, and scratch allocation/concurrency checks. The isolated scratch test
prints zero force-sized allocations for T1/4/8/11 and atom/pool transitions,
including N5000=195000 and N10000=390000 half pairs. These are fresh correctness
checks, separate from Task 1's retained allocation profiles and timings.

## Rebuilt binaries and scientific smoke

Both release artifacts were built from this worktree before use and copied
out of shared Cargo output. `binary-hashes.json` records:

| Artifact | SHA-256 |
| --- | --- |
| Fresh debug CLI | `cc63f7796cb78bc41c4f719a649398d5348be4d61fab182ec581bfe9371313d2` |
| Fresh release CLI | `1b39924b2ab9141748c49a1936df9a57cda8a70d3b879c5008da33996e0113f8` |
| Fresh release scientific example | `b7a9e1882cd0449f2624a00a6a0569c841d264266a68c32a7fc70e6078d83ceb` |
| Smoke raw stdout | `394797f8789becbbcad3fa556f29154851875bf8be25951e4d779dd87ad8dac8` |
| Smoke structured JSON | `05176371bd2e4aedc8653dc61b025d9daf34700f52b496f71fe8e510363775f0` |

The fresh release CLI is byte-identical to the immutable accepted
`target/b02-evidence/block/modena` built from measured production `0cea614`.
This supports retaining the historical timing evidence at its original epoch;
it does not relabel those measurements as newly executed. The new scientific
binary and smoke have their own source/hash provenance.

At 18:38:25 UTC, the schema-2 smoke passed the unchanged 67 observation gates.
Its NVE energy bands were 4.53571794478691e-5, 1.1338612830091854e-5 and
2.8347420262253564e-6 at dt 0.004/0.002/0.001; both halving ratios passed
4±0.3 (4.000240605049538 and 3.999874671202429). Six complete raw series of
13,667 values each remain in `scientific-smoke/validation.stdout`.
`validation.json` retains the six analyses and two three-seed replica groups,
with recorded burn-in and finite-series uncertainty. These observations and
estimates do not certify converged ensembles or material accuracy.

## Scope and retained measurements

[B02's guide](B02_PAIR_SCRATCH.md) retains the Task 1 source epochs, all 380 CLI
runs (314 measured, 66 warmups), one-thread file parity, allocation/static
profiles and rejected dense variant. The repeated N5000/four-thread runtime
benefit supports that declared case. Large/default-pool results remain
uncertain, including the negative N10000/11T series; durations, paired ratios
and ratios of medians remain distinct. Zero warmed force allocations coexist
with retained O(TN) scratch and sometimes higher process RSS. Instrumented
Rust-requested peak bytes and OS RSS are different measurements.

Task 2 added no performance samples and did not alter B09 or earlier B02
evidence. This record describes local verification only. Task/whole-branch
review and PR/CI state are reported separately by the controller; it makes no
claim of current CI success or final branch approval.

## Final saved-continuation guard follow-up

The final whole-branch review found that the timing validator checked the
declared physical signature but omitted actual saved identity, DOF and
controller state. Fix `940c95b` updates only the benchmark driver and tests;
it changes no Rust production source or force arithmetic. The earlier broad
357/35/five-worker/34-tool/67-observation results above remain attached to
their original `360defb` epoch and were not rerun for this tooling-only fix.

The validator now independently checks N saved masses of 1, integer species
0, integer atomic numbers 0, charges 0 and element labels `Ar`; zero
constraints; three conserved COM degrees of freedom; and explicit null saved
integrator, thermostat and barostat state for the fixed velocity-Verlet NVE
deck. Native classical atomic numbers are **unspecified zeros**, not 18
inferred from `Ar` labels. The general Rust v2 reader defaults an older
missing atomic-number field to an empty buffer; this fixed current benchmark
instead requires the explicit N-zero buffer written by its pinned binaries.
The other omitted identity/DOF/control fields are required by the Rust reader.

The success fixture now includes every continuation field and the normalized
LJ signature. Before the fix, both agreeing child processes with saved mass
1→2 or a missing `conserved_com_dof` still passed; the new regressions failed
as expected. Afterward, **eight focused driver tests** passed in 6.861 s and
the **full 37-test benchmark-tool suite** passed in 15.519 s, without skips.
This includes missing fields, malformed buffers, boolean/integer confusion,
wrong native Z, incorrect DOF/control state, and all previous subprocess,
portability and exact-file gates. Formatting and diff checks passed.

An additional focused probe used the actual immutable Rust CLI `check`
command: it loaded the complete v2 fixture through normal continuation,
rejected the conflicting saved mass, and rejected the missing COM-DOF field
at v2 parsing. The probe integrated no steps and produced no performance
sample. The fixture is synthetic test data, not a scientific trajectory.

Evidence is retained separately in `target/b02-evidence/final-fix/`, including
RED/GREEN logs, `rust-reader-proof.json`,
`retrospective-continuation-audit.json` and `production-source-proof.json`.
The new validator hash is
`1f7967cac640b1f6de1b636263964138f5c6dfd51ee353b5f36745c58b86fcfc`.
It validated all **380** original complete files and all **48** one-thread
checkpoint/thermo hash pairs; all 1,943 original evidence-file hashes remain
unchanged. Original loaded-driver hashes/source epochs and every timing row
remain separate and intact. No incompatibility required weakening a gate.

All 48 core/CLI production source/manifest hashes match the whole-branch
review base `8d20fed`; this is a comparison to that base, not a new claim that
the earlier documentation-corrected file matches measured `0cea614` source.
The preserved and shared release CLI bytes still have SHA-256
`1b39924b2ab9141748c49a1936df9a57cda8a70d3b879c5008da33996e0113f8`.
No rebuild, kernel change, timing campaign or broad-validation rerun was
needed. [Review closure](B02_REVIEW.md) records the approved scoped rereview,
PR publication and the actual remote CI limitations separately.

# B03 verification: rejected scalar packed-layout experiment

The measured serial-scalar SoA, AoSoA4 and Cluster4 prototypes fail the frozen
promotion criterion in every declared case. Keep production unchanged. Both
protocols retain all 96 comparisons, including unfavorable results; no SIMD or
parallel packed implementation was measured. [The complete experiment guide](B03_PACKED_LAYOUT.md)
publishes every case's full-step medians, paired directions, setup, RSS and
retained capacity, plus separate production-auto controls.

All three tasks and the final whole-branch review approved the experiment,
including two focused Task 2 fix rounds. [PR #57](https://github.com/ricardoleal20/modena/pull/57)
is open, assigned and labeled against B02. [Review closure and decisions](B03_REVIEW.md)
record the approved source, deferred Minor and actual remote CI limitation.
No merge was performed. The other seven research cases retain their own
acceptance gates; B04 is already active in a separate worktree.

## Fixed criteria and observed outcome

The [binding design](superpowers/specs/2026-11-06-b03-packed-layout-design.md)
fixes these gates before measurement:

| Gate | Evidence and outcome |
| --- | --- |
| Logical identity, permutation and masks | Reviewed Task 1 tests compare exact logical-ID roundtrips and actual consumed `(i,j,scale_bits)` masks against independently constructed half-pair sets. Padding, exclusions/scales, additive forces, coordinate/box changes, reorder/resize and same-epoch list replacement pass. Production System/topology/restart representation is unchanged. |
| Independent forces/energy/virial | Component force abs 1e-10 + rel 1e-12; energy/virial abs 1e-8 + rel 1e-12. Task 1 Rust scalar-half-pair and independent Python all-pairs checks pass. Task 2 direct-oracle-v2 audits all 1,344 primary invocations, including all 2,688 force-ready/final states. Maximum absolute F/U/W discrepancies are 1.6361911825413245e-14 / 1.6952981241047382e-9 / 5.791662260890007e-9. These are historical audited results, not Task 3 physics reruns. |
| Full-step and timestep convergence | Reviewed Task 1 tests compare complete Verlet trajectories with serial reference and retain B06's fixed physical duration 2, dt 0.004/0.002/0.001, halving-ratio criterion 4 +/- 0.3. Original B06 ratios 4.000240605049538 / 3.999874671202429; supported-cutoff experiment ratios 4.000119199096041 / 3.9999937828743466 for all modes. Task 3's scoped benchmark regression bundle repeats the existing Rust test and reports the same values; it does not replace or relabel Task 1's independent Python audit. |
| Performance in each declared case | Fraction saved = `1 - median(candidate full-step seconds) / median(its matched serial baseline seconds)`. Require >=5%, favorable direction in >=4/5 pairs, then a separate five-pair reproduction. All 96 primary and 96 snapshot-free cases fail. Zero reproduction campaigns trigger; zero candidates qualify for promotion. No pooled median across cases/protocols supplies acceptance. |
| Memory and regressions | Retained packed capacity and actual process high-water RSS are assessed per case. Cluster4 has large memory and full-step regressions; no arbitrary new memory threshold was introduced. See the guide's [resource assessment](B03_PACKED_LAYOUT.md#memory-and-regression-assessment). |

The complete case space is N 864/2048/5000/10000 × requested pools 1/4/8/0 ×
no-rebuild/forced × SoA/AoSoA4/Cluster4. Actual pools are 1/4/8/11. Each
candidate has its own interleaved-serial A/B warmup and five alternating AB/BA
pairs. Separate production-auto runs never enter scalar acceptance. Per
protocol, 192 scalar warmup + 960 measured + 32 production warmup + 160
production control invocations total 1,344, all valid and exit 0. The two protocols
total 2,688 invocations; the separate exploratory pilot retains 32. No adverse row
was removed. Steps(no-rebuild/forced) were frozen after pilot at 864:700/190,
2048:200/40, 5000:90/20, 10000:45/20, dt 1e-6, seed 17.

Primary best fraction saved is -1.5415778% (SoA, N 2048, requested 4, forced;
0/5 favorable). Snapshot-free best is -0.9116640% (AoSoA4, N 5000, requested 1,
forced; 0/5 favorable). Worst Cluster4, N 10000, default 11, forced costs 3.58294x
primary / 3.53042x snapshot-free. Each outcome is a rejection in its own case,
not an inference about another layout implementation or material.

## Source and executable epochs

Task 3 checks ran from the explicit B03 worktree
`/Users/ricardoleal20/.codex/worktrees/modena-b02/modena` at
`885a2ac26a3ac8a5f40f51a11dcd5095529ce5da`, with clean tracked source before
and after checks. Documentation changes in this closure commit do not change
that tested implementation epoch. B04 source was not used. No Rust source has
changed since reviewed Task 1 `66af0e4bcca329dd3445d1d3b202dd3055c6afc8`:
`git diff 66af0e4 -- Cargo.toml Cargo.lock crates` is empty. Production core tree
is `39b46c7496c6820219056fe38a821378a32633e8`, identical to B02 at `604039e`;
the entire `crates/modena-core` projection is unchanged. This verifies the
production paths, not merely selected file names.

| Epoch | Identity and meaning |
| --- | --- |
| Measured executable, both protocols and pilot | SHA256 `5001bb6b8bdfddc60e74e9fb9bee391e177b1a9992851f47e9ae1a3c8f53c2dc`; ordinary release rustc 1.97.1, aarch64-apple-darwin, LLVM 22.1.6. Each directory retains a frozen copy; the mutable shared target executable is not provenance. |
| Primary driver | Loaded SHA256 `5764f6b110841d181f6ce90f670ecca020441aeae5e9d334451fe08df965ef82`. Launch HEAD `f13aa727a507f543717af3872202324859adcb78` had new driver/tests untracked; copies/hashes and a separately timed postlaunch status observation are retained. Later commit `4b63e9b` preserves measured bytes without relabeling launch as clean. |
| Direct numerical audit | Frozen startup source SHA256 `c1a312eef14a74aa5769b09c374c459d198b88397923ccd0f64db18293f03d78`, epoch direct-oracle-v2. Primary/audit startup Python identity was not captured; a timestamped postlaunch Python 3.14.5 observation is explicitly later evidence. |
| Snapshot-free driver | Clean launch HEAD `0c54db355b830c49e687e97e2e9aa8502bbaa5d1`; loaded SHA256 `fc3baeb1cc309667fadd6b0160bf6de4b682ed3a35f89c5380c40ad1853df58b`; Python 3.14.5 captured at startup. Same measured executable/configuration/steps, separate protocol. |
| Post-timing strict saved-state types | Source SHA256 `d673ba1e01aa7cd0656a70f5e007546180ac630f6ffb834c458f8917438e9b11`, corresponding to `9497353`; metadata-only certification covers all 4,032 primary saved states. It is not the source that ran timing or the earlier physics audit. |
| Review fixes and final tested tooling | Task 2 fixes `e3ef012`/`885a2ac`; current driver SHA256 `1b36b053ae0c2022ef2441765fd929f1ca510fda0928b2738078a46ea6119a03`; test SHA256 `1711f2f71f7199d865e3bb5fc798230599c81614ab51f57cd17b5084e5470269`; CI YAML SHA256 `a65ce79cab6787f13b9686e0d3bdc13abd50b1120db2a4a742271723fa343175`. No performance rows came from this tooling epoch. |
| Fresh contract-only executable | Task 2 CI-consumer fix built SHA256 `822b0005de5566ea855b2e955a9a4a0f8253f9e4a3dbc2c13c8b24ed78e3a98b` in an isolated fresh tracked-source clone with explicitly recorded CI/docs overlay. Task 3 supplies this frozen executable and its actual provenance explicitly to the real Rust contract gate; it was never the measured executable. |

Unchanged Rust SHA256 values: `src/packed_pair.rs`
`a1ea00b373f2612746294d860240dee27cae438db6ef794770ae9850ddeb7702`,
`tests/packed_pair.rs`
`559558f890633c1bee298f64eac14965fe4ee667c2ec426ec5cc082352fb37e4`,
`examples/layout_step.rs`
`4ab17b29958d2c3fc7146ee5fa2089079a2eb5cd99218baa95fcffceb7db27f2`.
Paths are under `crates/modena-bench`. Complete tracked hashes are saved in the
Task 3 source receipt, rather than inferred from branch names.

## Fresh Task 3 checks, 2026-11-06

These checks verify the final benchmark/tooling diff; they are correctness/build
checks, not timings. Production is untouched, so no broad production scientific
campaign, GPU or MPI runtime check was rerun. The scoped Rust bundle naturally
includes its existing full-step/convergence test. No independent numerical gate
was weakened or replaced.

All commands used the B03 worktree above. Agent-shell commands use RTK;
preserved runner scripts and product subprocess argv invoke tools directly.
Fresh logs and command/environment/source receipts live in
`target/b03-evidence/task3-closure-20261106`. The shared target symlink resolves
to `/Users/ricardoleal20/Documents/Code/Misc/modena/target`; compilation instead
used isolated `CARGO_TARGET_DIR=/private/tmp/modena-b03-task3-target-20261106`,
`CARGO_NET_OFFLINE=true`, `RUSTFLAGS=-D warnings`. No shared-target binary or
historical manifest was overwritten.

The Python discovery used these explicit matching prerequisites:

```sh
rtk proxy env \
  B03_LAYOUT_STEP_BINARY="$PWD/target/b03-evidence/task2-fix2-20261106/ci-runtime/b03-layout-step-contract/layout_step" \
  B03_LAYOUT_STEP_PROVENANCE="$PWD/target/b03-evidence/task2-fix2-20261106/ci-runtime/b03-layout-step-contract/provenance.json" \
  B03_TEST_EVIDENCE_DIR="$PWD/target/b03-evidence/task3-closure-20261106/tooling-tests" \
  /Users/ricardoleal20/.pyenv/versions/3.12.8/bin/python -m unittest discover -s bench/tests -v
```

These local frozen paths are historical prerequisites, not files guaranteed in
a fresh clone. For fresh builds, use the complete [build/freeze/provenance recipe](B03_PACKED_LAYOUT.md#required-actual-rust-executable-test)
or the committed CI preparation step. The contract requires both matching
paths, checks real binary/compiler/source provenance, and never silently skips.
`B03_TEST_EVIDENCE_DIR` must be fresh when rerunning artifact-retaining tests.

| Actual command, with isolated environment above | Exact result | Log |
| --- | --- | --- |
| `rtk proxy /Users/ricardoleal20/.pyenv/versions/3.12.8/bin/python -m unittest discover -s bench/tests -v` | Exit 0; 47 tests, 21.367s, 0 failures/errors/skips. Includes all 10 layout tests and genuine interleaved-serial/Cluster4 native example invocations. | `bench-python312.log` |
| `rtk cargo test --locked --offline -p modena-bench --lib --test packed_pair --example layout_step -- --nocapture` | Exit 0; library 0 tests, integration 7 passed, example 4 passed; each reports 0 failed/ignored/measured/filtered. | `bench-rust.log` |
| `rtk cargo fmt -p modena-bench -- --check` | Exit 0; no output. | `bench-fmt.log` |
| `rtk cargo clippy --locked --offline -p modena-bench --lib --test packed_pair --example layout_step -- -D warnings` | Exit 0; 0 compiler/clippy warnings. | `bench-clippy.log` |
| `rtk cargo +1.90.0 check --locked --offline -p modena-bench --lib --test packed_pair --example layout_step` | Exit 0; 0 warnings; rustc 1.90.0, aarch64-apple-darwin, LLVM 20.1.8. | `bench-rust190.log` |
| `rtk proxy /Users/ricardoleal20/.pyenv/versions/3.12.8/bin/python /private/tmp/b03-task3-receipts.py` | Exit 0; all 3 historical manifest identities match; 21 selected members match their original hashes/sizes. All 192 case decisions recomputed directly from retained paired rows, including exact AB/BA order and medians; 0 timing passes. | `historical-receipts-and-cases.json` |
| `rtk git diff --cached --check`, `rtk git diff --check`, scoped source projection and local-link check | Exit 0; only the four owned documentation files differ; all 37 local link targets exist; reviewed Rust projection remains unchanged. | `document-checks.json` |

`commands.json` records exact argv, resolved environment, UTC start times, elapsed
check durations and exits. `source-before.json`/`source-after.json` establish
tested HEAD, core projection and unchanged tracked bytes through the checks.
The runner sources are copied beside those receipts. Python is 3.12.8;
Rust tests/clippy use 1.97.1. No warnings appeared in these check logs.

The 10 layout tests exercise real subprocess failures and cover malformed
schema/thread types, positive elapsed time, zero-time summary handling, wrong
metadata/states, nonfinite values, wrong saved-force/FNV data, direct-oracle
rejection of shared bias, secondary binary mismatch and output reuse. The new
strict-metadata, zero-elapsed and actual-Rust-contract test artifacts are retained
under `tooling-tests`; other regression fixtures use bounded temporary
directories. These fixtures are distinct from the retained campaign records. The
reviewed snapshot-row finding is closed by the documented exact subset contract
and real Rust executable test: snapshot.row equals stdout minus only
`pool_setup_seconds`, with exact shared JSON types/values. The unbuffered
snapshot writer observation remains for whole-branch review triage; no writer
change or new measured executable is included here.

Task 2's CI-consumer fix separately validated the exact YAML build/freeze block
and 10 covering tests in a fresh tracked-source clone with a recorded overlay.
That local check is not the whole CI job. Remote GitHub Actions is unverified
under the known account/billing lock; no remote green claim is made.

## Artifact receipts and limits

The historical controller independently read all 6,874 manifest members,
2,871,575,609 bytes, verifying sizes and SHA256. Its receipt is
`target/b03-evidence/controller-task2-integrity-9497353.json`.
Task 3 freshly checked the following manifest-file identities and 21 selected
members (including frozen binaries, row indexes, summaries and audit reports),
without rereading every saved state or rewriting any historical file:

| Directory under `target/b03-evidence` | Historical verified members/bytes | Unchanged manifest SHA256 |
| --- | --- | --- |
| `task2-campaign-20261106` | 4,074 / 2,791,589,229 | `19a7e7f57e42619fa43233b8c8d71a34931b2f8ac6c35b48a7d6748f0cdd56b8` |
| `task2-no-snapshot-20261106` | 2,700 / 12,372,708 | `294a92657dd0f5a81270105557f566f482a6b29635d93be4c80faa4fcaf9bc61` |
| `task2-pilot-20261106` | 100 / 67,613,672 | `6ad4740bc2f7974570f2572c94fc83ab3562cb2057c17912e64ff3d0271395f8` |

The new Task 3 directory has its own SHA256 manifest, commands/logs, document
copies, all 192 explicit case outcomes and the receipt sources. Generated
artifacts remain ignored. Original stdout/stderr bytes and losslessly compressed
snapshot JSON remain retained; both raw/compressed snapshot hashes are recorded.

Interpretation remains bounded:

- Primary RSS includes retained initial snapshot Values and instrumentation
  cache effects. Snapshot writing occurs outside the full-step timer; direct
  unbuffered File output is a source/telemetry inference about collection cost,
  not a measured syscall count. Snapshot-free rows provide a separate protocol,
  with no fresh full vectors or full-state SHA. Independently recomputed FNV is
  noncryptographic and does not substitute for the primary state audit.
- At N 10000, SoA/AoSoA4 retain 468.75 KiB packed capacity; Cluster4 retains about
  21,537–21,545 KiB. Baseline packed bytes 0 does not mean zero total memory.
  Across pools/policies primary baseline/Cluster4 RSS ranges 172.98–179.20 /
  215.64–225.72 MiB; snapshot-free ranges 17.16–31.47 / 53.64–86.52 MiB. These are
  all-row ranges; matched medians remain per case. Cross-session RSS/wall
  differences are not causal kernel speedups.
- Orthorhombic periodic, single-species reduced shifted-force LJ,
  sigma=epsilon=mass=1, cutoff 2.5, skin 0.3 is the declared domain. Truncated FCC
  N 5000/N 10000 fixtures contain vacancies and are not equilibrated ensembles.
  Unsupported inputs reject. No chemistry/engine parity, general layout benefit,
  SIMD/parallel packed performance or full GPU/MPI runtime claim follows.

Production storage, API, pair dispatch, deterministic path and defaults remain
unchanged. The delivered research result is reproducible rejection of these
scalar prototypes in the declared cases; production promotion is not supported.

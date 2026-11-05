---
status: passed
---

# B06 implementation verification

Local implementation on 2026-11-05, verified in the working tree based on `f166399070b2d6ebfd91f51e5934c132f04d46c0`. The user explicitly authorized “Aplica el B06” and subsequently requested its PR and integration. This record covers the [B06 diagnostic contract](B06_ENERGY_DIAGNOSTICS.md), separately from the earlier research report and historical scientific results. The results below were captured before creating the B06 commit and PR.

Environment: macOS 26.6.2, arm64; Rust 1.97.1; Python 3.14.5. The compiler-floor checks additionally use installed Rust 1.90.0. Cargo commands use `--locked --offline`; no optional hardware runtime validation is inferred from compilation.

## Automated checks

Commands below were executed locally through the required `rtk` proxy.

| Check | Result |
| --- | --- |
| `cargo test --locked --offline --workspace` | 296 passed, 1 ignored; 36 suites, 170.83 s |
| `cargo test --locked --offline -p modena-mcp` after adding JSON policy/real-unit assertions | 42 passed; 7 suites |
| Python client tests against the release CLI | 21 passed |
| `python3 -m unittest discover -s bench/tests -v` | 22 passed |
| `cargo clippy --locked --offline --workspace --all-targets -- -D warnings` | Passed |
| Same Clippy command with `--features modena-core/gpu` | Passed; compilation/lint only |
| `cargo +1.90.0 check --locked --offline --workspace --all-targets` | Passed |
| Same compiler-floor check with `--features modena-core/gpu` | Passed; compilation only |
| `cargo fmt --all -- --check` | Passed |
| `git diff --check` | Passed |

The MCP package rerun includes one new request-path test added after the workspace run, plus strengthened assertions for real-unit display versus reduced diagnostic energies. Counts are separate command results, not an aggregate. The ignored workspace item is the pre-existing illustrative Rust builder doc block marked `ignore`.

The benchmark-tool suite uses `/usr/bin/time -l`. Its kernel clock query is denied inside the sandbox; the suite passed when the automatic approval reviewer allowed the local test command outside that restriction. Timing scripts and failure assertions were retained. This test-suite run is not a performance comparison.

Focused regressions cover:

- Energy-offset invariance of the normalized metrics and selected verdict with a nonzero endpoint change; the legacy relative number intentionally differs.
- Analytic linear energy/slope/RMS and a bounded excursion returning to its starting energy.
- Every-step statistics with sparse logs, NVT/MC-NPT applicability, zero kinetic reference, invalid observations and selected/invalid threshold policies.
- CLI JSON output; Python transport, NVT and explicit policies; exact older-CLI flag rejection without rerunning a failed simulation.
- MCP full-rate status, finite last logged samples followed by an unlogged nonfinite failure, failed-step telemetry and anomaly evidence.
- Sampling means/variance against independent expectations, correlated versus effective samples, explicit and heuristic discard, energy offsets, constant/nonfinite series, capped autocorrelation tails, replica separation and actual CSV transport.
- Validation failure artifacts for malformed/scalar JSON and invalid sampling, retaining raw output.

An independent read-only code-review agent found three issues during implementation: older CLI compatibility, lost failure telemetry at sparse cadence and scalar JSON roots bypassing a validation artifact. All three were reproduced by regressions, corrected and retested. Its second review reported no Critical/Important findings and one public-comment omission; that comment now includes `run_failure`. Automated review does not substitute for scientific validation.

## Executed scientific smoke

```bash
cargo build --locked --offline --release -p modena-cli -p modena-core \
  --example scientific_validation --bin modena
python3 bench/validate.py --output target/b06-scientific-smoke
```

Result: **67/67 existing scalar observations passed**, with schema version 2, six analyzed series and two groups of three replicas. Scalar expectations and tolerances were not relaxed. Raw output is retained locally in [validation.stdout](../target/b06-scientific-smoke/validation.stdout); the structured result and provenance are in [validation.json](../target/b06-scientific-smoke/validation.json). These `target/` artifacts are build outputs, not tracked source files.

Provenance captured by the run at `2026-11-05T06:53:26.021088+00:00`, one Rayon thread:

| Artifact | SHA-256 |
| --- | --- |
| `scientific_validation` release executable | `d6ef67f94c0c166251ed7d3bfb8eb286af3680041ed5477cd8398dbdfbaa4373` |
| `bench/analyze.py` | `b09b9022527258bd797bef86c9d8b3b9e8b156042424130274b2490decccd6d4` |
| Raw `validation.stdout` (1,534,589 bytes) | `394797f8789becbbcad3fa556f29154851875bf8be25951e4d779dd87ad8dac8` |

The payload also records the base commit, dirty working-tree status, command, environment and raw stderr. The command was executed on the modified working tree; the base commit alone does not identify the new implementation.

### Two-body NVE fixture

All runs use the same smooth LJ Hamiltonian, initial state and reduced physical duration of 2. Initial kinetic energy is 0.01. They deliberately remain `uncalibrated`, without a universal advisory policy.

| Reduced timestep | Samples including initial state | Maximum absolute energy excursion | Population RMS | Energy/time slope |
| --- | --- | --- | --- | --- |
| 0.004 | 501 | 4.53571794e-5 | 1.83306776e-5 | 1.71949258e-6 |
| 0.002 | 1,001 | 1.13386128e-5 | 4.58359559e-6 | 4.33042105e-7 |
| 0.001 | 2,001 | 2.83474203e-6 | 1.14604815e-6 | 1.08668280e-7 |

Excursion ratios on timestep halving are **4.0002406050** and **3.9998746712**, within the existing `4 ± 0.3` criterion. This fixture demonstrates the expected second-order reduction of a bounded energy error; its nonzero fitted slopes do not establish a secular trend. It does not validate every potential, timestep or material.

### Observable sampling fixtures

Each stored series has 13,667 rows with stride 3; the explicit burn-in discards 334 rows and retains 13,333. The CSVR stored interval is 3 reduced time units, and rigid-rotation interval is 0.03. Scalar moments still use all 40,000 production steps per seed. Values below summarize the stored-series estimates, which differ slightly from the full-rate scalar moments.

| Observable | Seed | Mean | Estimated g | Estimated effective samples | Conditional mean standard error |
| --- | --- | --- | --- | --- | --- |
| CSVR kinetic energy | 7 | 5.99732243 | 1.033666 | 12,898.75 | 0.02150263 |
| CSVR kinetic energy | 42 | 6.00812808 | 1.156286 | 11,530.88 | 0.02311348 |
| CSVR kinetic energy | 91 | 5.98959212 | 1.073862 | 12,415.94 | 0.02212832 |
| Rigid-rotation kinetic energy | 7 | 1.53190565 | 19.572177 | 681.22 | 0.04743185 |
| Rigid-rotation kinetic energy | 42 | 1.50588911 | 19.523974 | 682.90 | 0.04671778 |
| Rigid-rotation kinetic energy | 91 | 1.51192247 | 15.767470 | 845.60 | 0.04209461 |

The observed CSVR replica-mean range is `[5.98959212, 6.00812808]`, with between-replica standard error `0.00537536`. The rigid-rotation range is `[1.50588911, 1.53190565]`, with between-replica standard error `0.00786197`. These min/max ranges are not confidence intervals. All six retained series had empty analyzer warning lists; this is not proof that slow modes are absent or that the ensemble is correct.

No new external LAMMPS, DFT, ML, GPU-runtime or MPI-scaling campaign was performed for B06. The measured outcomes establish bounded diagnostic behavior and regression checks; they do not establish general material accuracy, engine replacement or speedups.

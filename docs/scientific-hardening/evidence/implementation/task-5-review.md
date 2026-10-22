### Spec Compliance

- ❌ Issues found: the delivered GPU benchmark wrapper does not retain the required failure metadata when its executable is absent or an OS launch error occurs (`bench/gpu_ladder.py:32`, `bench/gpu_ladder.py:44`). The CPU benchmark, scientific exporter, bounded CI, documentation corrections, archive, and measured CPU campaign satisfy the reviewed Task5 requirements.
- ⚠️ Cannot verify from this task diff: all Task1–4 production invariants and fresh whole-workspace/Python/release/GPU-adapter/MPI2 validation. The coordinator owns those checks; their pending status is explicitly disclosed in `docs/scientific-hardening/RESULTS.md:93`, and is not a missing Task5 deliverable. The absence of concurrent work during historical timings is supported by the handoff report but cannot be independently reconstructed from the saved timing logs alone.

### Strengths

- `bench/benchmark.py:104` retains actual input decks/hashes, stdout/stderr, return codes, CPU/RSS/null values, thread count, precision, DOF and drift; `bench/benchmark.py:169` accurately labels whole-process timing and separate-process warmup. Executable fixture coverage and the real CLI smoke are concrete (`bench/tests/test_benchmark.py:20`, `bench/tests/test_benchmark.py:48`); the saved final run records five passing tests with no skips (`docs/scientific-hardening/evidence/raw/tooling-release.log:1`).
- Numerical gates are real observations with analytical expectations, multiple seeds, and stated tolerances (`crates/modena-core/examples/scientific_validation.rs:54`, `:130`, `:189`, `:213`, `:257`, `:304`). The exporter includes constrained rigid rotational sampling, dt convergence, direct-Ewald mesh references, full-cell strain and a stochastic JSON seam. All 67 saved rows independently satisfy their numerical gates (`docs/scientific-hardening/evidence/scientific-release/validation.json:572`); the initially failed coarse mesh remains preserved and acknowledged (`docs/scientific-hardening/RESULTS.md:62`).
- Saved benchmark evidence has 24 zero-returncode processes, six warmups and 18 measured runs. Recomputed all six median timing/RSS groups and all 24 deck hashes: they agree with the publication (`docs/scientific-hardening/RESULTS.md:112`, `docs/scientific-hardening/evidence/matched-lj/summary.json:1`). Clean source commit and binary hashes are recorded (`docs/scientific-hardening/evidence/matched-lj/benchmark.json:3`, `docs/scientific-hardening/evidence/scientific-release/validation.json:573`). Retained current CLI, LAMMPS, scientific-example and probe-library hashes match their respective records.
- The archive manifest matches all 46 archived files and their original audit sources byte-for-byte (`docs/scientific-hardening/evidence/original-audit-manifest.json:1`). Historical labels/whitespace are evidence, not current code defects. The current probe results are tied to the original probe source hash and release library hash (`docs/scientific-hardening/evidence/raw/original-probes-current.json:67`).
- Domain and compatibility claims are appropriately narrowed, including orthorhombic SPME, diagnostic constrained pressure, experimental acceleration, corrected real-charge normalization and lack of strict historical speedup (`README.md:8`, `docs/scientific-hardening/UNITS.md:3`, `docs/scientific-hardening/RESULTS.md:118`). The deferred tabulated header now correctly describes joint Hermite energy/force interpolation (`crates/modena-core/src/potential/tabulated.rs:5`). CI explicitly builds and runs bounded numerical validation and uploads its evidence (`.github/workflows/ci.yml:103`).

### Issues

#### Critical (Must Fix)

- None found.

#### Important (Should Fix)

- **I1 — GPU launch failures lose the required evidence.** `bench/gpu_ladder.py:32` calls `sha(binary)` unconditionally before an artifact exists. Running the documented optional tool without its prebuilt executable raises `FileNotFoundError` and leaves no `gpu.json`, return code, or retained stderr. The unprotected `subprocess.run` at `bench/gpu_ladder.py:44` similarly loses evidence on an OS launch error (including a missing `/usr/bin/time`). This fails the task requirement that failed benchmark runs retain machine-readable failure provenance; it does not invalidate the recorded CPU timings or kernels. Use the main benchmark/validator failure pattern: nullable unavailable hash, protected process launch, retained error/logs and nonzero return code, and a timer fallback when unavailable. Add a small executable GPU-wrapper missing-binary/failing-process fixture so this path cannot regress.

#### Minor (Nice to Have)

- **M1 — Criterion integration timing includes teardown that its scope omits.** `crates/modena-bench/benches/integrator_step.rs:25` takes ownership of `System` and `NeighborList`, then returns `()` after `black_box(e)`. Both allocations are destroyed inside the measured closure, so the sample includes buffer deallocation as well as the promised single steady step. `bench/BENCHMARKS.md:28` lists setup outside timing and step operations inside, but omits this teardown. Return the owned inputs/result for Criterion to dispose outside the timed routine, or explicitly document teardown as part of the measurement. No published numerical Criterion performance claim currently depends on this.
- **M2 — Nonfatal benchmark output noise is retained.** `docs/scientific-hardening/evidence/raw/criterion-force-smoke.log:3` and `docs/scientific-hardening/evidence/raw/criterion-step-smoke.log:3` report “Gnuplot not found, using plotters backend.” All 18 force and nine integration cases still report success; this is fallback noise rather than a failed physics/build check. State the fallback in the check report or explicitly configure the available plotting backend for future runs; preserve the original logs.

### Assessment

**Task quality:** Needs fixes.

**Reasoning:** The scientific checks and evidence provenance are carefully bounded and substantiate the published results. Fix the GPU wrapper failure-provenance path before treating all delivered benchmark tools as satisfying the failure-sensitive contract; the remaining findings are scope/reporting polish.

### Review checks and boundaries

- Reviewed the supplied three-commit diff package in targeted passes over tooling/CI, scientific example/Criterion, documentation, and evidence metadata. Historical raw logs were inspected as provenance, not regraded as current source. No changed source was separately reread except the GPU producer CSV output segment: its `main` hunk ends before emission, so a focused check of `crates/modena-core/examples/gpu_ladder.rs:165` through `:248` was necessary to verify the new wrapper’s row-selection contract. The header and numeric rows match; the summary row is correctly excluded.
- Focused outside-diff risk checks: archive preservation checked the 46 original paths named by the manifest; campaign/binary identity checked only the binaries and release library named by saved provenance. No repository crawl, Git mutation, suite rerun, build or benchmark campaign was performed.
- Static inspection establishes I1 and M1 directly; existing saved runs do not cover missing GPU executables. Only this review report was written. Coordinator should retain it in the durable evidence archive after fixes/review closure.

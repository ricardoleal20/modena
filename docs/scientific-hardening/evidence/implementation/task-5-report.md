# Task5 report: reproducible scientific evidence

## Immutable handoff

BASE3bf7f54c683115aabc366d59f3d7d33b36669297; branch
ricardoleal20/modena-scientific-hardening. No publication/merge. Commits:

05d7d8f6a0b811c09f701a42071862a49bc0cf96 :memo: Record: matched resource timings and scientific evidence
11d7683841786699e77d02d7c55c2d8af9cdac40 :white_check_mark: Add: reproducible scientific validation and benchmark tools

An additional report-only archival commit follows; it changes no measured code.
Root owns fresh full-workspace/Python/GPU-adapter/MPI2 checks and independent
whole-branch review after this handoff. This report does not claim those future
checks passed. Source11d7683 was clean for final timing and numeric campaigns;
subsequent evidence/CI-release-binary environment changes are not measured kernels.

## Scope and deliverables

Owned benchmark wrappers/runner/decks, stdlib tooling tests, Criterion workloads,
scientific exporter/validator, CI bounded smoke, README/ARCHITECTURE/results/unit
and water/legacy validation docs, durable evidence. Integration changes: GPU ladder
sets conserved COM metadata3 before init_velocities and declares resources/scope;
only deferred tabulated header corrected to describe joint Hermite interpolation.
No speculative production engine rewrites; no subagents.
Preserved other agents' changes.

`bench/benchmark.py` uses positive actual requested steps/cadences, generated matched
LJ shifted-force decks, explicit precision/threads/DOF/seed/parameters, separate
fresh-process warmup and measured repetitions. Retains stdout/stderr/inputs/hash,
returncode, wall/user/system CPU/RSS with unavailable null, drift and binary/compiler/
hardware/Python/commit/dirty/load provenance. Failure exits1 after writing artifact;
missing binary becomes127. Renderer converts numeric sizes/actual steps rather than
string integer formatting/hardcoded2000. CLI whole-process startup is included;
no trajectory/thermo file writers enabled. MODENA/LAMMPS RNG algorithms differ.
Only MODENA1t vs serial LAMMPS is marked resource-matched. GPU tool retains outer
resources and CSV integration rows with distinct in-process timing scope, f32 kernel/
f64 host, COM3 and different-Hamiltonian gather probe caveat.

Criterion18 force cases (9 workloads ×cached/rebuild) +9 integration cases:
LJ N256/864/2048, ionic SPME N64/216/512, molecular SPME N81/192/375.
Cached-force scope includes force clearing/model compute/neighbor displacement
checks; forced-rebuild scope adds rebuild. Neighbor-build bench isolates list
creation/build. Integration setup builds cached neighbors outside measured one-step
window; timed step includes kicks/force clearing/compute/checks. Molecular fixtures
benchmark electrostatics of SPC/E-like geometry/exclusions, not complete water or
rigid integration. RAYON_NUM_THREADS controls Criterion resources.

Scientific example exports actual numerical observed/expected/tolerance/basis:
CSVR stationary Gamma mean/variance for3seeds; constrained/unconstrained BAOAB
free-spectator Maxwell moments for3axes/3seeds; rigid triangle rotational Gamma
moments/residuals; LJ finite differences; NVE same-duration dt-halving bands;
SPME actual mesh15³/20³/direct-Ewald comparison and isotropic full-cell strain
virials; exact stochastic JSON13+17 continuation. 40000samples/seed plus1000burn;
optional --extended400000samples. Correlated moments are explicitly not formal
uncertainties; rigid bands use conservative OU inflation. CI runs smoke+tool+Python
regressions with explicit release binary, uploads numerical evidence on failure.

## RED/GREEN and observed checks

All shell commands use RTK; raw output captured through `rtk proxy`.
Supported Python3.12.14:
/Users/ricardoleal20/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/bin/python3
No installations/pytest dependency. Five executable unittest cases pass.

- Initial failure/success metadata fixture RED: new benchmark script absent and
  artifacts not produced; tooling-red.log. GREEN runs actual runner/subprocess
  fixtures exit17/0, verifies code/logs/step overrides/precision/DOF/Python and RSS
  never fabricated0. Both metadata tests execute renderer through success fixture.
- Validator nonfinite/out-of-band RED absent artifact; validation-tool-red.log.
  GREEN synthetic actual executable outputs measured3 vs expected1±.1 and NaN;
  runner fails and retains numerical/raw provenance. Missing binary additional
  meaningful RED traceback/no metadata in missing-validator-red.log; GREEN produces
  exit127, null binary hash and artifact with overall failure.
- Real release CLI RED: generated log/trajectory cadence0 rejected by strict
  validator (real-cli-tool-red.log and cli-smoke raw JSON/logs). Fixed to steps+1
  with file paths absent. GREEN actual release3steps/1warmup has drift and no CSV/
  XYZ writes. This is covered in5tests, no skip.
- Example absence recorded scientific-example-red.log; initial numerical58-row
  campaign56/58 passed, coarse kcut3.2 mesh failed declared.005 E/force limits
  (.00856/.01242 errors). Preserved scientific-smoke/validation.json unchanged.
  Acceptance uses predeclared.005 gates on kcut4.8/6.4 meshes15³/20³ and adds rigid
  mode rows; final67/67 pass with no tolerance weakening to hide coarse failures.
  A resolution-accuracy boundary is reported, not called fixed production code.
- Criterion compilation initially exposed NeighborList not Clone; bench-compile.log.
  Fixed benchmark untimed setup to build list, no engine clone addition.
  cargo bench -p modena-bench --no-run passes (bench-compile-final.log).
  cargo bench -p modena-bench --bench force_kernel -- --test passes18cases.
  cargo bench -p modena-bench --bench integrator_step -- --test passes9cases.
- cargo build --release --bin modena passed; CLI artifact used by tooling/benchmark.
  cargo build --release -p modena-core --example scientific_validation passed.
  cargo clippy -p modena-core --example scientific_validation -p modena-bench
  --benches -- -D warnings passed. cargo fmt --all -- --check passed.
  cargo check -p modena-core --features gpu --example gpu_ladder passed.
  No new actual-GPU/MPI suite duplicated; final root will exercise them.

Exact tool interfaces:

```sh
rtk proxy <Python3.12path> -m unittest discover -s bench/tests -v
rtk proxy <Python3.12path> bench/validate.py --output target/scientific-hardening-clean-numeric
rtk proxy <Python3.12path> bench/validate.py --extended --output target/scientific-extended
rtk proxy <Python3.12path> bench/benchmark.py --lammps /opt/homebrew/bin/lmp_serial --ladder 864 2048 20000 --threads 1 --steps 1000 --warmup-steps 200 --runs 3 --output target/scientific-hardening-matched
rtk proxy <Python3.12path> bench/gpu_ladder.py --ladder 864 --steps 100 --runs 3 --threads 1 --output target/gpu-new
```

No build/test concurrency during matched timing: sent root readiness notice and
pgrep confirmed no cargo/rustc/validation/probe processes before launch. Root
reported no other active heavy implementation/review. Kept all heavy work stopped
until24/24 subprocesses finished, then informed root. Timing command returned0.

## Actual numerical and timing results

Durable docs/scientific-hardening/RESULTS.md (Spanish) cites machine JSON/raw logs.
Final numerical provenance/source11d7683 clean,67observations passed. CSVR K
means6.000882/5.984838/5.987383 vs6±.085; variance5.990228/6.053033/5.984279 vs6±.25.
Rigid Krot means1.533384/1.503093/1.515532, variances1.533444/1.485694/1.506766 vs1.5;
constraint residual<1e-10. NVE energy-band ratios4.0002406/3.9998747 vs4±.3.
SPME maxforceerror.00233876/.000594538, full strain FD differences<1e-7.
Restart JSON seam maxerror0. All BAOAB Maxwell axis moments and LJ FD rows are
present individually in evidence/scientific-release/validation.json with tolerances.

Re-ran exact original audit probes source against current default release core.
Direct rustc linking first failed Apple LLVMbitcode compatibility; retained failure.
-C lto=thin -C codegen-units=1 succeeded. Run RAYON_NUM_THREADS1 gives CSVR Tmean
.99933785,Tvariance.16712408 vs1/1/6; freeatom var1.00073504 constrained/unconstrained;
barostatPE/cache-force errors0 at printed precision. Machine row tolerances and
library/source hashes in raw/original-probes-current.json; raw original bug
prediction labels remain unchanged and are not treated as current expectations.

Matched whole-process median ms/step MODENA1t/LMPserial:
N864 .601879/.309412 (ratio1.945),N2048 1.738660/.736351 (2.361),
N20000 17.326401/6.525943 (2.655). 1000steps,3reps; separate200step warmup;
24processes exit0,6warmup+18measured. AppleM3Pro11logical; Rust1.97.1;
LAMMPS22Jul2025Update6; matched shifted-force rc2.5/skin.3/a1.6796/dt.005/T1/3N-3DOF.
Per-run raw wall/user/systemCPU/RSS/drift/logs retained under evidence/matched-lj.
Important concern: N20000 MODENA RSSmedian427212800B vs LMP22003712B. Reported
as measured gap, no unmeasured attribution/no engine optimization in Task5.
Relative drifts valid for both engines, not superior-conservation claims.
Historical baseline preserved; different DOF/microstates/shared desktop preclude
strict before/after speedup. Multicore default comparisons are not fair1t speed wins.

## Honest limitations and archive

Real-charge ruling copied into plan scope/RESULTS/UNITS/WATER_MODEL: correct
sqrt(C/(sigma*epsilon)), old sigma² physical energy error, reduced/sigma1 unchanged,
old real electrostatic checkpoints at sigma!=1 normal identity mismatch rejects.
No historical migration claim, no corrected-model comparison to old real timing.
Unsupported quantum/reactive/manybody/fullHPC, full anisotropic/constrainedNPT and
unvalidated constraint pressure explicit. SPME orthorhombic; directEwald/pairs
triclinic, domain policies/box moves restricted, GPU f32 hostoffloadexperimental.
Legacy LAMMPS inconsistent-energy/exactconservation/one-drawprojection/relaxation-
explains-NVE claims withdrawn. Short regression sampling is not converged liquid
observables or formal uncertainty. JSON seam observes component stochastic state;
full physical-identity/file/extended-state coverage lives in prior task suites and
will be freshly verified by root, not mislabelled as this67-row exporter.

Byte-identical original audit sources/baseline/probes/logs copied with manifest;
Task1–4 implementation reports/reviews/ruling and relevant raw logs durable under
implementation; originals of modified BENCHMARKS/VALIDATION/WATER_MODEL archived.
No binaries/large trajectories tracked; retained text/JSON ~1MB with commit
rationale for independent review after ignored workspace cleanup. This report is
also copied durable; final review reports should be archived by root when available.


## Task5 fix round 1 (from immutable de8529a)

Read task-5-review.md and addressed all actionable findings in assigned scope.
I1: GPU wrapper catches unavailable/unreadable binary hashing (nullable hash and
hash-error text), creates gpu.json before launch, records actual launch_command
and selected resource_timer, catches OS launch errors as return code127 with
retained stdout/stderr/artifact, and returns failure1. Existing process exit codes
are preserved. New optional --timer PATH defaults /usr/bin/time; absent timer
falls back to direct launch with CPU/RSS null. No engine/bench kernel rewrite.

Meaningful executable RED/GREEN: added five tests in bench/tests/test_gpu_ladder.py.
Initial RED failed5/5, including the actual missing-binary SHA traceback. After
introducing only the timer selection interface, the second RED still failed5/5
with actual FileNotFoundError/PermissionError OS launches and absent artifacts.
Both raw RED logs retained. GREEN5/5,0skips tests actual missing executable,
exit17fixture with stdout/stderr, nonexecutable engine, nonexecutable timer, and
successful CSV-producing fixture with absent timer. All assert artifact/logs,
return codes, nullable missing hash or available hash, steps/thread/DOF and null
unavailable CPU/RSS. Tests launch real small executable fixtures, not GPU/CPU
scientific or performance campaigns. Final command:

rtk proxy /Users/ricardoleal20/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/bin/python3 -m unittest discover -s bench/tests -p test_gpu_ladder.py -v

5tests passed with no skips; exact outputs under
raw/task-5-fix-round-1-gpu-green-final.log. py_compile wrapper/test passed.
An explicit missing-GPU invocation returned1 and retained gpu.json with hashnull,
process127 and logs under evidence/task-5-fix-round-1-missing-gpu; outer command
stdout/stderr also retained. This artifact correctly labels dirty source de8529a
plus the uncommitted fix rather than claiming a clean future commit. Invocation:

rtk proxy <Python3.12path> bench/gpu_ladder.py --binary target/absent-gpu-fixture --engines cpu --ladder 864 --steps 7 --runs 1 --output docs/scientific-hardening/evidence/task-5-fix-round-1-missing-gpu

M1: explicitly documented that Criterion integration closure destroys owned cloned
System/NeighborList buffers before returning; measured scope includes teardown/
deallocation in addition to step. This is documentation only; no kernel changed,
no Criterion performance values published, no retiming required.
M2: check report/Spanish RESULTS now explicitly state nonfatal
“Gnuplot not found, using plotters backend” fallback in Criterion smoke output.
18force/9integration successes remain the prior observed result; raw bytes preserved.

Original measurements, engine binary hashes,67-row numerical campaign and24-run
matched timings remain untouched. Verified all24 recorded input hashes unchanged.
Scoped whitespace check excluding byte-preserved evidence passed. Task5 review
archived durable; updated report and all new RED/GREEN/raw failure evidence
archived. No builds/new CPU/GPU/MPI campaigns or subagents in this fix round.
Root full release checks remain pending. No new concern requiring production work.

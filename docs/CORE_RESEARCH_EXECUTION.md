# Remaining core research execution

User authorization 2026-11-06: continue and finish all remaining research cases,
with an active goal. Earlier B01/B06 merged; B09 PR55 and B02 PR56 open.
No further merge authorization inferred. Branches always start `ricardo/`.

Completion requires each backlog acceptance criterion supported by executed
evidence, reviewed changes and an assigned/labeled PR. Hardware/reference
requirements remain explicit; compile-only or synthetic fixtures do not close
scientific/runtime acceptance. Rejected optimization variants are legitimate
experimental outcomes, with measurements and reasons retained.

| Case | State | Dependencies / acceptance focus |
| --- | --- | --- |
| B03 | Delivered in open [PR #57](https://github.com/ricardoleal20/modena/pull/57), assigned/labeled, task and whole-branch reviews approved; all 96 cases reject in each protocol | [Verification and source epochs](B03_VERIFICATION.md), [all-case timing/resource evidence](B03_PACKED_LAYOUT.md). Production unchanged; no SIMD/parallel packed rejection. |
| B04 | Independently active in separate modena-b04 worktree; G1, actual initialized MACE G2 and certified stencil G3 independently approved; private base synchronized with B03; G4 pending | Periodic directed graph images, independent enumeration; keep classical half list distinct. Preserve the active plan and source isolation. |
| B05 | Fixed protocol and pure fixtures/metrics committed on ricardo/b05-spme-calibration; Task1 review approved, no force/cost acceptance | SPME convergence/cost against converged Ewald; declared calibration domain. |
| B07 | pending | Tensor strain derivatives, constraint impulses, periodic and timestep oracles. |
| B08 | pending | Molecular MC volume measure, periodic molecules, full rejection restoration and analytic volume distribution. |
| B10 | pending | Converged DFT provenance and domain-specific ML calibration; actual material/model/reference required for scientific closure. |
| B11 | pending | Resident GPU step and actual device transfer traces, f32 force/convergence and full-step timings. |
| B12 | pending | Persistent local MPI state, typed IDs, migration/halos; actual 2/4-rank trajectories and scaling. |

Source of truth: docs/research/2026-11-05-core-improvements/phase3_analysis/ACCEPTANCE-BACKLOG.md.
Execute per-case design, plan, implementation/review, validation and publication;
do not repeat closed B01/B02/B06/B09 work or overwrite prior evidence. Address
available independent work before asking for an external resource. No spending,
weight downloads or destructive actions inferred from the request.

## Selected domain and resource inventory (2026-11-06)

The user selected **FCC aluminium** for B10 before labels/model evaluation.
The known ASE venv has ASE/NumPy but no Torch/FAIR-Chem; pw.x is absent from
PATH. Local OpenMPI5.0.10 is available, and system metadata identifies Apple
M3 Pro with Metal. These inventory observations do not prove a current GPU
run, MPI trajectory, DFT convergence or ML inference. Fortran16.2 and make
are present; a Quantum ESPRESSO installation still needs preparation.

B04 has advanced from design preparation to independent work in its separate
worktree; its acceptance remains pending. B05 design preparation is available,
with implementation and acceptance pending. B03's packed prototypes use scalar
execution with an equivalent serial interleaved control. Production-auto remains
a separate control; results cannot reject unimplemented parallel/SIMD designs.

## B03 experiment delivery (2026-11-06)

Both frozen protocols retain 1,344 invocations each, plus 32 separate pilot
invocations. All 96 comparisons per protocol fail the fixed >=5% median full-step
gain and >=4/5 favorable-pair criterion; zero independent reproduction runs
trigger. Masks/identity, independent force/energy/virial and full-step/timestep
gates pass in the reviewed evidence. Large Cluster4 retained memory and timing
regressions reinforce rejection. The primary protocol retains snapshot Values
during steps; the separate snapshot-free protocol emits no vectors/full-state
SHA. No results are pooled across protocols and no candidate is promoted.

Final tested implementation epoch is `885a2ac26a3ac8a5f40f51a11dcd5095529ce5da`.
Task 3 freshly passes Python 3.12.8 benchmark discovery 47 tests (including 10 layout
tests and mandatory real Rust executable contracts), scoped benchmark Rust
7 integration + 4 example tests, fmt/clippy and Rust 1.90. Production core remains
tree `39b46c7496c6820219056fe38a821378a32633e8`, and Rust sources are unchanged
since reviewed `66af0e4`. Historical measured executable `5001bb6` and driver/audit
epochs stay intact; the fresh CI-contract executable `822b000` is a distinct epoch.
See [B03 verification](B03_VERIFICATION.md) for exact commands, hashes, counts,
resources and preserved artifact receipts. These checks do not establish a
remote CI result; the known Actions account/billing lock remains.

All three tasks and the final whole-branch review approved B03, delivered in
assigned/labeled open PR #57 against B02. [Review closure](B03_REVIEW.md) records
all controller decisions, the nonblocking writer Minor, workflow archive and
actual Community billing failure before any job step. Main CI's base-main
filter excludes the stacked PR. No additional merge was performed. B04 is
already active; do not restart it or treat B03 as completing the eight-case goal. B05/B07/B08/B10/B11/B12 still require their actual source,
implementation, validation, independent review and PR evidence. B10 remains
FCC aluminium; build-only GPU/MPI checks and unvalidated DFT/ML fixtures cannot
close their scientific/runtime requirements.

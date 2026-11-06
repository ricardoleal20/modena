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
| B03 | Task1 reviewed at66af0e4; campaign pending, ricardo/b03-packed-layout based604039e | Internal packed layouts and cluster masks; include conversion/rebuild, exact logical identity and scalar oracle. |
| B04 | pending | Periodic directed graph images, independent enumeration; keep classical half list distinct. |
| B05 | pending | SPME convergence/cost against converged Ewald; declared calibration domain. |
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

B04/B05 read-only design preparations are ready; their implementations and
acceptance remain pending. B03's packed prototypes use scalar execution and
require an equivalent serial interleaved control. Production-auto remains a
separate control; results cannot reject unimplemented parallel/SIMD designs.

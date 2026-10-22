# Task 2 independent spec and quality review

Review target: `fdbe886886e0eb143f99276c46debcbea836e6d0..9fa3b9e6fefa326a16a473ccade05709d7e76585`.

## Spec compliance: compliant

The implementation covers the concrete Task 2 requirements:

- Joint cubic Hermite interpolation derives radial force from the energy polynomial, uses force-derived endpoint slopes, and extends the lower tangent consistently. Shifted-force cutoff handling is retained. The independent cubic fixture and existing derivative properties cover the new contract.
- Direct reciprocal virial uses the isotropic fixed-kappa strain derivative. Direct and SPME background terms add `3 U_bg`. Tests differentiate energy while scaling both coordinates and cell, and isolate the background with an independently calculated expectation.
- Normalized indexed pair policies retain independently enabled interactions, validate bounds and contradictory duplicates, and scale pair energy/force/virial. Coulomb model/list reconciliation checks overlaps and visits each special pair once. The correction outside the real cutoff consistently implements the documented bare-fraction difference convention; tests cover both sides of the cutoff.
- Required concrete/composite/backend cutoffs are exposed. Triclinic neighbor search uses all-pair enumeration and the complete cell participates in neighbor and reciprocal cache identity. Both audited skew fixtures are represented.
- Domain and composite force models honor accumulation; unsupported domain policies/box changes and GPU LJ policies are explicit errors. Focused caller inspection found force clearing before integration and replica recomputation; the manual domain integration fixture is updated accordingly.
- SPME owns plans, grids, worker scratch, and influence factors under an instance mutex; mesh changes rebuild the workspace and same-mesh cell/kappa changes refresh influence data. No global influence slot remains. Direct Ewald returns an immutable Arc cache snapshot. The allocation and concurrent-model tests exercise these contracts.
- Scope remains scalar isotropic pressure; tensor/constraint stress, runner guards, and production policy wiring appropriately remain later tasks. Not rewriting pair scratch without measured justification matches the conditional plan requirement.

## Evidence and verification limits

Read the supplied diff, Task 2 brief/report, accepted spec, and plan Global Constraints. The initial large diff output was truncated; only the missing ranges were then read in smaller chunks. Independently checked that HEAD is the reviewed commit and the tracked worktree was clean. Checked the unchanged reciprocal-vector generator specifically for skew-cell completeness: its integer bounds use direct-basis norms, so the full-cell cache does not conceal an incomplete enumeration. Checked integrator, Nose-Hoover, replica, and initial runner callsites specifically for accumulation integration risk; GPU host readback also adds into existing force/virial.

Test assertions and numerical expectations in the diff substantiate the reported scientific coverage. The allocation binary has one test and counts global allocations/reallocations around warmed complete passes, including the parallel FFT case. Cache correctness comparisons use independently constructed model instances. The reported default count arithmetic is consistent: 125 covering checks plus 4 water tests equals 129.

Per the review brief, no broad suite, Clippy matrix, or redundant test rerun was performed. Historical RED executions, reported GREEN process exits, optional actual GPU adapter execution, and the measured zero-allocation results cannot be independently certified from source/report alone; no raw execution logs were supplied in the task evidence directory. This is a provenance limit, not an identified implementation failure. No focused experiment was needed after the two concrete code-inspection risks were resolved.

## Strengths

- Independent analytical and finite-difference expectations check physics rather than simply reproducing implementation arithmetic.
- Reconciliation avoids allocating a merged policy per force pass and avoids applying duplicate Coulomb weights twice.
- Cache ownership and lock lifetimes make cross-simulation state isolation evident from the code, while preserving reusable storage.
- Backend restrictions are explicit, and the report distinguishes allocation evidence from an unmeasured throughput claim.

## Findings

### Critical

None.

### Important

None.

### Minor

1. `crates/modena-core/src/potential/tabulated.rs:4-6` still says the conservation contract is whatever the table says and makes the caller responsible for a self-consistent table. The new joint Hermite interpolation guarantees `dU/dr = -F` between nodes even when the supplied endpoint values are an imperfect approximation to an underlying potential. Update this header to distinguish the guaranteed interpolation derivative contract from the caller's responsibility for the physical accuracy of the input table. This is documentation only and does not block approval.

## Quality verdict: Approved

No blocking correctness, scientific-scope, or maintainability finding. The minor header correction can be made alongside later documentation cleanup. Review created only this ignored evidence file; no source, index, or HEAD changes.

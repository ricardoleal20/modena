# Task 1 scoped re-review — fix round 1

FIX_BASE `c6fef76ae58500bd545930224cda0980c607fed0` → HEAD `dd22f471a40781085c9eb375e0afcc1f7be239a9`.

Read the complete fix diff and the Fix round 1 implementation-report addendum. Scope was the two original findings and new breakage introduced by the fix. No suite rerun, source/index/HEAD modification, or broad re-review was performed.

## Original findings

### Important: systematic rigid-water rotational cooling — ADDRESSED

Dynamical SHAKE now supplies the saved pre-drift water positions to the coupled Newton solver. Its residuals use the evolving predicted distances; its Jacobian columns and position impulses use fixed pre-drift vectors. I checked the three Jacobian rows against the mass-weighted OH1/OH2/HH incidence relations, including their signs and distinct off-diagonal entries. They are correct for the generally nonsymmetric reference-direction solve. The existing Cramer solver handles a general matrix; its former symmetric label was merely documentation.

This fixes the specific mechanism identified in the first review. With the existing correction-only velocity update and final tangent projection, the reference-direction solve implements a RATTLE drift. Standalone initialization still selects current directions through None; it retains the previous geometric projection behavior. Reference slicing is covered by the existing atom-index and reference-buffer validation, and iteration-limit/singular-solve/final-residual failures remain errors.

The new regressions have independent physical oracles:

- A free equilateral rotor rotating about a principal axis has constant K=1.5. Tests exercise both BAOAB at gamma=0 and velocity Verlet at three timesteps, measuring maximum drift over the entire trajectory. Reported RED reproduces the independently observed cooling; GREEN maximum relative errors range from approximately 1e-11 to 1e-9, comfortably below the 1e-7 bound.
- Subtracting instantaneous COM translation leaves three nonlinear rotational modes, so K_rot has Gamma(3/2,T) statistics. Mean=variance=1.5 at T=1 is correct. The centered squared-energy estimator has canonical variance 13.5; the report and tests account for energy autocorrelation instead of treating consecutive samples as independent. Two friction/timestep cases and three seeds cover both weak and stronger thermal coupling. Reported weak-friction mean changes from 1.26588 before the fix to 1.52568 afterward; the second case gives 1.50568. Both post-fix variances also pass their stated correlation-aware limits.

The statistical tests support the scoped constrained-mode repair. They do not prove arbitrary molecular or orientation-distribution correctness, and the report explicitly avoids that claim.

### Minor: obsolete parallel layout comment — ADDRESSED

RigidWater documentation now describes selected contiguous triples with unrelated atoms permitted. It distinguishes initialization projection from reference-direction dynamical RATTLE and removes the incorrect parallel-chunk claim.

## New breakage in the fix diff

None found. The fix changes only the water position projection and its documentation, with dedicated regression coverage. It leaves public signatures, stochastic state, pressure calculations, and component parameters unchanged.

## Verification evidence

The implementation report records both new tests failing before the production fix, followed by 39 integration cases and four constraint unit cases passing, including the rigid-water/SPME trajectory. Default workspace Clippy, formatting, and diff checks are reported green. These results were evaluated against the code and independent physical expectations; I did not rerun them.

## Spec and task-quality status

**Spec: compliant for this scoped fix; both original findings addressed.**

**Task quality: Approved.** No remaining Critical, Important, or Minor findings in this fix diff.

The original Task 2/3 integration checks remain controller obligations, including explicit constrained-pressure limitations and constrained-NPT/box-rescaling rejection. They are unchanged by this fix and were not re-reviewed here.

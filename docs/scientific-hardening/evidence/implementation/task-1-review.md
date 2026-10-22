# Task 1 independent review

Reviewed immutable BASE `8078e1d70301d040d1113fc9125be09fa4cb7786` → HEAD `c6fef76ae58500bd545930224cda0980c607fed0`, the task brief, accepted spec, complete diff package, and implementation report. Read-only review except this report. The full core suite was not rerun; the reported 166-pass result, formatting, and Clippy results are implementation evidence, not independently rerun results.

## Spec compliance

**Issues found.** Most specified interfaces and audited repairs are delivered, but the supported constrained-water Langevin path still has systematic dissipation in its drift projection. The free-atom regression is necessary and well derived, but does not validate the constrained modes. This blocks approval of the broad constrained-sampling acceptance.

Compliant portions:

- CSVR uses the independently correct chi-square factor, signed transition, sum-of-squares evaluation, and special handling for one DOF. Conditional and stationary moments and three fixed seeds provide meaningful validation.
- Projection APIs return Result and callers propagate errors. Index, mass, geometry, tolerance, reference-buffer, singular-solve and final-residual checks are substantially improved.
- System exposes explicit constraints and conserved COM reductions with manual-system 3N defaults; Langevin retains its stochastic COM modes.
- Triclinic nearest-image search has a valid inverse-row bounding argument; tests include a translation beyond an adjacent-cell stencil. Scaling consistently rebuilds the basis and inverse.
- CSVR/Langevin checkpoints preserve exact named RNG state; NHC includes its extended positions and momenta. Kind/parameter comparisons, JSON float roundtripping and failure-before-mutation checks are appropriate.
- Replica exchange moves systems together with cached forces and neighbors, applies the destination/source temperature factors, tracks last-step PE, and rejects absent/invalid model cutoff and invalid configuration parameters. Its acceptance sign is correct; the report correctly notes that the old executable sign was already equivalent.
- The independently derived rigid-dimer pressure assessment supports the explicit missing-constraint-virial limitation. Deferring its general implementation and requiring constrained-NPT rejection is allowed by the brief.

## Findings

### Important — Rigid-water drift projection systematically cools rotational modes

Locations: `crates/modena-core/src/integrator.rs:208`, `:227`, `:274`; underlying water position solve `crates/modena-core/src/constraints.rs:378` and `:392`.

BAOAB calls `project_drift` twice per step, adding each position correction divided by its drift interval and then projecting velocities tangent. This is correct for a RATTLE drift when the position correction uses the constraint directions at the beginning of that drift. The water Newton solver instead repeatedly uses the current predicted geometry's directions. For a freely rotating equilateral rigid triple, this normalizes the predicted triangle radially, and the following tangent projection decreases angular velocity on every drift. There is no force or thermostat heat removal to justify that loss.

I compiled one standalone Rust probe against the existing built core rlib (no checkout changes and no suite rerun). Three unit masses form an equilateral triangle of circumradius 1, centered at (50,50,50) in a cubic box of side 100. All three constrained distances are sqrt(3). Initial velocities describe unit angular velocity about z. ForceModel returns zero, and Langevin has gamma=0, target T=1, seed=7. Initial K=1.5; exact free rigid dynamics keeps K constant. Running to time 10 gives:

| dt | Final kinetic energy | Relative drift |
|---|---:|---:|
| 0.01 | 1.4285714285016147 | -4.7619% |
| 0.005 | 1.463414634099744 | -2.4390% |
| 0.001 | 1.4925373149746837 | -0.4975% |

This is systematic first-order dissipation, not bounded second-order shadow-energy error. For this symmetric rotor, one projected drift of duration h changes omega to omega/sqrt(1+h^2*omega^2); after two drifts per BAOAB step, the observed decay is K(t)=K(0)/(1+t*dt/2) for the stated initial omega. A positive, weak Langevin friction must compete with this artificial cooling, so the unconstrained spectator's correct T/m variance cannot establish canonical constrained rotational sampling.

The current-normal water solve predates this task, but retaining it in the newly corrected constrained BAOAB path leaves the accepted sampling requirement unresolved. Implement the coupled water drift correction using the saved pre-drift constraint directions, or explicitly reject this water/integrator combination until validated. Add an independently expected rigid-water rotor conservation/timestep regression and a constrained-mode equilibrium statistic (rotational kinetic energy/equipartition at more than one friction and timestep). Keep the existing unrelated-free-atom regression.

### Minor — Water projection comments retain obsolete parallel-layout claims

Location: `crates/modena-core/src/constraints.rs:48` (RigidWater documentation).

The comment still describes contiguous triples as enabling one disjoint chunk per parallel ray, while the implementation now iterates selected triples sequentially and allows unrelated atoms. Update the comment to describe the actual indexing contract; this avoids preserving the previous all-atoms-are-waters assumption in documentation.

## Cannot verify in Task 1

These are explicit cross-task controller checks, not additional Task 1 failures:

- Task 2 must implement every concrete cutoff provider and verify actual non-2.5 cutoffs through replica execution. Its triclinic neighbor fallback must recover the audited skew-cell pair; correct min_image alone is insufficient.
- Task 3 must set DOF before initial velocity rescaling/output, project initial constraints, preserve System metadata on restart, and subtract COM only for dynamics preserving the initialized zero momentum. Check this with constrained and unconstrained initialization tests, including Langevin.
- Task 3 must serialize the physical identity alongside component state. VelocityVerlet uses the allowed null stateless default, so changed constraints/kind cannot be rejected by that hook alone; runner signature tests must reject these mismatches.
- Task 3 must reject invalid live dt, mass, box, thermostat parameters and unsupported zero-temperature/zero-DOF NHC combinations before stepping. The component snapshot finite checks do not replace execution validation. Singular/nonfinite triclinic cells also require rejection before nearest-image use.
- Task 3 must reject constrained NPT/box rescaling and clearly identify constrained pressure as incomplete. The rotor pressure test documents a known missing term; it does not validate pressure production.
- Reported RED chronology is consistent with the new tests and old code in the diff. The complete historical command logs were not independently replayed; a final diff cannot prove when tests were run.

## Strengths

The report distinguishes analytical expectations, actual tested properties, and deferred integration work. CSVR tests avoid the previous implementation-derived oracle. State restoration is narrow and mostly mutation-safe. Mechanical replica swaps preserve force/cache identity, and the acceptance-sign explanation avoids an unnecessary regression. The missing pressure term is demonstrated with an independently interpretable rotor instead of being hidden behind an invented stress correction.

## Task quality

**Needs fixes.** Resolve the Important rigid-water drift/sampling issue and its regression coverage before treating Task 1 constrained sampling as complete. No Critical findings.

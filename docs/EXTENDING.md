# Extending MODENA

The Rust library lives in `crates/modena-core`. Begin with an existing
implementation close to the method you want to add, and read its tests.
Public APIs may change during early development.

## A pair potential

Implement [`PairPotential`](../crates/modena-core/src/potential/mod.rs).
`energy_force(r2, ti, tj)` returns pair energy and `F_radial / r`. The pair
vector is `r_i - r_j`, so the force on atom `i` is that vector multiplied by
the returned scalar. `ti` and `tj` are species indices. Supply the cutoff,
its square, and the number of supported species when it exceeds one.

Use [`PairForceModel`](../crates/modena-core/src/force_model.rs) to evaluate
the potential over the CPU neighbor list. Check `dU/dr = -F_radial`, cutoff
behavior, species symmetry, and a short NVE trajectory. Tabulated models must
derive energy and force from the same interpolation; the existing implementation
uses joint cubic Hermite interpolation.

## A general force model

Implement [`ForceModel`](../crates/modena-core/src/force_model.rs) for a model
that does not fit pairwise evaluation. `compute()` must:

1. Read the current positions, box, species, and other relevant state.
2. **Accumulate** into `system.forces` and `virial.w`; the caller clears them
   once before a complete force evaluation.
3. Return this model's potential-energy contribution.

Declare `neighbor_cutoff()` when a driver needs a neighbor list.
`conserves_total_momentum()` defaults to false; opt in only when the model
preserves total momentum for the supplied system. This affects the eligible
center-of-mass kinetic modes. Finite charged SPME is an example where the
discretized model does not preserve continuous translation symmetry.

[`CompositeForceModel`](../crates/modena-core/src/composite.rs) combines model
contributions. Invalidate cached neighbor or reciprocal data when positions,
box geometry, or model state require it. Forces, energy, and virial must describe
the same state. Test composition against separately evaluated terms. The scalar
virial must account for explicit box dependence; an atom-based `r·F` sum is
insufficient for reciprocal electrostatics.

## Integrators and configuration

[`Integrator::step`](../crates/modena-core/src/integrator.rs) starts with
`system.forces` holding forces at the current positions and leaves forces at
the advanced state. Preserve this cache contract and report observables for
that state. Constraint handling and stochastic ensemble sampling need their
own derivation and tests; a final temperature near the target is insufficient.

If a method has stochastic or extended state, include it in checkpoint
continuation and physical identity checks. Compare a split run with an
uninterrupted run. Declare unsupported combinations explicitly.

A library implementation is not automatically available from TOML. To expose
one, add its typed input to
[`io/input.rs`](../crates/modena-core/src/io/input.rs), reject invalid and
unsupported combinations, then wire construction in
[`runner.rs`](../crates/modena-core/src/runner.rs). Add a small end-to-end input
test. Update Python and MCP schemas and capability descriptions where relevant.

MPI and GPU implementations have additional restrictions and examples in
`crates/modena-core/examples`. Their Cargo features enable library code; the
ordinary TOML runner uses CPU/Rayon.

Document evidence and remaining limitations in
[VALIDATED_SCOPE.md](VALIDATED_SCOPE.md), and follow
[CONTRIBUTING.md](../CONTRIBUTING.md) for verification and pull requests.

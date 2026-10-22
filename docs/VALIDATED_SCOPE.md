# Validation scope and known limitations

MODENA is an early development engine. This document describes implemented
methods, their focused regression coverage, and the limits of the release.
A passing test or a completed run does not certify arbitrary inputs or
converged material properties.

## Methods with focused coverage

| Area | Implementation and evidence |
| --- | --- |
| Pair interactions | Lennard-Jones, Buckingham, Morse, soft-sphere, and tabulated models. Force/energy and cutoff checks cover selected parameters. Tabulated energy and force share cubic Hermite interpolation, including the lower endpoint tangent extension; the legacy selector is still named `linear`. |
| NVE | Velocity-Verlet, with short energy/momentum tests and a timestep convergence check. Pair and direct-Ewald paths support orthorhombic and triclinic periodic cells. |
| NVT | Stochastic velocity rescaling, Langevin BAOAB, and Nosé–Hoover chains. Coverage includes conditional/stationary kinetic moments, free-atom Maxwell variance, selected constrained sampling, and extended-Hamiltonian checks. |
| Molecular models | Bond/angle/dihedral terms, explicit topology weights and exceptions, typed XYZ, SHAKE/RATTLE, and rigid-water projection. Force gradients, constraint errors, and selected rigid-rotation statistics are checked. |
| Electrostatics | Direct Ewald and orthorhombic SPME. Checks include force/energy gradients, strain virials, background/exclusion terms, and mesh refinement against direct Ewald for selected charged configurations. |
| Pressure coupling | Unconstrained isotropic Monte Carlo NPT. An ideal-gas volume distribution and accepted/rejected move cache integrity are tested. |
| Restart | Checkpoints include stochastic and extended state and physical identity. Tests compare uninterrupted and split runs for selected integrators/controllers and reject incompatible continuation. |

Start with [the argon examples](../examples/) and
[chemical workflows](../examples/chemical/README.md).
Representative tests are
[`sampling_hardening.rs`](../crates/modena-core/tests/sampling_hardening.rs),
[`rigid_sampling.rs`](../crates/modena-core/tests/rigid_sampling.rs),
[`force_hardening.rs`](../crates/modena-core/tests/force_hardening.rs),
[`barostat_hardening.rs`](../crates/modena-core/tests/barostat_hardening.rs), and
[`runner_hardening.rs`](../crates/modena-core/tests/runner_hardening.rs).

## Limits that affect interpretation

- **Constrained pressure is diagnostic.** It omits constraint impulse stress.
  Constrained NPT and full anisotropic coupling are rejected. Berendsen coupling
  is for equilibration; the legacy piston is experimental.
- **Liquid and chemical validation is incomplete.** Fixed-volume water and
  molecular examples exercise code paths, not converged density, structure,
  transport, or pressure. Broad chemical parameter sets and corrected real-unit
  ionic/water comparisons with a reference engine remain unvalidated.
- **Mesh accuracy depends on resolution.** SPME supports orthorhombic cells.
  A finite charged mesh can exert net force under a common translation; it
  retains all center-of-mass kinetic modes. Pair/direct-Ewald dynamics remove
  three modes only when the protocol conserves initialized zero momentum;
  Langevin retains them. Compare degrees of freedom before comparing temperature.
- **Replica exchange needs broader validation.** Swap and force-cache tests
  cover mechanics, not general convergence or equilibrium distributions.
- **MPI and GPU are restricted experimental backends.** MPI pair decomposition
  rejects special-pair policies and box changes. GPU uses f32 LJ offload with
  host readback. The ordinary TOML runner uses CPU/Rayon; features do not select
  a distributed or GPU simulation.
- **No general performance claim follows from the tests.** GPU tests may skip
  without an adapter. MPI needs a real multirank launch. A successful build or
  skip provides no hardware accuracy or scaling evidence.

NVE drift is one numerical diagnostic. For thermostatted or pressure-coupled
runs, system total energy changes through energy exchange; it is not the
conserved quantity. Small drift cannot establish correct ensemble sampling,
a correct force law, or converged observables.

## Units and restart compatibility

Real inputs use Å, ps, kcal/mol, amu, K, elementary charge, and atm. Thermodynamic
output uses the selected units; checkpoint mechanics and conservation energies
remain reduced. Tabulated input is supplied in reduced units.

The corrected real charge normalization is
`q* = q_e sqrt(332.06371 / (sigma_Å epsilon_kcal/mol))`.
The earlier factor changed Coulomb energy by sigma². Old real electrostatic
results with sigma different from one describe a different Hamiltonian.
Normal continuation checks identity and rejects incompatible checkpoints;
state-only import is an explicit protocol change, not migration of old physics.
See [unit normalization](scientific-hardening/UNITS.md).

## Evidence and reproduction

Run `cargo test --locked --workspace` for current regression checks. For the
bounded numerical harness:

```bash
cargo build --locked --release --bin modena
cargo build --locked --release -p modena-core --example scientific_validation
python3 bench/validate.py --output target/scientific-smoke
```

The harness needs Python 3.11 or newer. Its output records expected values,
tolerances, and independent bases. Correlated samples and bounded checks are
not formal uncertainty estimates.

[Scientific results](scientific-hardening/RESULTS.md) and
[verification records](scientific-hardening/VERIFY.md) preserve earlier
measurements with commit and environment provenance. Their passing results
belong to the recorded revisions. New release checks must be recorded separately;
CI configuration alone is not evidence that remote CI has run.

The workspace declares Rust 1.90 to cover the locked optional dependency graph.
The default and GPU targets passed local checks with Rust 1.90 for version
0.1.1; the numerical regression suite used Rust 1.97.1. CI also checks the
minimum compiler. See the [release verification](releases/v0.1.1-verification.md).
A build on one platform does not establish support for every target or backend.

## Historical material

[Benchmark methodology](../bench/BENCHMARKS.md) and the scientific records link
raw data and earlier reports. Timings, error bounds, and explanations apply to
their recorded inputs, resources, accuracy settings, and revisions. They are
not re-established by this documentation update.

[`ARCHITECTURE.md`](../ARCHITECTURE.md), [`ROADMAP.md`](../ROADMAP.md), and older
planning notes include design intent and past status. Use this document, the
current source, and the tests for release scope.

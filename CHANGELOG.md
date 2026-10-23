# Changelog

## Unreleased

- Refresh Cargo dependencies, workflow actions, and the optional MCP SDK;
  adapt benchmark and GPU APIs while retaining the Rust 1.90 compiler floor.
- Group compatible Dependabot updates and check both Cargo workspaces and the
  optional Python tooling in CI.
- Complete contribution and commit guidance, issue and pull request templates,
  conduct, security, and maintainer policies.
- Document a proposed design system and scientific interface guidance.
- Include MIT notices and discovery metadata in distributable packages.
- Correct MCP response field names to the protocol's camelCase wire format.
- Reject overflowing structure counts and invalid cost estimates with recoverable
  MCP errors; reserve queued and running capacity atomically and validate entire
  sweeps before launching jobs.
- Rust callers of `RunManager::estimate_run_cost` now handle a `Result` so invalid
  input is reported explicitly.

## 0.1.1 — 2026-10-22 (prerelease)

### Numerical methods and workflows

- Correct stochastic velocity rescaling and constrained Langevin sampling;
  track kinetic degrees of freedom and report constraint projection failures.
- Refresh force, energy, and virial caches after changes of box or replica state.
- Interpolate tabulated energy and force with a shared Hermite polynomial.
- Use strain derivatives for Ewald/SPME virials and handle triclinic pair and
  direct-Ewald geometry.
- Add unconstrained isotropic Monte Carlo NPT with rejected-move rollback.
  Reject constrained NPT and full anisotropic coupling.
- Preserve stochastic and extended ensemble state in checkpoints; validate
  physical identity during continuation and provide explicit state-only import.
- Correct real-unit charge normalization. Earlier real electrostatic runs with
  a reference length other than one use a different Hamiltonian and are not
  interchangeable through normal continuation.
- Add typed molecular input, separate LJ/Coulomb topology weights, stricter input
  checks, and consistent Python/MCP units and metadata.
- Add focused numerical regression tests and benchmark provenance tools.
- Add independent CSVR kinetic-energy and signed-velocity regression checks.
- Respect custom Cargo target directories in the Python client, clean up
  temporary inputs on errors, and preserve CLI failures before parsing output.

### Documentation and build metadata

- Replace broad completion claims with an early development scope and explicit
  limitations.
- Document the CLI quickstart, contribution process, support, and Rust extension
  contracts; add software citation metadata.
- Use the public repository URL in Cargo metadata.
- Raise the declared Rust minimum to 1.90 to cover resolved optional GPU
  dependencies, with a CI build check for that compiler version.
- Select the MPI test executable from Cargo's artifact messages in CI.

## 0.1.0

Initial development version. Earlier reports and benchmark measurements remain
in the repository with their original provenance; they do not validate later
versions.

# Contributing

Bug reports and pull requests are welcome through the
[GitHub issue tracker](https://github.com/ricardoleal20/modena/issues) and
[pull requests](https://github.com/ricardoleal20/modena/pulls).
For a substantial change, open an issue first to explain the intended behavior
and how it can be checked. The APIs and input format may change during early
development.

## Reporting a problem

Include the MODENA version or commit, operating system, Rust version, command,
and a small input file that reproduces the problem. For a numerical discrepancy,
include units, timestep, cutoff, ensemble settings, initial conditions, and the
expected result with its source. Small logs and the exact error are more useful
than a full trajectory.

## Making a change

Keep each pull request focused. Explain the behavior it changes and how you
checked it. Run the relevant checks from the repository root:

```bash
cargo fmt --all -- --check
cargo clippy --locked --workspace --all-targets -- -D warnings
cargo test --locked --workspace
```

For changes to physics, add a check against an independently derived result or
reference calculation. Force/energy gradients, conservation, and ensemble
statistics test different properties; use the ones relevant to the change.
State tolerances, units, and any sampling limits. A test that copies the same
formula as the implementation is insufficient evidence for a physics claim.

Update [validation scope](docs/VALIDATED_SCOPE.md) when coverage or limitations
change. Keep performance claims tied to a reproducible command, hardware,
compiler, system size, and accuracy settings. Existing benchmark reports are
historical records.

Optional backends need separate checks. GPU tests may skip without an adapter;
MPI tests need a system MPI installation and multiple ranks. Describe which
paths actually ran. See [extending the engine](docs/EXTENDING.md) for the force
accumulation contract and configuration wiring.

Contributions are made under the repository's [MIT license](LICENSE).

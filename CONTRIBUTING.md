# Contributing

Thank you for helping improve MODENA. Documentation, small reproducible bug
reports, numerical checks, and focused code changes are welcome. Read the
[code of conduct](CODE_OF_CONDUCT.md) and [validation scope](docs/VALIDATED_SCOPE.md)
before contributing. The APIs and input format may change during early development.

Use the [issue tracker](https://github.com/ricardoleal20/modena/issues) for bugs,
numerical discrepancies, and proposals. For a substantial change, discuss the
intended behavior and how to check it in an issue before investing in an
implementation. See [support](SUPPORT.md) for questions and [security](SECURITY.md)
for sensitive reports.

## Reporting a problem

Include the MODENA version or commit, operating system, Rust version, command,
and a small input file that reproduces the problem. For a numerical discrepancy,
include units, timestep, cutoff, ensemble settings, initial conditions, and the
expected result with its source. State the observed values, comparison tolerance,
and whether the result depends on a seed, thread count, or optional backend.
Small logs and the exact error are more useful than a full trajectory. Remove
credentials and private data before posting.

## Development setup and checks

Fork the [repository](https://github.com/ricardoleal20/modena), clone your fork,
and create a branch for your change. Install Rust 1.90 or newer with `rustfmt`
and `clippy`. The default build needs no system MPI or GPU runtime. Run these
checks from the repository root:

```bash
cargo fmt --all -- --check
cargo clippy --locked --workspace --all-targets -- -D warnings
cargo test --locked --workspace
```

Python tooling requires Python 3.11 or newer. The regression suites use the
standard library and need a built CLI for integration checks:

```bash
cargo build --locked --release --bin modena
MODENA_BINARY="$PWD/target/release/modena" python3 -m unittest discover -s python/modena/tests -v
python3 -m unittest discover -s bench/tests -v
```

Set `MODENA_BINARY` to your executable if you use a custom Cargo target directory
or a platform-specific executable suffix. For changes to numerical behavior,
also run the bounded scientific harness:

```bash
cargo build --locked --release -p modena-core --example scientific_validation
python3 bench/validate.py --output target/scientific-smoke
```

Use the harness's `--binary` option if the scientific example is in a custom
target directory or has a platform-specific executable suffix.

For GPU changes, run:

```bash
cargo test --locked -p modena-core --features gpu --test gpu_pair -- --nocapture
```

GPU tests can return successfully after skipping when no suitable adapter is
available. Report whether the hardware path actually ran. MPI changes need a
system MPI installation and a real multirank launch of the `mpi_decomp` test;
a single-process run skips the distributed check. Follow the artifact-selection
and two-rank launch in [CI](.github/workflows/ci.yml). Enabling either feature
does not select that backend in the ordinary TOML runner.

## Numerical and performance evidence

For changes to physics, add a check against an independently derived result or
reference calculation. Force/energy gradients, conservation, and ensemble
statistics test different properties; use the ones relevant to the change.
State tolerances, units, and any sampling limits. A test that copies the same
formula as the implementation is insufficient evidence for a physics claim.
For statistical comparisons, describe equilibration, sample counts, correlation,
and uncertainty assumptions. A short run or small NVE drift alone cannot
establish correct ensemble sampling or converged material properties.

Update [validation scope](docs/VALIDATED_SCOPE.md) when coverage or limitations
change. Keep performance claims tied to a reproducible command, hardware,
compiler, system size, and accuracy settings. Existing benchmark reports are
historical records.

## Pull requests and review

Keep each pull request focused. Explain the problem, resulting behavior, and
checks you ran, including skipped or unavailable checks. Add regression coverage
when behavior changes and update affected examples or documentation. See
[extending the engine](docs/EXTENDING.md) for the force accumulation contract
and configuration wiring.

Use [Conventional Commits](COMMIT_RULES.md) for new commits and pull request
titles. Existing history does not need to be rewritten. The maintainer reviews
changes and may request revisions; response times depend on availability.
See [governance](GOVERNANCE.md) for decisions and release responsibility.

By submitting a contribution, you agree to license it under the repository's
[MIT license](LICENSE). Keep third-party attribution and license notices intact.
No contributor license agreement or sign-off is required.

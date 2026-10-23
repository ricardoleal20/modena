# MODENA

[![CI](https://github.com/ricardoleal20/modena/actions/workflows/ci.yml/badge.svg)](https://github.com/ricardoleal20/modena/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-087F78.svg)](LICENSE)

MODENA (Molecular Dynamics Engine for Numerical Atomistics) is a molecular
dynamics engine written in Rust. It includes a library and a command-line
program that reads TOML input and writes XYZ trajectories and CSV thermodynamic
logs.

This is an early development release for classical pair, ionic, and small
molecular simulations. It includes bonded interactions, rigid constraints,
Ewald/SPME electrostatics, and NVE/NVT integration. Numerical checks cover
selected models and workflows; they do not establish converged liquid
properties or general agreement with a reference engine.
See [validation scope and known limitations](docs/VALIDATED_SCOPE.md).

## Build and run

Install Rust 1.90 or newer, clone the repository, and build the CLI:

```bash
git clone https://github.com/ricardoleal20/modena.git
cd modena
cargo build --locked --release --bin modena
cargo run --locked --release --bin modena -- check -i examples/argon_nve.toml
cargo run --locked --release --bin modena -- run -i examples/argon_nve.toml
```

The example runs 256 Lennard-Jones particles for 2,000 NVE steps in reduced
units. It writes `examples/argon_nve_traj.xyz` and
`examples/argon_nve_energy.csv`; output paths are relative to the input file.
The CLI reports final energies, temperature, pressure, and NVE energy drift.
Drift is one diagnostic of a simulation's numerical behavior.

To install the CLI from a source checkout onto your Cargo binary path:

```bash
cargo install --locked --path crates/modena-cli
modena --help
```

[argon_nvt.toml](examples/argon_nvt.toml) gives a Langevin example.
[Chemical workflow examples](examples/chemical/README.md) cover molecular
input, restart, ions, isotropic Monte Carlo NPT, and fixed-volume water.
These short runs are examples, not validated material parameter sets.
The default build needs no MPI installation or GPU runtime.

## Development

```bash
cargo test --locked --workspace
cargo fmt --all -- --check
cargo clippy --locked --workspace --all-targets -- -D warnings
```

- [Contributing](CONTRIBUTING.md): reporting bugs and proposing changes.
- [Extending the engine](docs/EXTENDING.md): Rust traits and input wiring.
- [Validation scope](docs/VALIDATED_SCOPE.md): evidence, limitations, and historical reports.
- [Changelog](CHANGELOG.md) and [support](SUPPORT.md).
- [Commit rules](COMMIT_RULES.md), [code of conduct](CODE_OF_CONDUCT.md), and
  [governance](GOVERNANCE.md).
- [Security policy](SECURITY.md): reporting vulnerabilities and running untrusted inputs.
- [Design system](DESIGN.md): visual guidance for future interfaces.
- [Maintainer guide](docs/MAINTAINING.md): release checks and repository settings.

The [Python package](python/README.md) configures runs through a CLI subprocess.
`modena-mcp` provides an additional MCP interface. MPI and GPU backends are
optional library code with separate examples; enabling their Cargo features
does not select them in the TOML runner.

## License and citation

MODENA is licensed under the [MIT license](LICENSE). To cite the software, use
[CITATION.cff](CITATION.cff) and include the version or commit used in your work.

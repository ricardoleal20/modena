# modena-cli

The command-line interface to **MODENA**, a molecular dynamics engine in Rust.
The executable is named `modena`. It reads TOML configurations and writes XYZ
trajectories and CSV thermodynamic logs.

From a [source checkout](https://github.com/ricardoleal20/modena), with Rust 1.90
or newer:

```sh
cargo install --locked --path crates/modena-cli
modena --help
modena check -i examples/argon_nve.toml
modena run -i examples/argon_nve.toml
```

Example paths refer to files in that checkout. Output paths are relative to
the input file. See the
[quickstart](https://github.com/ricardoleal20/modena#build-and-run) and
[validation scope](https://github.com/ricardoleal20/modena/blob/main/docs/VALIDATED_SCOPE.md)
before interpreting results. MODENA is in early development.

Licensed under MIT. The complete notice is included in `LICENSE`.

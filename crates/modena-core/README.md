# modena-core

Rust library for **MODENA**, the Molecular Dynamics Engine for Numerical
Atomistics. This is an early development library for classical pair, ionic,
and small molecular simulations. Public APIs and input formats may change.

The default build requires Rust 1.90 or newer and no MPI or GPU runtime.
Optional `mpi` and `gpu` features provide library backends with separate
system requirements. Feature flags alone do not select a backend in TOML runs.

See the [repository quickstart](https://github.com/ricardoleal20/modena#build-and-run),
[extension guide](https://github.com/ricardoleal20/modena/blob/main/docs/EXTENDING.md),
and [validation scope](https://github.com/ricardoleal20/modena/blob/main/docs/VALIDATED_SCOPE.md).
Numerical checks cover selected models and workflows; they do not establish
general material accuracy or converged liquid properties.

Licensed under MIT. The complete notice is included in `LICENSE`.

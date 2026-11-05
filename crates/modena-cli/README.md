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
modena run --diagnostics-json -i examples/argon_nve.toml
```

Example paths refer to files in that checkout. Output paths are relative to
the input file. See the
[quickstart](https://github.com/ricardoleal20/modena#build-and-run) and
[validation scope](https://github.com/ricardoleal20/modena/blob/main/docs/VALIDATED_SCOPE.md)
before interpreting results. MODENA is in early development.

Energy diagnostics include every integration step in reduced units. NVE
advisories require explicit `run.energy_error_thresholds` and otherwise remain
`Uncalibrated`; NVT/NPT have no NVE verdict. The optional JSON flag exposes the
whole report and a separate `run metadata = {JSON}` line (core ensemble, DOF,
pressure convention and optional provider identity). The energy diagnostics line
retains its existing format. See the [B06 metric definitions](https://github.com/ricardoleal20/modena/blob/main/docs/B06_ENERGY_DIAGNOSTICS.md).

## Native LJ snapshot

`modena snapshot` reads exactly one strict JSON request from stdin and writes
one JSON response to stdout. It computes energy and forces through the existing
native shifted-force LJ pair model, without integration, velocities, simulation
configuration or output files. Errors have nonzero status and diagnostics on
stderr, without a partial response on stdout.

```sh
modena snapshot <<'JSON'
{
  "protocol": "modena-native-lj-v1",
  "units": {"length": "angstrom", "energy": "eV", "force": "eV/angstrom"},
  "atomic_numbers": [18, 18],
  "positions_angstrom": [[0, 0, 0], [4, 0, 0]],
  "cell_angstrom": [[20, 0, 0], [0, 20, 0], [0, 0, 20]],
  "pbc": [true, true, true],
  "sigma_angstrom": 3.4,
  "epsilon_ev": 0.0103,
  "cutoff_angstrom": 8.5,
  "skin_angstrom": 0.3
}
JSON
```

The response has exactly `protocol`, `units`, `energy_ev` and
`forces_ev_per_angstrom` (N rows of three values). Units/protocol match the
request; no stress property is exposed. The example returns energy about
-0.00897848103 eV and first-atom x force about +0.00560861341 eV/Angstrom.
The pair convention is U_sf(r) = U(r) - U(rc) + (r-rc) F(rc), with
F(r) = -dU/dr and zero interaction at/above rc.

All fields are required; unknown fields and trailing JSON are rejected.
Only homogeneous atomic numbers 1..118, finite N×3 positions, full PBC and
positive finite diagonal orthorhombic cells are supported. Sigma/cutoff must be
positive, epsilon/skin nonnegative, all finite. Cutoff plus skin must be at
most half each box length. Minimum-image neighbor pairs with squared distance
below 1e-14 Angstrom² are rejected rather than silently skipped. Input is
limited to 8 MiB, 100000 atoms and a dense native neighbor grid of 1000000 cells;
overflow/nonfinite evaluation fails explicitly. There are no optional ASE or
NumPy dependencies for this command. The optional
[`modena.ase.ModenaCalculator`](../../python/README.md#modena-as-an-optional-ase-calculator)
consumes this protocol.

## Optional ASE EMT provider

Install the optional calculator in a local environment and discover its exact pins:

```sh
python3 -m venv .venv-ase
.venv-ase/bin/python -m pip install 'ase>=3.26,<4'
modena provider-info --python "$PWD/.venv-ase/bin/python" --timeout-ms 30000
```

`provider-info` writes one strict core JSON descriptor: installed ASE/NumPy
versions, EMT source SHA-256, bridge protocol, task, composition, units and
capabilities. Copy its exact `ase_version` and `model_sha256` into `[model]`
with `type="ase_emt"`, `composition="total"`, explicit per-atom
`atomic_numbers=[29,29,29,29]` for a one-cell FCC Cu structure, and the selected
`python` executable. Set `units="real"` and structure `mass=63.546` explicitly.
No chemical numbers or masses are inferred from cosmetic element labels.
The [Python example](../../python/README.md#optional-ase-emt) renders a complete
input using discovered pins.

This first provider supports seven standard metals (Al, Ni, Cu, Pd, Ag, Pt, Au),
fully periodic fixed-cell unconstrained NVE/NVT and total energy/forces/stress.
PBC is exactly `[true,true,true]`; only `asap_cutoff=false` is supported.
Charges, topology/electrostatics additions, constraints and barostats are rejected.
The provider manages neighbors and starts a fresh Python process for every
snapshot; this is a parity workflow with startup latency, not a throughput or
material-validation claim. `timeout_ms` defaults to 30000 per evaluation.
Human run output names the installed provider/version; `--diagnostics-json`
includes its full core identity in metadata. Native runs need no Python/ASE.

Licensed under MIT. The complete notice is included in `LICENSE`.

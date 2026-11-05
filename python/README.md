# MODENA Python wrapper

The `modena` Python package provides typed configuration objects and runs the
MODENA molecular dynamics CLI through a subprocess. The base package needs
Python 3.11 or newer and has no runtime dependencies. The Rust engine is a
separate executable and is not included in the Python distribution.

## MODENA as an optional ASE calculator

Install the optional dependency with `python -m pip install './python[ase]'`
from a source checkout, and build/install the matching CLI. Import this adapter
explicitly; base `import modena` does not import ASE or NumPy.

```python
from ase import Atoms
from modena.ase import ModenaCalculator

atoms = Atoms("Ar2", positions=[[0, 0, 0], [4, 0, 0]], cell=[20, 20, 20], pbc=True)
atoms.calc = ModenaCalculator(
    sigma_angstrom=3.4, epsilon_ev=0.0103, cutoff_angstrom=8.5,
    skin_angstrom=0.3, binary="/absolute/path/to/modena", timeout=30.0,
)
print(atoms.get_potential_energy())  # eV; about -0.00897848103
print(atoms.get_forces())           # eV/Angstrom; first x about +0.00560861341
atoms.calc.set(epsilon_ev=0.02)     # invalidates cached properties
```

This evaluates the native **shifted-force** LJ model using one CLI process per
snapshot. Both potential and radial force go to zero at the cutoff. It supports
homogeneous atomic numbers 1..118 (one species), full periodicity, finite
positions and positive diagonal orthorhombic cells. Atomic numbers specify
identity; explicit sigma/epsilon specify the interaction without a parameter
database. Sigma/cutoff must be positive, epsilon/skin nonnegative and all values
finite. Cutoff plus skin must fit within half every cell length. Constraints,
nonzero initial charges or magnetic moments, mixed species, skew/rotated cells
and partial PBC are rejected. Coincident or singular minimum-image pairs
(squared distance below 1e-14 Angstrom² in the neighbor list) are rejected.

Only `energy` and `forces` are implemented. Stress requests raise ASE's
`PropertyNotImplementedError`; no tensor stress is inferred from scalar virial.
Atom positions/numbers/cell/PBC and calculator parameter changes invalidate the
cache, including position changes below ASE's usual tolerance. Constraints,
charge and spin are checked even on cached accesses. Failures
clear results. Use `calc.set(...)` for changes to model parameters, `binary`
or `timeout` (seconds, positive and finite). The executable is a single path or
PATH name, passed without a shell; spaces within an executable path work.
With `binary=None`, the existing wrapper's PATH/local build discovery is used.

The [snapshot CLI protocol](../crates/modena-cli/README.md#native-lj-snapshot)
caps input at 8 MiB, atom count at 100000, the grid at 1000000 cells, and
retained pairs within cutoff plus skin at 1000000 (up to 16 MB of requested
pair storage on a 64-bit host). Sparse inputs are limited by their actual
neighbor count. This initial adapter is a scientific interchange path; it does not
establish bulk material accuracy or claim efficient large-scale ASE dynamics.

MODENA is in early development. See the
[validated scope and known limitations](https://github.com/ricardoleal20/modena/blob/main/docs/VALIDATED_SCOPE.md)
before interpreting simulation results. Short examples demonstrate workflows;
they do not establish converged material properties.

## Install from a source checkout

Install Rust 1.90 or newer, Python 3.11 or newer, and Git. From a checkout:

```bash
git clone https://github.com/ricardoleal20/modena.git
cd modena
cargo install --locked --path crates/modena-cli
python3 -m venv .venv
. .venv/bin/activate
python -m pip install ./python
```

Make sure Cargo's binary directory (usually `~/.cargo/bin`) is on `PATH` so
the wrapper can find `modena`. On Windows, activate the environment with
`.venv\Scripts\Activate.ps1` in PowerShell. These instructions install the
wrapper and CLI from this repository; they do not depend on a PyPI release.

For wrapper development, replace the last command with
`python -m pip install -e ./python`. To keep the Rust binary in the checkout,
use `cargo build --locked --release --bin modena` and pass its absolute path
to `check()` and `run()` with `binary="/path/to/modena/target/release/modena"`
(`modena.exe` on Windows). When importing the package directly from this source
tree, automatic discovery also checks release/debug build directories, including
`CARGO_TARGET_DIR`; installing the CLI on `PATH` or passing `binary=` works
regardless of where the wrapper is installed.

## Configure, validate, and run

```python
from modena import Config, Fcc, LennardJones, Run, VelocityVerlet, check, run

cfg = Config(
    structure=Fcc(lattice_constant=1.6796, replicas=[2, 2, 2]),
    model=LennardJones(sigma=1.0, epsilon=1.0, cutoff=2.5),
    integrator=VelocityVerlet(dt=0.002),
    run=Run(steps=100, initial_temperature=0.8, seed=42),
)

print(cfg.to_toml())
print(check(cfg))
result = run(cfg)
print(result.steps_rows[-1])
print(result.conservation)
print(result.energy_diagnostics)
```

This example uses Lennard-Jones reduced units. `check()` asks the engine to
validate the rendered configuration; it does not certify the physical model.
`run()` returns `RunResult` with thermodynamic rows as a list of dictionaries,
the scoped conservation verdict, an optional timing string, the rendered TOML
and `energy_diagnostics` containing full-rate reduced-energy/time statistics.
The appended optional `metadata` dictionary contains the core ensemble, DOF,
pressure convention and `provider` descriptor (`None` for native models).
Older CLIs without the metadata line return `metadata=None`.
NVE is `Uncalibrated` unless `Run(energy_error_thresholds={"warning": ...,
"drifting": ...})` selects case-specific tolerances. NVT/NPT report
`NotApplicable`; cold NVE has an unavailable normalization. The verdict does
not certify accuracy or sampling. Older CLIs remain supported with a legacy
verdict and `energy_diagnostics=None`. See the repository's
[B06 diagnostic contract](https://github.com/ricardoleal20/modena/blob/main/docs/B06_ENERGY_DIAGNOSTICS.md).

Input, trajectory, and checkpoint paths resolve relative to the current
directory, or to `base_dir=` when supplied to `check()` or `run()`. The wrapper
collects thermodynamic data through a temporary CSV and leaves a configured
`run.log_file` unchanged. Explicit trajectory or checkpoint outputs can still
write files. Both operations raise `RuntimeError` when the CLI reports failure;
`check()` has a 120-second subprocess timeout and `run()` has a one-hour timeout.

The optional `mcp` extra declares an MCP client dependency. The public
`check()` and `run()` functions use the CLI subprocess; the extra does not
switch their transport. The separate Rust MCP executable and server interface
are described in the
[main repository](https://github.com/ricardoleal20/modena).

## Optional ASE EMT

The default package remains dependency-free. Install the extra into your local
Python environment with `python -m pip install './python[ase]'` from the checkout;
it selects `ase>=3.26,<4`. Discover exact installed pins for each run:

```python
import sys
from pathlib import Path
from modena import AseEmt, Config, Fcc, Run, VelocityVerlet, check, provider_info, run

# Keep the venv executable path: resolving its symlink can select base Python.
python = sys.executable
identity = provider_info(python, timeout_ms=30000)
cfg = Config(
    structure=Fcc(3.61, [1, 1, 1], element="Cu", mass=63.546),
    model=AseEmt(
        ase_version=identity["ase_version"],
        model_sha256=identity["model_sha256"],
        atomic_numbers=[29, 29, 29, 29], composition="total", python=python,
    ),
    integrator=VelocityVerlet(0.0001),
    run=Run(steps=2, initial_temperature=10.0, log_every=1), units="real",
)
Path("cu_emt.toml").write_text(cfg.to_toml())
print(check(cfg))
result = run(cfg)
assert result.metadata["provider"] == identity
```

All three functions accept `binary=` for an explicit CLI. `provider_info()`
returns the unmodified core JSON from `modena provider-info`; exact version and
source-hash pins, explicit per-atom Z and `composition="total"` are required.
The hash identifies installed EMT source/model data, not learned weights.
`AseEmt` defaults to `python="python3"`, `pbc=[True,True,True]`,
`asap_cutoff=False`, `timeout_ms=30000` and `transport="one_shot"`. Select
`AseEmt(..., transport="persistent")` to reuse one Python process and EMT
calculator per provider. Transport does not change physical identity or
compatible checkpoint signatures. The engine enforces seven standard
metals (Al, Ni, Cu, Pd, Ag, Pt, Au), real units, full PBC, unconstrained fixed-cell
NVE/NVT and total composition. Charges, topology/electrostatics additions,
constraints and all barostats are rejected. Masses remain explicit mechanical
inputs; labels and force-field types never infer Z or mass.

`check()` and `run()` resolve a `model.python` path containing a separator
against `base_dir=` (default cwd), preserving virtual-environment symlinks.
Bare executable names remain PATH lookups. Other input/output path behavior is
unchanged. One-shot launches a fresh Python process per snapshot; persistent
reuses a session. The provider manages its own neighbors. Persistent failures
leave the session terminal; drop/failure kills and reaps the direct child.
Cancellation is checked at safe setup/step boundaries; in-flight calls wait
under the evaluation deadline plus bookkeeping. This short Cu workflow is an
interface reference. [B09 evidence](../docs/B09_VERIFICATION.md) records local
full-run startup savings and the limits of both adapter directions.

## Development and support

From the repository root, after building the CLI:

```bash
MODENA_BINARY="$PWD/target/release/modena" python3 -m unittest discover -s python/modena/tests -v
```

Set `MODENA_ASE_PYTHON=/absolute/path/to/venv/bin/python` to execute the optional
actual ASE transport tests; absent this variable they explicitly report skips.
Use a Python 3.11+ interpreter for both default and optional suites. Complete
core/CLI/MCP/direct-worker commands and parity tolerances are in the repository's
[B01 guide](https://github.com/ricardoleal20/modena/blob/main/docs/B01_EXTERNAL_PROVIDER.md#optional-reference-tests)
and [verification record](https://github.com/ricardoleal20/modena/blob/main/docs/B01_VERIFICATION.md).
`MODENA_BINARY` selects the executable for this test suite. For application
code, use `binary=` or put the CLI on `PATH`. Use the
[issue tracker](https://github.com/ricardoleal20/modena/issues) for reports and
questions, and read the
[contribution guide](https://github.com/ricardoleal20/modena/blob/main/CONTRIBUTING.md)
for the full development checks.

This package is licensed under the MIT license. The included `LICENSE` file
contains the same notice as the repository's root license.

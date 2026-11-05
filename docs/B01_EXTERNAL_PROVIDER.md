# B01: pinned ASE EMT external provider

`ase_emt` runs the installed official `ase.calculators.emt.EMT` through the
ordinary TOML runner. The optional Python environment is separate from the
default Rust build. This first reference checks transport, unit conversion,
derivatives and continuation against direct EMT calls; it does not establish
copper material accuracy or agreement with another dynamics engine. Executed
versions, tolerances and local results are in [B01 verification](B01_VERIFICATION.md).

## Supported configuration

Use `units="real"`, `composition="total"`, full periodicity and an
unconstrained fixed cell with NVE or NVT. Supported atomic numbers are Al=13,
Ni=28, Cu=29, Pd=46, Ag=47, Pt=78 and Au=79. `asap_cutoff=false` is the only
supported cutoff convention. A complete explicit `atomic_numbers` vector is
required, in structure order; element display labels and force-field `species`
do not supply chemical identity.

Set mechanical masses yourself, in amu. FCC/XYZ require `mass`, and rocksalt
requires both masses. The existing random-alloy `SpeciesSpec.mass` defaults to
1.0 when omitted; after parsing the runner cannot distinguish that default
from an explicit 1.0. Supply every alloy species mass explicitly too. No mass
is inferred from Z.

Non-real units, partial/nonperiodic PBC, unsupported Z, charges, spin inputs,
topology/electrostatic additions, constraints/water, pressure coupling/NPT and
composition other than `total` are unsupported and fail validation/protocol
checks. Spin and charge arrays are not accepted worker fields. Missing energy,
force or stress properties fail; they are never replaced with zeros. The total
provider receives no automatic native pair, bonded or electrostatic terms.

## Reproduce a short Cu workflow

From the repository root with Rust 1.90+ and Python 3.11+, install an optional
local environment and build this checkout's CLI:

```bash
rtk proxy python3 -m venv .venv-ase
rtk proxy .venv-ase/bin/python -m pip install 'ase>=3.26,<4'
rtk proxy cargo build --locked --release --bin modena
rtk proxy mkdir -p target/b01-example
rtk proxy target/release/modena provider-info \
  --python "$PWD/.venv-ase/bin/python" > target/b01-example/provider.json
```

The supported installation range is not a reproducibility pin. Discovery
reports the actual installed ASE and NumPy versions, bridge protocol, units,
capabilities and EMT source SHA-256. EMT embeds its parameters in source, so
this hash identifies that model data; it is not an ML weights hash. Preserve
the venv executable path: dereferencing its symlink can select base Python.

Generate an input using the discovered exact pins:

```bash
rtk proxy .venv-ase/bin/python - <<'PY'
import json
from pathlib import Path
import sys

out = Path("target/b01-example")
identity = json.loads((out / "provider.json").read_text())
q = json.dumps  # JSON string escaping is valid for these TOML basic strings.
(out / "cu.toml").write_text(f'''units = "real"
[structure]
type = "fcc"
lattice_constant = 3.61
replicas = [1, 1, 1]
element = "Cu"
mass = 63.546
[model]
type = "ase_emt"
python = {q(sys.executable)}
ase_version = {q(identity["ase_version"])}
model_sha256 = {q(identity["model_sha256"])}
atomic_numbers = [29, 29, 29, 29]
composition = "total"
pbc = [true, true, true]
asap_cutoff = false
timeout_ms = 30000
transport = "one_shot" # default; use "persistent" for one process/calculator
[integrator]
type = "velocity_verlet"
dt = 0.0001
[run]
steps = 2
initial_temperature = 10.0
seed = 42
log_every = 1
log_file = "cu.csv"
trajectory_file = "cu.xyz"
trajectory_every = 1
restart_file = "cu.checkpoint"
restart_every = 1
''')
PY
rtk proxy target/release/modena check -i target/b01-example/cu.toml
rtk proxy target/release/modena run --diagnostics-json \
  -i target/b01-example/cu.toml > target/b01-example/run.stdout
```

Outputs resolve relative to `cu.toml`. `dt` is ps, initial temperature K,
length Å, energy kcal/mol, force kcal/(mol Å), and rendered pressure atm.
Pressure uses standard atmospheres with the exact-SI conversion
`4184/(6.02214076e23*1e-30*101325) = 68568.42296625092` atm per
kcal/(mol Å³). The earlier 69477.8 display scale is corrected; see the
[independent pressure evidence](B01_VERIFICATION.md#corrected-pressure-units-and-final-fix-evidence).
Checkpoint mechanics and energy diagnostics remain reduced. This two-step
workflow exercises the interface; it is not an ensemble or stability study.

For a short NVT transport check, replace the integrator table with:

```toml
[integrator]
type = "langevin_baoab"
dt = 0.0001
gamma = 1.0
target_temperature = 10.0
seed = 42
```

Use a separate output directory to preserve the NVE artifacts. The executed
four-Cu NVT example reports 12 kinetic DOF and `not_applicable` conservation,
with three full-rate samples; this does not validate canonical sampling.

## Evaluation and failure behavior

Each snapshot starts a new Python process and constructs fresh Atoms and
`EMT(asap_cutoff=False)`. Protocol `modena-ase-emt-v1` transfers positions in Å,
cell vectors as rows, explicit Z and PBC; responses contain eV, eV/Å and
tensile-positive eV/Å³ stress. The real-unit runner uses reference length 1 Å
and energy 1 kcal/mol; the library also supports explicit `UnitFactors`.
The exact conversion is `eV_to_kcal_mol = 1.602176634e-19 * 6.02214076e23 / 4184`.
Scalar configurational virial is `W = -V trace(stress)`; pressure metadata is
`unconstrained_provider_stress`, with the existing kinetic contribution.
Cold-Cu CLI pressure is independently checked against direct official ASE stress
using exact eV, Angstrom and standard-atmosphere SI definitions.

The provider manages its own periodic neighbors. MODENA does not construct a
classical graph or invent an external cutoff. Identity, property presence,
shape, symmetry and finite converted values are checked before accumulator
updates; failed results leave existing forces/virial intact. Construction
checks the requested ASE version and model hash, and evaluations require the
complete discovered identity, including NumPy, to remain equal.

Transport executes the selected interpreter without a shell. Defaults are
`python="python3"`, `timeout_ms=30000` and `transport="one_shot"`. One-shot
calls use temporary regular files to avoid pipe deadlocks. Output bounds are
polled at 16 MiB stdout and 64 KiB stderr, so transient overshoot is possible.
Diagnostics are limited to 4096 UTF-8 bytes plus short error prefixes.

B09 adds explicit `transport="persistent"` in TOML, Python `AseEmt`, and MCP
configuration. A provider owns one process and one EMT calculator, including
the identity handshake; repeated evaluations reuse both. Changed positions,
cell and numbers invalidate ASE's cache exactly, including changes smaller
than its normal cache tolerance. The Rust provider still requires the pinned
per-atom numbers in every snapshot. Persistent messages have a four-byte
big-endian length followed by strict JSON, with a separate
`modena-ase-emt-session-v1` envelope, increasing request IDs and a fresh
256-bit OS-random nonce echoed by the response. A buffered reply cannot match
the nonce of a future snapshot. The unchanged
`modena-ase-emt-v1` physical identity remains inside each response. Frames are
limited to 16 MiB, lifetime stderr to 64 KiB (only 4096 bytes retained), and IO
queues are bounded. Dedicated IO threads allow the request deadline to cover
both pipe writes and reads even when the worker stops consuming input.
Borrowed snapshots serialize scaled coordinate triples directly into one
capped frame buffer; no per-snapshot JSON trees or complete transport copies
precede the size check. IO threads also validate response correlation and
enforce fatal stdout/stderr failures while the model is idle. Terminalization
disconnects both queues and releases queued frames immediately.

Drop kills and reaps the direct worker child. Timeout, EOF, malformed or stale
replies, changed identity and invalid/missing properties terminate and reap
the persistent process; that provider instance rejects later calls without
respawning. Failed accumulation also terminates it, preserving all existing
force/virial accumulators. Invalid snapshots are rejected before transmission.
There is no retry, stale result fallback or general calculator selection.
These bounds are not an operating-system sandbox.

MCP cancellation is observed before setup, between provider calls,
and before each integration step. One in-flight call may finish under its
configured deadline plus ordinary bookkeeping; this is not a hard global
wall-time guarantee. One-shot startup latency is expected; no throughput claim
follows from transport support alone.

## Checkpoints and public API migration

`System.atomic_numbers: Vec<u8>` is separate from species; new native atoms use
zero (unspecified). `ContinuationState.atomic_numbers` is serialized with a
default empty vector for older v2 files. An empty old vector is compatible only
with unspecified native Z; a nonempty vector must have exactly N entries.
External continuation requires exact Z and validates the runtime provider
before accepting cached mechanics. External signatures include runtime
identity, PBC, task, composition and cutoff convention. Changed Z, version,
hash, NumPy identity or convention is rejected. Interpreter location and
timeout and process transport are settings that do not change physical identity;
compatible checkpoints can continue across one-shot and persistent runs. Native
`modena-physical-v2` identity is preserved except for the corrected real-pressure
factor. A fixed-cell real checkpoint may migrate only when its entire signature
matches after replacing precisely the old pressure factor at index 5 of the
seven-entry `output_unit_factors` array. The old value must exactly equal
`epsilon/(sigma*sigma*sigma)*69477.8`; all other factors, provider identity/Z,
ensemble and mechanical inputs remain exact. A new checkpoint writes the
corrected factor. Old real pressure-coupled checkpoints cannot continue:
`restart_mode="state_only"` explicitly begins a new protocol (including updated
pressure targets/compressibility) and still validates the selected provider.
Reduced-unit continuation is unchanged.

Rust consumers must update struct literals for the new System/continuation
fields and `RunMetadata.provider: Option<ProviderIdentity>`. `RunMetadata` and
`RunSummary` are `Clone`, no longer `Copy`; clone or borrow when needed.
`AseEmtSpec` literals also require `transport`, using
`modena_core::external::AseEmtTransport::OneShot` (or `Default::default()`)
for the compatibility default and `Persistent` for a session.
Exhaustive matches must handle `ModelSpec::AseEmt` and
`PressureStatus::UnconstrainedProviderStress`. The existing by-value
`RunEvent::Finished { summary }` contract is retained. Exhaustive error matches
must also handle `ModenaError::Cancelled`. `Simulation::from_config_with_stop`
accepts a borrowed `AtomicBool` for cancellable setup; existing constructors
remain available. Clean cancellation and successful repeat runs remain usable.
Any run error makes that Simulation terminal: further runs and checkpoint
writes return an error without touching an existing checkpoint. The buffers may
be partially advanced; construct a fresh instance from the last good checkpoint.
`current_step()` records completed steps; Failed events name the attempted step.

Path-like relative `model.python` values (for example `./venv/bin/python`) resolve
against the input directory or MCP `work_dir` in the shared runner. Absolute
paths and bare PATH names (including `.emt-python`) remain unchanged, and
executable/venv symlinks are preserved. A leading dot alone does not make a
basename a relative path. Because the configured interpreter is a String, a
resolved path that cannot be represented as UTF-8 is rejected with a Config
error before launch; it is never converted lossily. Standalone `provider-info` has no input directory and uses its process
working directory.

`modena provider-info [--python EXECUTABLE] [--timeout-ms MILLISECONDS]` emits
one JSON descriptor on success and no JSON on failure. `run --diagnostics-json`
retains `energy diagnostics = {JSON}` and additionally emits a separate
`run metadata = {JSON}` line. Native metadata has `provider:null`. Python
exports `AseEmt` and `provider_info`; appended `RunResult.metadata` is optional
and remains `None` with older CLIs. See the [Python guide](../python/README.md)
for typed configuration/path handling and [MCP guide](../crates/modena-mcp/README.md)
for schema and status transport. MCP external cost estimates explicitly fail;
native pair timing/memory estimates cannot be applied to Python EMT. Async
status carries identity through the lifecycle and completion; the current
integration test does not guarantee observing a nonterminal state.

## Optional reference tests

```bash
rtk proxy env MODENA_ASE_PYTHON="$PWD/.venv-ase/bin/python" \
  cargo test --locked -p modena-core --test external_provider -- --nocapture
rtk proxy env MODENA_ASE_PYTHON="$PWD/.venv-ase/bin/python" \
  cargo test --locked -p modena-cli -p modena-mcp -- --nocapture
rtk proxy env MODENA_ASE_PYTHON="$PWD/.venv-ase/bin/python" \
  MODENA_BINARY="$PWD/target/release/modena" .venv-ase/bin/python \
  -m unittest discover -s python/modena/tests -v
rtk proxy .venv-ase/bin/python crates/modena-core/tests/test_ase_emt_worker.py
```

Rust optional checks print `SKIP` and return successfully when
`MODENA_ASE_PYTHON` is absent; those harness passes do not execute ASE. Python's
wrapper suite reports the optional reference as skipped. The direct worker
suite imports ASE and must be invoked with the optional environment. ASE
3.29.0 with NumPy 2.5.4 emits an upstream shape-assignment DeprecationWarning
from `ase/atoms.py` during the independent direct reference; it is retained,
not broadly filtered. Package inspection must include
`src/external/ase_emt_worker.py` in `cargo package --list -p modena-core`.

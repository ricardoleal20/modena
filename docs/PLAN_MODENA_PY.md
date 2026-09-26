# Plan — `modena-py` (Python client) and the Rust config language

Status: **planned** (2026-09-26). Drives the current `modena` CLI + keeps the
engine as the single source of truth for physics and config validation.

## 0. Why and the two deliverables

1. **`modena-py`** — a Python package to *configure* and *drive* MODENA runs
   from Python, so scientific pipelines, notebooks, and agent workflows can
   build systems, pick models/integrators/ensembles, run, and read back
   thermo + conservation verdicts without shelling out by hand.
2. **The Rust "modena-language"** — a typed, fluent, programmatic configuration
   API in Rust. Today Rust configures MODENA only via TOML
   (`InputConfig::from_str` → `Simulation::from_config`). Plan: a
   `SimulationBuilder` DSL that produces a validated `InputConfig` with the
   same checks (units conversion, cadence ≥ 1, species range, neighbor_skin),
   so Rust applications can configure runs in code exactly as a Python user
   configures them in a dict.

**Single-source-of-truth principle.** All config validation lives in Rust
(`InputConfig::from_str` → `validate()`). Both the Python layer and the Rust
builder must route through it — never re-implement validation in Python.

## 1. Current state (verified 2026-09-26)

- TOML schema: `units`, `[structure]` (fcc | xyz | random_alloy | rock_salt),
  `[model]` (lennard_jones | buckingham | morse | soft_sphere | tabulated, with
  per-type tables + overrides), `[electrostatics]` (Ewald), `[topology]`
  (bonds/angles/dihedrals), `[constraints]` (SHAKE/RATTLE),
  `[integrator]` (velocity_verlet | langevin_baoab | nose_hoover_chain),
  `[thermostat]`, `[barostat]`, `[run]` (`steps`, `initial_temperature`,
  `seed`, `neighbor_skin`, `trajectory_file`, `trajectory_every`, `log_file`,
  `log_every`).
- Rust programmatic seam today: `InputConfig` fields are `pub` (constructible
  in Rust), and `Simulation::from_config(&InputConfig, base_dir)` builds
  + runs. **No fluent builder exists.**
- Python today: only the demo pattern (dict → hand-written TOML → subprocess).
  No package, no schema mirror, no event access.
- `modena-mcp` exposes the run over stdio MCP (launch, monitor events,
  interrupt, sweep) — a standard interface any Python MCP client speaks.

## 2. `modena-py` — transport options (decision)

| Option | Event/interrupt | Work | Verdict |
|---|---|---|---|
| **A. Subprocess the CLI** (`modena run -i cfg.toml`) | none live (parse CSV/XYZ at end) | hours | Ship for v0.1 (zero deps, works everywhere) |
| **B. MCP client** (official `mcp` SDK over stdio of `modena-mcp`) | full: streamed `Thermo`, interrupt, sweep, all verdicts | ~2-4 days | Ship for v0.2 — this is the agent-native claim |
| **C. PyO3 native bindings** (`modena._engine`, maturin) | in-process, fastest | ~1-2 weeks | v0.3 (optional) |

Recommended: **A now, B next** (they share the same `modena.config` schema
mirror); C only if a driver needs in-process structs (rare with the CLI).

### Python package layout

```
python/
├── pyproject.toml            # (name "modena", uv/pip installable, zero deps)
├── modena/
│   ├── __init__.py           # from .config import (Fcc, Ld..., Argon...)
│   ├── _toml.py              # low-level: dataclass->TOML renderer (TOML 1.0)
│   ├── config.py             # typed mirrors of the schema (dataclasses)
│   ├── run.py                # subprocess driver (A): run()/check(); parses
│   │                         #   thermo CSV + final summary + conservation
│   ├── mcp_client.py         # (v0.2) async client (B) over stdio MCP
│   └── postprocess.py        # read XYZ/CSV into arrays (numpy optional)
└── tests/                    # parity tests: python cfg TOML == hand TOML
```

### Python API sketch (v0.1, transport A)

```python
from modena.config import Fcc, LennardJones, VelocityVerlet, Run
from modena.run import check, run

cfg = LennardJones(sigma=1.0, epsilon=1.0, cutoff=2.5)
r = run(
    structure=Fcc(lattice_constant=1.6796, replicas=[2, 2, 2], element="Ar"),
    model=cfg,
    integrator=VelocityVerlet(dt=0.005),
    run=Run(steps=2000, initial_temperature=0.8, seed=42, neighbor_skin=0.3),
)
print(r.final_conservation)   # ConservationVerdict: Healthy | Warning | Drifting
print(r.perf_us_per_step)
assert r.frame("traj.xyz").positions.shape == (32, 3)
```

Design rules:

- **Parity by construction**: `config.py` renders TOML *exactly* like the
  schema (key names/order via a single renderer), so `InputConfig::from_str`
  is the validator — Python never re-implements validation or unit conversion.
- Dataclasses with `=None` optional fields match serde defaults
  (`neighbor_skin=0.3`, `trajectory_every=100`, `log_every=10`, …).
- `modena.pyproject` defines no runtime deps for transport A
  (stdlib `subprocess` + `csv` + `tomllib`); numpy optional.
- `mcp_client.py` (v0.2) drives `modena-mcp` for streaming/agent workflows —
  the exact surface an LLM agent uses, so "Python client" and "agent-native"
  share one interface.

### Milestones

- **M-PY-0** (½ d): `python/modena` skeleton + `_toml.py` renderer +
  `config.py` (Fcc/Xyz/RandomAlloy/RockSalt × 5 models × 3 integrators ×
  thermostat/barostat × topology/constraints) + `run.py` `check`/`run` +
  parity tests that render the decks in `bench/validation/*/modena.toml` and
  diff against them byte-for-byte. Acceptance: `modena-py` reproduces the
  SHAKE/Ewald decks.
- **M-PY-1** (½ d): `modena.pyproject`, README install (`pip install -e ./python`)
  + a `examples/python_driver.py` (the demo from the audit, packaged).
- **M-PY-2** (2-4 d): `mcp_client.py` streaming client (events, interrupt,
  sweeps) + tests against a real `modena-mcp` session.
- **M-PY-3** (1-2 wk, optional): PyO3 `_engine` for in-process runs.

## 3. The Rust "modena-language" — `SimulationBuilder`

Goal: configure runs from Rust code with the same ergonomics as the Python
dict API and the same validation as TOML.

### Design (new: `crates/modena-core/src/config/builder.rs`)

```rust
let sim = Simulation::builder()
    .fcc(FccBuilder::new(1.6796).replicas([4, 4, 4]).element("Ar"))
    .model(ModelBuilder::lj().sigma(1.0).epsilon(1.0).cutoff(2.5))
    .integrator(IntegratorBuilder::velocity_verlet().dt(0.005))
    .run(RunBuilder::new().steps(2000).temperature(0.8).seed(42).skin(0.3))
    .build()?; // runs InputConfig::validate() + from_config
```

Implementation notes:

- Thin: every `*Builder` produces the existing `InputConfig` types (fields are
  already `pub`), then calls `InputConfig::validate()` before construction — so
  **one** validator serves TOML, Python, and Rust.
- Model/type-pair tables + overrides use the same `*TypeParams`/`*Override`
  structs as the TOML path; type safety where the schema has it
  (`sigma > 0`, `n_types` vs species, cadences ≥ 1).
- Optional `LatticeBuilder` sugar (fcc/rock_salt/random_alloy) returning the
  builder so Rust tests/examples stop hand-rolling `System`s.
- `build()` returns `Simulation` (or `Simulation::builder()` on
  `InputConfig`, keeping the seam discoverable).

### Milestones

- **M-RL-0** (1 d): `SimulationBuilder` + `RunBuilder` + `ModelBuilder`
  (all 5 potentials incl. per-type + overrides) + `validate()` re-use +
  one integration test (builder-built run == hand-TOML run).
- **M-RL-1** (½ d): structure/electrostatics/topology/constraints/
  thermostat/barostat builders + docs (book page "Configuration in Rust").
- **M-RL-2** (½ d): parity must — builder output `InputConfig` serializes
  back to TOML byte-identical to the hand-written deck for all
  `bench/validation` cases.

## 4. Sequencing vs the paper

| Item | Before first paper? | Rationale |
|---|---|---|
| M-PY-0 + M-PY-1 (Python driver) | ✅ yes | Reproducibility section: pipelines call `modena` from Python; small, low-risk |
| M-RL-0 (Rust builder) | ✅ yes | Makes the "config square" (TOML = Python = Rust) a one-line paper claim |
| M-PY-2 (MCP client) | 🟡 nice | It *is* the Feature-A demo; but agents already use `modena-mcp` directly |
| M-PY-3 (PyO3) | ❌ after | effort ≫ benefit for the paper |
| Rust SIMD/hot-loop work | see perf roadmap | separate PR(s), measured |

## 5. Risks

- **Drift between schema and Python mirror**: mitigated by byte-parity tests
  (render the shipped decks and diff) run in CI.
- **Two transports (A, B) duplicating logic**: kept to parity by sharing
  `config.py`; `mcp_client.py` only adds transport.
- **PyO3 resists .toml-path**: not started until M-PY-3 decides.

## 6. References

- Schema source of truth: `crates/modena-core/src/io/input.rs`.
- Runner seam: `Simulation::from_config` (`crates/modena-core/src/runner.rs`).
- `modena-mcp`: stdio MCP server (`crates/modena-mcp/`).
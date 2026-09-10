# MODENA — Feature A (Agent-Native Interface) + Feature B (Conservation Contract) Design

**Status:** Design draft for user review (brainstorming-skill architectural path).
**Date:** 2026-09-10
**Builds on:** `docs/NOVELTY_PROPOSAL.md`

---

## 1. Goal & scope

Integrate two features from the novelty proposal into MODENA, plus a shared
refactor that both require:

- **Feature B — Conservation contract.** Make "the reported potential energy is
  the exact integral of the reported force" an *enforced, tested engine
  contract*, verified by property-based tests. Low risk; pure upside.
- **Feature A — Agent-native interface.** Expose MODENA's simulation primitives
  as typed, validated tools an LLM agent can call (build system, run probe,
  launch sweep, read thermo, detect anomalies) — not by shelling out to the
  CLI, but by calling the engine's library API behind a control surface.
- **Shared enabler — observable runner.** Both features need the run loop to
  emit step-level state and accept interruption. Today `Simulation::run()` is a
  blocking run-to-completion loop writing only to files. This is the one core
  change both features share; we design it once.

**Out of scope (deferred):** the Tauri GUI (designed in parallel as a separate
surface over the same runner API), real ML force fields (Feature C), WASM build
(Feature D), Verlet skin (roadmap M7).

### Non-goals (honest framing)

- We do **not** claim novelty in differentiable MD, equivariant ML potentials,
  or GPU MD (see NOVELTY_PROPOSAL §1).
- The agent interface is the *primary* novelty claim; Feature B is the
  *methodology/trust* pillar that supports it. Both must be demonstrable.

---

## 2. The shared enabler: an observable, interruptible runner

### Problem

`Simulation::run()` (`crates/modena-core/src/runner.rs:95`) loops `1..=steps`,
calls the integrator, and writes to file writers. An agent and a GUI both need:

1. **Streaming** — step-level `Energies` + step number delivered to a callback
   or channel, not only to disk.
2. **Interruption** — stop a run cleanly mid-step (agent cancels; user hits
   pause).
3. **Multiple concurrent runs** — an agent sweeps launch N runs; a GUI may show
   one foreground + background sweeps. Each needs its own loop instance.
4. **Progress + status** — running / paused / done / failed, with current step.

### Design

Make the run loop drive an **event channel** and check a **cancel token**:

```rust
/// Emitted by the runner as it steps. Consumed by the MCP layer, the GUI,
/// and (optionally) the file logger.
#[derive(Clone, Debug, Serialize)]
pub enum RunEvent {
    Started { atoms: usize, dt: f64, steps: usize },
    Thermo { step: usize, time: f64, energies: Energies },
    Frame { step: usize, time: f64, system_snapshot: SystemSnapshot },
    Finished { summary: RunSummary },
    Failed { error: String },
}
```

```rust
pub struct RunHandle {
    /// Events streamed from the running simulation.
    pub events: Receiver<RunEvent>,
    /// Signal the run to stop after the current step.
    pub cancel: Sender<()>,
    /// Join the runner thread (blocking) to retrieve the final summary.
    pub join: JoinHandle<Option<RunSummary>>,
}
```

`Simulation` gains a `run_streaming(&mut self, cancel: Receiver<()>, sink: Sender<RunEvent>)`
that yields a `RunEvent::Thermo` every `log_every` steps (and `Frame` every
`trajectory_every`), checks `cancel` each step, and returns the summary.

`run()` (blocking, file-based) stays as a thin wrapper over
`run_streaming` + a file-writing sink, so the CLI and existing tests are
unchanged.

**Concurrency model:** each run is a `std::thread` (or, later, a rayon task).
The agent's sweep = a small pool of these. MODENA's low RAM (~3–6 MB/run) is
what makes a pool of, say, 8–16 concurrent runs affordable — this is the
"low-RAM enables dense agent streams" claim from the proposal, made concrete.

**Why not async (tokio):** the force kernel is CPU-bound rayon work, not I/O.
Threads + channels keep the core dependency-free; the MCP layer may use tokio
externally but the core stays sync + std.

**Risk:** `System` is `Clone` (added for benchmarks) but cloning every frame for
`Frame` events is expensive. Resolution: `SystemSnapshot` is a cheap
positions-only `Vec<f32>` (f32 for transfer; positions only, not forces),
emitted at `trajectory_every`, not every step. The GUI/MCP can subsample
further.

---

## 3. Feature B — the conservation contract

### The contract

> **For every `PairPotential`, the energy `U(r)` returned by `energy_force`
> must be a function whose radial derivative equals `−f_over_r · r`, i.e.
> `dU/dr = −F(r)` where `F(r) = f_over_r · r`.** Equivalently: the potential
> energy the engine reports is the exact integral of the force the engine
> applies. This makes the NVE Hamiltonian `H = KE + PE` the true conserved
> quantity.

This is *already true* for `LennardJones` (shifted-force, analytic). Feature B
makes it **enforced** rather than incidental.

### How it's enforced

**3.1 Trait-level documentation contract.** Update `PairPotential::energy_force`
doc to state the invariant explicitly and note that it is property-tested.

**3.2 Property-based tests** (`crates/modena-core/tests/property_consistency.rs`,
using `proptest`):

For a registered potential (start with LJ; extensible to future potentials),
for random configurations of N atoms in a box with r < cutoff:

1. **Force = −∇U (finite difference).** For each atom i, perturb its position
   by ±h in each axis and check
   `F_i,k ≈ −(U(+h) − U(−h)) / (2h)` to within `h²` truncation + tolerance.
2. **Pair self-consistency.** For a single pair at random r < cutoff,
   `dU/dr ≈ −F(r)` via finite difference on the `energy_force(r2)` output directly.
3. **No NaN/Inf.** For random r² in `(0, cutoff²]` and random species,
   outputs are finite. (Rejects degenerate inputs.)
4. **Cutoff continuity.** Force at `r = rc` is ~0 (shifted-force property);
   energy at `r = rc` is continuous (limit from below equals the value).
5. **NVE drift bound.** A short NVE run of M steps has relative energy drift
   `< C · dt⁴` (symplectic velocity-Verlet is 4th-order in the shadow
   Hamiltonian). This generalizes the existing hard-coded conservation test
   into a parametrized property.

**3.3 A `ConservationReport` type** the runner can emit:

```rust
pub struct ConservationReport {
    pub e0: f64,
    pub e_final: f64,
    pub relative_drift: f64,    // (E_final - E0) / |E0|
    pub verdict: ConservationVerdict,  // Healthy / Warning / Drifting
}
```

`verdict` thresholds (LJ, dt=0.005): `<1e-4` Healthy, `<1e-3` Warning, else
Drifting. This is the structured signal an agent (and the GUI's conservation
gauge) interprets — *drift is not an error, it's data*. This is the "trust
pillar" the agent reasoning trace references ("drift 4.7e-6 → healthy").

### Testing strategy

- `proptest` cases run fast (N ≤ 64, ≤ 200 steps) so they fit in CI.
- Existing `tests/conservation.rs` stays; the property test generalizes it.
- CI runs `cargo test` (includes proptest) — add `proptest` to dev-deps.

### Effort

Low. ~1 day. No core algorithm change — LJ already satisfies the contract.

---

## 4. Feature A — the agent-native interface (`modena-mcp`)

### The gap (recap)

LLM science agents (Coscientist, ChemCrow) drive *experiments/chemistry*; for
MD they generate LAMMPS/GROMACS input scripts blindly and parse text output. No
MD engine exposes **typed simulation primitives** as an agent-callable surface.
MODENA's low RAM makes dense parallel agent streams tractable. This is the
primary novelty.

### Design

A new crate `crates/modena-mcp` exposing MODENA primitives as **MCP tools**
(consumable by any MCP-capable agent — Claude, etc.). MCP is the natural fit
because it's the emerging standard for tool-calling and Cherry Studio already
speaks it.

**4.1 The tool set.** Each tool = typed JSON in/out, validated against the
existing `InputConfig` types.

| Tool | Input | Output | Notes |
|---|---|---|---|
| `list_potentials` | — | `[{name:"lennard_jones", params:[σ,ε,cutoff]}]` | Discovery |
| `build_system` | structure spec (fcc/sc/xyz, lattice/density, replicas) | `{atoms, density, box}` | Validates; returns derived counts |
| `estimate_run_cost` | system + steps | `{est_wall_s, est_rss_mb}` | From benchmark model; helps agent budget |
| `run_probe` | full config (structure+model+integrator+run) + steps (≤500) | `{run_id, final_thermo, conservation}` | Short synchronous probe; returns immediately |
| `run_async` | full config + steps | `{run_id}` | Launches a streaming run; agent polls |
| `get_run_status` | `run_id` | `{state, step, latest_thermo, conservation}` | Poll a run |
| `get_run_events` | `run_id`, since_step | `[{step, thermo, ...}]` | Streamed events |
| `cancel_run` | `run_id` | `{state:"cancelled"}` | Clean stop |
| `launch_sweep` | base config + swept param + range[] | `{sweep_id, run_ids[]}` | Parallel runs (T ladder, density scan) |
| `get_sweep_results` | `sweep_id` | `[{run_id, params, final_thermo, conservation}]` | Aggregated |
| `detect_anomaly` | `run_id` or `sweep_id` | `{type, where, evidence}` | Engine-side heuristics: drift spike, T jump (phase transition), NaN |
| `checkpoint` / `resume` | `run_id` | path/id | Pause + later resume |

**4.2 Architecture.**

```
  MCP client (agent) ──JSON-RPC──▶ modena-mcp server (tokio)
                                        │
                          ┌─────────────┴─────────────┐
                          ▼                           ▼
                   RunManager (registry)      Tool handlers
                   - spawns RunHandle          - validate → call core
                   - tracks run_id → handle    - format JSON responses
                   - cancel / status
                          │
                          ▼
                   modena-core (Simulation::run_streaming)
```

`RunManager` owns the `HashMap<RunId, RunHandle>` and the run thread pool. It
maps MCP tool calls to `Simulation` operations. The MCP server is a thin
transport; all logic is `RunManager` + `modena-core` — so the GUI (Tauri) can
reuse `RunManager` directly as a Rust library, sharing one code path with the
agent. **One core, two surfaces** (the diagram from brainstorming).

**4.3 Anomaly detection — the agent's eyes.** This is what makes agent
"discovery" real rather than script-generation:

- **Drift spike:** `ConservationReport.verdict` degrades mid-run → alert.
- **Phase transition (sweep):** across a T-ladder, detect a discontinuity /
  inflection in `⟨PE⟩(T)` or `C_v(T)` → flag "possible transition at T≈x".
- **Divergence:** NaN/Inf or KE blow-up → mark run failed with cause.

These run **engine-side** (cheap, deterministic heuristics), so the agent
doesn't have to re-derive them from raw numbers. The agent's reasoning then
operates on *structured findings* ("transition detected at T≈1.08"), not
floating-point arrays.

**4.4 Validation & safety.** Every tool input is validated against
`InputConfig` (rejects negative σ, cutoff ≤ 0, N too large, etc.) before any
run starts. `estimate_run_cost` lets the agent (or the MCP layer) refuse runs
that would exceed a configured atom-step budget. This is the "safe, typed
primitives" claim — contrasted with blind LAMMPS script generation.

### The demonstration (what the paper shows)

An LLM agent, given a goal like *"find the melting temperature of the LJ
liquid at ρ=0.8442"*, autonomously:
1. builds the system via `build_system`,
2. launches a T-sweep via `launch_sweep` (e.g. T = 0.8, 0.9, 1.0, 1.1, 1.2),
3. polls `get_sweep_results`, calls `detect_anomaly`,
4. receives "transition detected near T≈1.08",
5. refines: `launch_sweep` T = 1.05…1.12,
6. reports the melting T with confidence + conservation evidence.

All steps are *typed tool calls* over a shared core — no shell, no script
generation, no text parsing. This is the hero demo for the paper.

### Effort

Medium. `modena-mcp` crate + `RunManager` + tool handlers + an MCP SDK
(`rmcp` or hand-rolled JSON-RPC over stdio). Depends on the runner refactor
(§2) and is cleaner if Feature B's `ConservationReport` exists first.

---

## 5. Integration order

```
  [Feature B: contract + property tests]        ← de-risks first, ~1 day
            │  (adds ConservationReport)
            ▼
  [Runner refactor: run_streaming + RunHandle]  ← shared enabler, ~2-3 days
            │
      ┌─────┴──────┐
      ▼            ▼
 [modena-mcp]   [Tauri GUI]                     ← parallel surfaces
   (Feature A)   (separate spec)
```

B first (lowest risk, immediately useful), then the runner refactor, then A
and the GUI fan out in parallel over the same `RunManager`/`run_streaming`.

---

## 6. Testing & verification

- **Feature B:** proptest suite in CI; existing conservation test
  generalizes; `ConservationReport` unit-tested with synthetic drifts.
- **Runner refactor:** `run()` behavior unchanged (CLI/tests green);
  `run_streaming` tested with a recording sink + an early-cancel; verify cancel
  stops within one step.
- **Feature A:** each tool handler unit-tested with mock `RunManager`; an
  integration test runs a real short sweep through the MCP surface and checks
  aggregated results; anomaly heuristics tested with constructed inputs
  (a synthetic drift spike, a synthetic T-jump).
- **Paper reproducibility:** the hero demo (agent finds melting T) scripted as
  a deterministic integration test with a fixed agent transcript, so the
  result is reproducible without a live LLM.

---

## 7. Open questions for the user

1. **MCP transport:** stdio (local agent, simplest) vs. HTTP/SSE (remote agent,
   Tauri can share)? Recommend **stdio first** (Cherry Studio / local Claude),
   HTTP later.
2. **Agent SDK:** use the Rust `rmcp` crate, or hand-roll minimal JSON-RPC?
   Recommend `rmcp` if mature enough — **[VERIFY]** its current state.
3. **Run pool size default:** how many concurrent runs? Recommend a configurable
   cap (default 8), gated by `estimate_run_cost`.
4. **Anomaly detection scope for v1:** just drift + NaN, or include the
   phase-transition sweep detector? Recommend **drift + NaN + basic transition
   inflection** for v1 (the latter is the demo's payoff).

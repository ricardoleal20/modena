# MODENA — Architecture & Methods

This document describes MODENA's algorithms at the precision needed for the
Methods section of a paper: the force/energy convention, neighbor search,
integration, parallel decomposition, and the two novelty features (the
agent-native interface and the conservation contract). It is a *methods*
document, not a benchmark report — for measured numbers see
`bench/BENCHMARKS.md` (single-core vs LAMMPS), `bench/scaling/SCALING.md`
(parallel scaling), and `bench/validation/VALIDATION.md` (non-LJ physics
validation).

## 1. Design overview

MODENA is a classical molecular-dynamics engine written in Rust, structured as
a Cargo workspace:

| Crate | Role |
|---|---|
| `modena-core` | Potentials, neighbor search, integrators, thermostats/barostats, electrostatics, bonded + constraints, parallel decomposition, the `ForceModel`/`Integrator` trait seams, and the streaming runner. |
| `modena-cli` | Thin binary that parses a TOML `InputConfig` and drives `Simulation::run`. |
| `modena-mcp` | Agent-native MCP server exposing the streaming runner to LLM agents (Feature A). |
| `modena-bench` | Criterion micro-benchmarks (force kernel, integrator step, neighbor build). |

The two central abstraction seams are the `ForceModel` trait (computes all
forces + PE + virial for a configuration) and the `Integrator` trait (advances
the system one step, calling a `ForceModel`). A new potential — analytical,
many-body, or a learned ML surrogate — implements `ForceModel` and drops into
the engine without touching the integrator, thermostat, or I/O. This is the
extension point a NequIP/MACE-style model will use.

## 2. Force / energy convention (the conservation contract)

**Every pair potential in MODENA returns `(U, f_over_r)` such that the reported
potential energy is the exact integral of the reported force:**

```
dU/dr = −F(r),   where F(r) = f_over_r · r
```

This is enforced at the `PairPotential` trait level and verified for every
registered potential by property-based tests (`tests/property_consistency.rs`:
force = −∇U by finite differences, no NaN/Inf, continuity at the cutoff). The
consequence is that the NVE Hamiltonian `H = KE + PE` is the *true* conserved
quantity under the symplectic integrator — energy drift is a genuine numerical
signal, not a truncation artifact. This is the conservation contract (Feature B;
see §7).

### 2.1 Shifted-force form

All shipped pair potentials (Lennard-Jones, Buckingham, Morse, soft-sphere, and
tabulated) use the **shifted-force** form (Toxvaerd, 1993), which makes both the
energy and the force continuous (→0) at the cutoff:

```
U_sf(r) = U(r) − U(rc) + (r − rc)·F(rc)
F_sf(r) = F(r) − F(rc)
```

where `F(rc)` is the *raw* radial force magnitude at the cutoff. With
`f_over_r = F_sf / r`, the identity `dU_sf/dr = −F_sf` holds exactly. The
force-shifted form is what gives velocity-Verlet its good NVE conservation under
a truncated potential: the force is continuous at `rc`, so there is no impulse
as a pair crosses the cutoff.

**LAMMPS correspondence.** MODENA's shifted-force LJ corresponds to LAMMPS
`pair_style lj/smooth/linear` (also force-shifted, `F(r) − F(rc)`). On the
matched FCC reference configuration the two engines agree on potential energy to
~8×10⁻⁶ (`bench/BENCHMARKS.md` §2), confirming the potential forms are
numerically equivalent for the LJ benchmark. Note the energy zero convention
differs: MODENA's `U_sf` is the integral of the shifted force, so its NVE
Hamiltonian is exactly conserved; LAMMPS's `lj/smooth/linear` uses a shifted
energy that is not the exact integral of its shifted force. This is the root
cause of MODENA's tighter energy conservation in the LJ benchmark.

### 2.2 Multi-species parameters

Multi-species systems use per-type-pair `(σ, ε)` tables combined with the
Lorentz–Berthelot rule by default (`σ_ij = (σ_i+σ_j)/2`, `ε_ij = √(ε_i·ε_j)`),
with optional explicit per-pair overrides. The cutoff is uniform across all
pairs (the standard assumption for cell-list neighbor search). The per-pair
`U(rc)` and `F(rc)` are precomputed at construction, so the inner loop does one
table lookup per pair.

## 3. Neighbor search

A cell (Verlet) list with a **skin** (M7). The box is binned into cells of side
≥ `cutoff`; pairs are enumerated within each cell and its 13 forward neighbor
cells (a half list — see §4). The skin extends the bin radius by a configurable
margin so the list need not be rebuilt every step: rebuild occurs only when any
atom has moved more than half the skin since the last build, checked cheaply via
a per-atom displacement bound. Minimum-image distances are computed against the
*global* box (the `SimBox::min_image` operation), supporting orthorhombic
periodicity. Triclinic boxes are supported at the `SimBox` seam level; the
parallel decomposition (§5) and barostats are orthorhombic-only.

## 4. Half-list force evaluation (serial)

`PairForceModel` evaluates the `PairPotential` over a **half neighbor list**:
each pair `(i, j)` is stored once, and Newton's third law writes both
`f_i += f_over_r·r_ij` and `f_j −= f_over_r·r_ij`. This halves the pair work
relative to a full list. The loop is parallelized with rayon: each thread folds
into its own `3N` force accumulator, and the per-thread accumulators are summed
by a reduce. This keeps the pair loop race-free without atomics. The scalar
pressure virial is accumulated as `W = Σ_ij r_ij·F_ij` (each pair once);
pressure (reduced) = `(2·KE + W)/(3V)`.

**M3 accumulate contract.** When bonded interactions or Coulomb terms are
present, a `CompositeForceModel` runs several `ForceModel`s in sequence. Each
model *accumulates* into `system.forces` and `virial` (via `add_force`); the
caller (the integrator or the composite) zeroes forces + virial once before the
pass. This lets a pair model compose with bonded / electrostatic models without
one clobbering the other.

## 5. Parallel decomposition (M9)

Spatial domain decomposition for distributed-memory (MPI) and GPU execution.

### 5.1 Ownership and the ghost halo

The orthorhombic periodic box is partitioned into a `nr[0]×nr[1]×nr[2]` grid of
subdomains — one per rank (`RankGrid`). Each rank *owns* the atoms whose centers
fall in its subdomain; atoms in the 26 periodic neighbor subdomains are imported
as read-only *ghosts* so the local pair kernel sees a complete neighborhood. The
grid is constrained (`RankGrid::balanced`) so every subdomain edge is `≥ cutoff`,
which makes a **single shell** of neighbor subdomains a sufficient ghost halo: a
pair with true minimum-image distance `< cutoff` has both partners present as
owned-or-ghost on the owning rank.

### 5.2 Work partition — bit-identical to serial

Pair work is partitioned by ownership of the **lower-global-id atom**: rank `r`
evaluates pair `(i, j)` (with `gid(i) < gid(j)`) iff atom `i` is owned by `r`.
The partner `j` is either owned by `r` (a *local* pair — Newton's third law
writes both `f_i` and `f_j`) or a *ghost* on `r` (its force is written to the
local ghost-force buffer and reduced back to its owner). **Each pair has exactly
one evaluator and each atom exactly one force accumulator**, so the decomposed
force/energy/virial summed over ranks is bit-identical to the serial half-list
evaluation — the correctness anchor verified in `tests/decomp.rs` (single- and
multi-rank == serial, halo sufficiency, post-migration match, decomposed-NVE
conservation) and across two real MPI ranks in `tests/mpi_decomp.rs`.

### 5.3 MPI halo exchange (`mpi` feature)

The distributed driver (`mpi_comm.rs`) performs the LAMMPS-style exchange each
force evaluation:

1. **Forward exchange** — each rank sends its owned atoms (position, species,
   global id, packed into one `f64` buffer) to its 26 neighbor ranks; received
   atoms become local ghosts.
2. **Local force eval** — a neighbor list over owned + ghosts (minimum-image
   against the global box); each owned atom `i` evaluates pairs `(i, j)` with
   `gid(i) < gid(j)`. Force on the owned `i` is applied directly; force on a
   ghost `j` is accumulated into a ghost-force buffer.
3. **Reverse exchange** — ghost-force buffers are sent back to the ghosts'
   owners, which sum them into their owned forces (Newton's third law across
   rank boundaries).
4. **All-reduce** — potential energy and virial are summed across ranks.

The point-to-point exchange is **deadlock-free by canonical ordering**: with
each neighbor, the higher-ranked process sends first and the lower-ranked
receives first, so no pair of blocking send/receive calls ever wait on each
other. The `LocalCommunicator` simulates `P` ranks in one process (forces
written directly to a shared `System`, so the reverse reduce is a no-op) — used
by the in-process example and the decomp unit tests, which need no MPI runtime.

**Atom migration.** After atoms cross a subdomain face, ownership is recomputed
from current positions; the migration count and re-partition are verified by
`tests/decomp.rs::ownership_migration_after_displacement`.

### 5.4 GPU pair kernel (`gpu` feature)

`GpuPairForceModel` offloads the LJ pair evaluation to a WGSL compute shader via
`wgpu`. Each thread computes one atom's force (the all-pairs form, f32): the
shifted-force energy/force matching the CPU path, half-pair energy and half-pair
virial (each pair visited twice across the two threads), and orthorhombic
minimum-image `d − box·round(d/box)`. The constructor returns `None` when no
adapter is present, so headless CI self-skips gracefully. The kernel computes in
f32 against the CPU's f64, so it is validated to single-precision agreement
(`tests/gpu_pair.rs`) and checked for NVE conservation under integration (the f32
drift tolerance is looser than the f64 CPU path but still far below any
non-conservative failure).

## 6. Integration & ensembles

**Velocity-Verlet (NVE).** The symplectic half-kick / drift / force / half-kick
form. `system.forces` holds `F(t)` on entry and `F(t+dt)` on exit, so the runner
computes forces once before the loop and the first half-kick uses the cached
`F(t)` — one force evaluation per step. With `constraints` set, SHAKE is applied
to the post-drift positions and RATTLE to the post-kick velocities (the
symplectic RATTLE scheme on the constraint manifold); tolerance 10⁻¹⁰, capped at
100 iterations.

**Langevin BAOAB (NVT).** The Leimkuhler–Matthews BAOAB splitting — a
stochastic, config-space-isotropic integrator with favorable sampling properties
at large `dt`.

**Thermostats / barostats (M5).** Velocity rescaling, Berendsen, and a Nosé–Hoover
chain for NVT; Berendsen and Parrinello–Rahman barostats for NPT. Replica
exchange (parallel tempering) is implemented as an outer driver over
independent replicas.

## 7. Novelty features

### 7.1 Feature A — agent-native interface (`modena-mcp`)

The streaming runner (`Simulation::run_streaming`) emits a step-level
`RunEvent` stream (`Started`, `Thermo`, `Frame`, `Finished`, `Failed`) and
accepts an early-stop flag. The `modena-mcp` crate exposes this over MCP, so an
LLM agent can launch a run, observe thermodynamic samples and trajectory frames
as they are produced, and interrupt. File output (CSV log, XYZ trajectory) is
written exactly as in the blocking CLI path, in addition to the callback — so a
run observed by an agent produces the same on-disk artifacts as a CLI run. The
blocking `run()` is `run_streaming` with no stop flag and a no-op sink.

### 7.2 Feature B — conservation contract

The force/energy consistency invariant (§2) is encoded as a trait-level
contract and verified by property-based tests for every potential, not just the
production LJ. The `ConservationReport` classifies a run's energy drift
(`Conserved` / `Drifting` / `Diverging`) with explicit thresholds, so a drift
signal is surfaced as a first-class verdict rather than buried in a log. The
contract is what makes the parallel decomposition's "bit-identical to serial"
claim meaningful: if the decomposed forces were not conservative, NVE
conservation would break even though the one-shot force match passed.

## 8. Reproducibility

- **Build:** `cargo build --release` (binary `target/release/modena`).
  `RUSTFLAGS="-D warnings"`; `rust-toolchain.toml` pins the toolchain.
- **LJ vs LAMMPS:** `bench/run_comparison.sh` runs the matched pairs and emits
  `bench/results/comparison.csv`; methodology + results in
  `bench/BENCHMARKS.md`. Requires `lmp_serial` on PATH.
- **Parallel scaling:** examples under `crates/modena-core/examples/` driven by
  `mpiexec`; results in `bench/scaling/SCALING.md`.
- **Non-LJ validation:** `bench/validation/VALIDATION.md` + matched decks.
- **Tests:** `cargo test --workspace` (default); `--features mpi` and
  `--features gpu` exercise the parallel backends. `cargo fmt --check` and
  `cargo clippy --all-targets -- -D warnings` (plus the `mpi`/`gpu` feature
  variants) are CI gates.
- **CI:** `.github/workflows/ci.yml` — fmt, clippy (default+gpu), MPI job
  (`mpiexec -n 2`, `--oversubscribe` on single-slot runners), GPU compile+skip,
  build+test on ubuntu+macos, bench-compile.

All benchmark numbers in the repo are measured on the recorded hardware, never
extrapolated; configurations too slow to run are skipped and documented as such.

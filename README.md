# MODENA

**Molecular Dynamics Engine for Numerical Atomistics**

A fast, extensible molecular-dynamics framework written in Rust. MODENA is
designed as a LAMMPS-class engine built around a small set of trait seams so
that new physics (many-body potentials, bonded terms, long-range electrostatics)
and new acceleration strategies (ML surrogate potentials, GPU kernels, domain
decomposition) can be added without touching the integrator, thermostat, or I/O.

> **Status — feature-complete (v0.1).** Five shifted-force pair potentials
> (Lennard-Jones, Buckingham, Morse, soft-sphere, tabulated) with multi-species
> Lorentz–Berthelot mixing; bonded interactions (bonds/angles/dihedrals) with
> SHAKE/RATTLE constraints; Coulomb + Ewald long-range electrostatics; NVE/NVT
> (Berendsen, velocity-rescale, Nosé–Hoover chain) and NPT (Berendsen,
> Parrinello–Rahman) ensembles; replica exchange; MPI domain decomposition and a
> WGSL GPU pair kernel; an agent-native MCP interface; and a force/energy
> conservation contract enforced by property-based tests. Validated against
> LAMMPS on LJ and Buckingham potentials (<0.001% energy agreement after
> convention correction); parallel scaling characterized on MPI and GPU. See
> `ARCHITECTURE.md` for the methods-level description and `bench/` for the
> measured validation and scaling reports.

> **Naming.** The working name is **MODENA**. The alternate brand **RHEA**
> (Rust High-performance Engine for Atomistics) is held in reserve; the crate is
> structured so a rename is mechanical.

---

## Quickstart

```bash
# build (release recommended)
cargo build --release

# run the canonical LJ-argon NVE benchmark (256 atoms, ρ=0.8442)
cargo run --release -- run -i examples/argon_nve.toml

# validate an input file without running
cargo run --release -- check -i examples/argon_nvt.toml

# run the test suite (includes NVE energy-conservation + NVT thermalization)
cargo test
```

Example output (NVE, 256 atoms, 2000 steps):

```
MODENA: finished 2000 steps in 1.07s
  final  KE = 158.83   PE = -1309.07   E_total = -1150.24
  final  T = 0.4136    P = -1.841      W = -1992.64
  perf   533 µs/step  (1875 steps/s)
```

`E_total` stays flat at −1150.24 across the whole run (drift ≈ 2×10⁻⁵).

## Workspace layout

```
modena/
├── crates/
│   ├── modena-core/        # the engine library (potentials, neighbor search,
│   │                       # integrators, thermostats/barostats, electrostatics,
│   │                       # bonded + constraints, parallel decomposition, runner)
│   ├── modena-cli/         # the `modena` binary (TOML input → Simulation::run)
│   ├── modena-mcp/         # agent-native MCP server (Feature A)
│   └── modena-bench/       # Criterion micro-benchmarks
├── examples/               # argon_nve.toml, argon_nvt.toml, ...
├── bench/                  # LAMMPS comparison + scaling + validation reports
├── ARCHITECTURE.md         # methods-level description (paper precision)
└── ROADMAP.md
```

## Architecture (trait seams)

The engine is intentionally small and built around a small set of traits.
Implementing a trait is how you extend MODENA — no core code needs to change.

| Trait | Responsibility | Shipped implementations |
|-------|----------------|--------------------------|
| `PairPotential` | short-range pair U/F | `LennardJones`, `Buckingham`, `Morse`, `SoftSphere`, `Tabulated` — all shifted-force, multi-species (Lorentz–Berthelot + overrides) |
| `ForceModel` | all forces + PE + virial | `PairForceModel` (rayon half-list), `CompositeForceModel` (pair + bonded + Coulomb), `GpuPairForceModel` (`gpu` feature, WGSL) |
| `Integrator` | time stepping | `VelocityVerlet` (NVE, + SHAKEL/RATTLE), `LangevinBAOAB` (NVT), `NoseHooverChain` (NVT) |
| `Thermostat` / `Barostat` | ensemble control | Berendsen, velocity-rescale, Nosé–Hoover chain / Berendsen, Parrinello–Rahman |

The **`ForceModel`** trait is the ML seam: a learned potential (NequIP/MACE
style, or a coarse-grained tabulated force) implements `compute()` and drops
into the existing integrator/thermostat/IO pipeline without touching core code.
A `CompositeForceModel` runs several `ForceModel`s in sequence, each
*accumulating* into the shared force/virial arrays, so a pair model composes
with bonded and electrostatic models.

### Why shifted-force potentials

A truncated potential with only an *energy* shift has a force discontinuity at
the cutoff, which leaks energy every time a pair crosses it. **Every** shipped
pair potential uses the shifted-force form (both U and F → 0 at the cutoff),
which is what lets velocity-Verlet conserve energy cleanly under truncation.
This is also the basis of the **conservation contract** (Feature B): the
reported potential energy is the exact integral of the reported force, so the
NVE Hamiltonian is the true conserved quantity and energy drift is a genuine
numerical signal — verified by property-based tests for every registered
potential. See `ARCHITECTURE.md` §2.

## Input format (TOML)

```toml
units = "reduced"

[structure]
type = "fcc"               # or "xyz" with file = "..." and box_size
lattice_constant = 1.6796
replicas = [4, 4, 4]
element = "Ar"

[model]
type = "lennard_jones"
sigma = 1.0
epsilon = 1.0
cutoff = 2.5

[integrator]
type = "velocity_verlet"   # or "langevin_baoab" (gamma, target_temperature, seed)
dt = 0.005

# Optional, only for velocity_verlet:
# [thermostat]
# type = "berendsen"        # or "velocity_rescaling", "none"
# target_temperature = 1.0
# tau = 0.5

[run]
steps = 2000
initial_temperature = 0.8
seed = 42
trajectory_file = "traj.xyz"
trajectory_every = 100
log_file = "energy.csv"
log_every = 10
```

**Units:** reduced (σ = ε = m = 1, T in ε/kB) or real units (Å, kcal/mol, amu,
K). Multi-species systems use per-type-pair `(σ, ε)` tables with
Lorentz–Berthelot combining and optional explicit overrides.

## Programmatic use

```rust
use modena_core::*;

let mut system = lattice::fcc(1.6796, [4, 4, 4], "Ar");
// ... init velocities ...

let lj = LennardJones::new(1.0, 1.0, 2.5)?;
let model = PairForceModel::new(lj);
let mut neighbors = NeighborList::new(lj.cutoff());
let mut integrator = VelocityVerlet::new();

let mut virial = Virial::default();
model.compute(&mut system, &mut virial, &mut neighbors)?;
for _ in 0..2000 {
    let e = integrator.step(&mut system, &*model, &mut neighbors, 0.005)?;
}
```

## Verification

`cargo test` (~120 tests) runs:
- Analytic energy/force checks for every pair potential (minimum, cutoff-zero,
  repulsive sign, multi-species cross interactions).
- **Conservation contract** (Feature B): property-based tests that
  `F = −∇U` by finite difference, outputs are finite, and the shifted-force
  form is continuous at the cutoff — for *every* registered potential, not just
  LJ.
- **NVE energy conservation** (< 10⁻³ relative drift) and momentum conservation.
- **NVT thermalization** (Langevin BAOAB, Nosé–Hoover) toward the target
  temperature.
- Bonded interactions + SHAKE/RATTLE constraint satisfaction.
- Parallel decomposition: decomposed force/energy/virial **bit-identical to
  serial** (in-process) and across two real MPI ranks (`mpiexec -n 2`).
- GPU pair kernel: single-precision agreement with the f64 CPU path.

**LAMMPS validation** (`bench/`): LJ and Buckingham potentials agree with LAMMPS
to <0.001% on energy after the force-shift/energy-shift convention correction,
with near-identical NVE drift. Parallel scaling (MPI strong/weak, GPU vs CPU) is
reported in `bench/scaling/SCALING.md`; methodology and results in
`bench/BENCHMARKS.md` and `bench/validation/VALIDATION.md`.

## License

MIT.

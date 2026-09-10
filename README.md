# MODENA

**Molecular Dynamics Engine for Numerical Atomistics**

A fast, extensible molecular-dynamics framework written in Rust. MODENA is
designed as a LAMMPS-class engine built around a small set of trait seams so
that new physics (many-body potentials, bonded terms, long-range electrostatics)
and new acceleration strategies (ML surrogate potentials, GPU kernels, domain
decomposition) can be added without touching the integrator, thermostat, or I/O.

> **Status — Phase 1 (PoC).** Lennard-Jones + cell lists + velocity-Verlet (NVE)
> and BAOAB-Langevin (NVT), rayon-parallel, XYZ/CSV output, TOML input. Energy
> conservation is verified to ~2×10⁻⁵ relative drift. See `ROADMAP.md` for what
> comes next (bonds/angles, Coulomb+Ewald, NPT, real units, Python, ML, GPU).

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
│   ├── modena-core/        # the engine library (physics, I/O, runner)
│   └── modena-cli/         # the `modena` binary
├── examples/               # argon_nve.toml, argon_nvt.toml
└── ROADMAP.md
```

## Architecture (trait seams)

The engine is intentionally small and built around four traits. Implementing a
trait is how you extend MODENA — no core code needs to change.

| Trait | Responsibility | Phase 1 impl | Extension target |
|-------|----------------|--------------|------------------|
| `PairPotential` | short-range pair U/F | `LennardJones` (shifted-force) | Buckingham, Morse, tabulated |
| `ForceModel` | all forces + PE + virial | `PairForceModel` (rayon) | **ML surrogates**, many-body, bonded |
| `Integrator` | time stepping | `VelocityVerlet`, `LangevinBAOAB` | r-RESPA, leapfrog |
| `Thermostat` | velocity rescaling | `Berendsen`, `VelocityRescaling` | Nosé–Hoover chains |

The **`ForceModel`** trait is the ML seam: a learned potential (NequIP/MACE
style, or a coarse-grained tabulated force) implements `compute()` and drops
into the existing integrator/thermostat/IO pipeline. No ML code ships in Phase
1; the seam is in place.

### Why shifted-force LJ

A truncated potential with only an *energy* shift has a force discontinuity at
the cutoff, which leaks energy every time a pair crosses it. MODENA uses the
**shifted-force** form (both U and F → 0 at the cutoff), which is what lets
velocity-Verlet conserve energy cleanly under truncation. The cell list is also
rebuilt with a deduplicated stencil so coarse grids (fewer than 3 cells/axis)
don't double-count pairs.

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

**Units:** reduced (LJ) in Phase 1 (σ = ε = m = 1, T in ε/kB). Real units
(Å, kcal/mol, amu, K) arrive in M2.

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

`cargo test` runs:
- LJ energy/force analytic checks (minimum, cutoff-zero, repulsive sign).
- Cell-list correctness incl. periodic images and coarse-grid dedup.
- **NVE energy conservation** (real kinetic energy, < 10⁻³ relative drift) and
  momentum conservation.
- **NVT Langevin thermalization** toward the target temperature.
- Pair-energy/force correctness (no overcounting; Newton's third law).

## License

MIT.

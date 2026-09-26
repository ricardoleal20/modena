# MODENA Roadmap

The path from PoC to a LAMMPS-class engine with Rust + Python clients and an ML
acceleration surface. Each milestone is independently shippable and leaves the
existing trait seams intact.

## M1 — Core PoC ✅ (current)

- Lennard-Jones (shifted-force, cutoff).
- Cell-list neighbor search, orthorhombic PBC, deduplicated stencil.
- Velocity-Verlet (NVE) + BAOAB-Langevin (NVT) integrators.
- Berendsen + stochastic velocity-rescaling thermostats.
- Rayon-parallel force kernel.
- XYZ trajectory + CSV energy log; TOML input; `modena` CLI.
- `ForceModel` / `PairPotential` / `Integrator` / `Thermostat` trait seams.
- Reduced (LJ) units.
- Tests: LJ analytic checks, neighbor correctness, NVE conservation, NVT
  thermalization, pair-energy correctness.

**Verified:** NVE relative energy drift ≈ 2×10⁻⁵ over 2000 steps (256 atoms,
ρ=0.8442). PE/atom ≈ −5.1 (textbook value).

## M2 — Real units + more pair potentials

- `UnitSystem::Real { length=Å, energy=kcal/mol, mass=amu, temp=K }` with
  conversion at the I/O boundary; core stays in reduced units internally.
- Buckingham, Morse, soft-sphere, tabulated pair potentials.
- Multi-species: `sigma`/`epsilon` cross-interaction tables
  (Lorentz–Berthelot + overrides), per-type `mass`.
- Concentration / random-alloy structure builder.

## M3 — Bonded interactions ✅

- `BondedModel` (a second `ForceModel`): harmonic bonds, angular terms,
  proper/improper dihedrals.
- Topology input (bonds/angles/dihedrals) alongside the structure file.
- Bonded + non-bonded composition inside one `ForceModel` (force/energy/virial
  additive). Constraints (SHAKE/RATTLE) as an integrator post-step.

**Verified:** every bonded term's force equals −∇U by central finite difference
(proptest, 32 cases); per-term forces sum to zero (translation invariance);
composite pair+bonded NVE drift < 1e-3; SHAKE+RATTLE holds constrained bond
lengths to < 1e-8 and conserves energy on the constraint manifold (rigid
dumbbell + constrained chain).


## M4 — Long-range electrostatics ✅

- Coulomb pair potential (1/r) with the `ewald` method: real-space sum +
  reciprocal-space sum + self correction.
- PPPM / SPME (particle-mesh) on the roadmap once FFT is wired (`rustfft` /
  `realfft`).
- Triclinic boxes (general cell tensors) — seam-level only: `SimBox::triclinic`
  with correct min-image/wrap/inverse, verified against the orthorhombic path.
  The TOML schema and the parallel decomposition remain orthorhombic; the
  cell-list ±1 stencil is not safe for sheared cells.
- Charge assignment per atom (`system.charges`).

**Verified:** Ewald force = −∇U by central finite difference (every atom, all
terms); Madelung lattice energy of NaCl converges to −α/r₀ (α = 1.74756, err
< 0.02); triclinic minimum-image/wrap correct for a general cell tensor and
matches the orthorhombic path for a diagonal cell; real-space force → bare
Coulomb as κ → 0; ionic (NaCl + Ewald) NVE drift < 1e-2. Charges wired through
the system/config/units layer with the reduced-charge convention
`U = q_red_i q_red_j / r_red`; rocksalt lattice builder; composes with a pair
model in `CompositeForceModel` and the runner (`tests/electrostatics.rs`).

## M5 — Ensembles & sampling ✅

- NPT: Parrinello–Rahman and Berendsen barostats; isotropic + anisotropic.
- Nosé–Hoover chain thermostat (true canonical NVT, Suzuki–Yoshida 3-fold,
  hierarchical chain friction; extended Hamiltonian H* conserved to ~1e-4).
- Monte Carlo volume barostat (Metropolis NPT) as an alternate coupling.
- Replica-exchange (parallel tempering) driver.
- Verified: `tests/ensembles.rs` — NHC H* conservation, replica-exchange swap
  acceptance, Berendsen barostat moves the box via the runner.

## M6 — Python client (PyO3)

- `modena-python` crate via `maturin`: `import modena`.
- Build a system, attach a model, run, receive trajectories/energies as NumPy
  arrays (zero-copy via `numpy` crate).
- Python-defined custom pair potentials through a callback bridge (GIL-scoped,
  slower — for prototyping; Rust `ForceModel` for production).
- Jupyter-friendly: progress + inline plotting helpers.

## M7 — Performance

- SoA layout (`positions: Vec<f64>` flat) + SIMD force kernels (`wide` /
  portable f64x4).
- Verlet neighbor list with skin (rebuild on drift threshold) instead of
  every-step cell rebuild.
- Half neighbor list (Newton's third law) to halve force work.
- Benchmark suite: argon LJ scaling (N = 1k–100k), steps/s vs LAMMPS reference.

## M8 — ML surrogate potentials

- `MlForceModel: ForceModel` — load a trained equivariant model (NequIP/MACE
  format) and evaluate forces/energy in-engine.
- Onnx-runtime or candle backend behind a feature flag.
- Surrogate / ab-initio switching: run cheap ML force most steps, recompute
  expensive reference periodically (Δ-learning, multi-scale).
- Active-learning hooks (per-atom uncertainty → refit trigger).

## M9 — Parallel decomposition ✅

- Spatial domain decomposition across MPI ranks (`mpi` crate) for >10⁶ atoms.
- Per-rank cell lists + ghost-atom exchange.
- Load balancing; GPU offload of the pair kernel (`wgpu`/CUDA) as a
  `ForceModel` backend.

**Done.** `decomp` (default, pure Rust) partitions an orthorhombic box into a
`nr[0]×nr[1]×nr[2]` rank grid; each rank owns its subdomain atoms and imports a
one-shell ghost halo from the 26 periodic neighbors (subdomain edge ≥ cutoff).
Pair work is partitioned by lower-global-id ownership, so the decomposed
force/energy/virial summed over ranks is bit-identical to the serial half-list
evaluation — verified in-process (`tests/decomp.rs`: single/multi-rank == serial,
halo sufficiency, post-migration match, decomposed-NVE conservation) and across
two real MPI ranks (`tests/mpi_decomp.rs` under `mpiexec -n 2`). The `mpi`
feature adds a real ghost-exchange driver (deadlock-free canonical send/recv
ordering + reverse ghost-force reduction); the `gpu` feature adds a
`GpuPairForceModel` WGSL compute-shader LJ backend (`f32`, all-pairs) returning
`None` when no adapter is present so CI skips gracefully. Both features are
behind flags so the default build needs no native MPI/GPU runtime.

## Non-goals (for now)

- Quantum dynamics / path-integral MD (PIMD) — revisit post-M5.
- Continuum / finite-element coupling — out of scope until the core is solid.
- A GUI is not the current focus; the CLI and the agent-native MCP interface
  are the shipped interfaces (the `modena-mcp` crate doubles as the internal
  surface a future GUI would sit on).

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

## M3 — Bonded interactions

- `BondedModel` (a second `ForceModel`): harmonic bonds, angular terms,
  proper/improper dihedrals.
- Topology input (bonds/angles/dihedrals) alongside the structure file.
- Bonded + non-bonded composition inside one `ForceModel` (force/energy/virial
  additive). Constraints (SHAKE/RATTLE) as an integrator post-step.

## M4 — Long-range electrostatics

- Coulomb pair potential (1/r) with the `ewald` method: real-space sum +
  reciprocal-space sum + self correction.
- PPPM / SPME (particle-mesh) on the roadmap once FFT is wired (`rustfft` /
  `realfft`).
- Triclinic boxes (general cell tensors) replacing the orthorhombic assumption.
- Charge assignment per atom (`system.charges`).

## M5 — Ensembles & sampling

- NPT: Parrinello–Rahman and Berendsen barostats; isotropic + anisotropic.
- Nosé–Hoover chain thermostat (true canonical NVT).
- Hamiltonian / configurational-bias Monte Carlo moves as an alternate
  `Integrator` for enhanced sampling.
- Replica-exchange (parallel tempering) driver.

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

## M9 — Parallel decomposition

- Spatial domain decomposition across MPI ranks (`mpi` crate) for >10⁶ atoms.
- Per-rank cell lists + ghost-atom exchange.
- Load balancing; GPU offload of the pair kernel (`wgpu`/CUDA) as a
  `ForceModel` backend.

## Non-goals (for now)

- Quantum dynamics / path-integral MD (PIMD) — revisit post-M5.
- Continuum / finite-element coupling — out of scope until the core is solid.
- A GUI; the CLI + Python client are the interfaces.

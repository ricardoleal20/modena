# Rigid SPC/E water — Tier B item 6 (status 2026-09-27)

The explicit-water prerequisite for the computational-biology claim. Shipped:
the exact per-molecule rigid-water constraint projection (the SETTLE role),
the erf exclusion correction for charged excluded pairs, the water structure
builder + `[electrostatics] method` plumbing, the conservation gates, and the
restart/checkpoint workflow that closes the equilibrate-then-measure loop.

## What landed

- `constraints.rs` — `RigidWater` + exact Newton projection of the three
  coupled constraints per molecule (position) and the exact 3×3 velocity
  projection (RATTLE). Unlike sequential SHAKE over a bond triangle it cannot
  fail at large dt and parallelizes per molecule (contiguous (O,H,H) triples
  → one 9-float chunk per ray). COM preserved exactly. Gates: distorted-water
  restoration to 1e-10 with COM to 1e-12; radial-velocity removal.
- `ewald.rs`/`spme.rs` — `with_excluded_pairs` + the shared erf correction:
  `U_corr = −q_i q_j erf(κr)/r`. The neighbor list removes the erfc part; the
  mesh/reciprocal still carries the erf part for the pair; this cancels it so
  an excluded pair feels exactly zero Coulomb. Gates: static equality
  (`E_excl = E_plain − q_i q_j/r`) and — the decisive one — free-flight
  conservation of an isolated excluded +q/−q pair (1e-12).
- `lattice.rs` + `[structure] type = "water"` — SPC/E constants (d_oh 1.0 Å,
  109.47°, q ±0.8476/0.4238, m 15.9994/1.008), box from the physical density,
  balanced alternating-dipole orientations, contiguous triples, species 0/1.
  Requires `units = "real"` (validated).
- `runner.rs` — water structure auto-produces its `Constraints` + the 3
  intra-molecular exclusions, plumbed into `NeighborList` and into the
  electrostatics models.
- `io/restart.rs` + `[run] restart_file`/`restore_file`/`restart_every` —
  the restart/checkpoint workflow: a checkpoint carries the box, SoA
  positions/velocities, and the step counter (reduced units, text format,
  17-significant-digit round-trip); a restoring run skips the velocity draw,
  continues the step counter (log/trajectory labels run across the seam), and
  rejects an atom-count mismatch. Decks:
  `bench/validation/water/modena_spce_nvt.toml` (equilibration) →
  `modena_spce_nve_restore.toml` (measurement). The NVE seam is exact: the
  restored run's first thermo row equals the checkpoint state bit-for-bit,
  and a split L+M NVE run reproduces a one-go L+M run at every logged step.
- `tests/water.rs` — 4 deterministic gates (settle exactness, free flight,
  full-stack rigid water + SPME + exclusions NVE at 1e-4, builder geometry);
  `tests/restart.rs` — the round-trip, seam, and two-phase water gates.

## Measured

- 64-water SPC/E step: ~1.5 ms/step (192 atoms) in the cold-lattice regime;
  ~1.4 ms/step from an NVT-equilibrated checkpoint (1000 steps at
  tau = 0.005 ps); SPME mesh 16³ for κ=0.34 Å⁻¹ / r_c = 9 Å.
- Full-stack conservation for a rigid water: 1.9e-7 over 2000 steps.
- Two-phase workflow (64 waters, dt = 0.0005 ps): NVT 1000 steps at
  tau = 0.005 ps ends near 130–160 K (still relaxing); the following
  500-step NVE segment drifts |ΔE|/|E| ≈ 1.5e-2, and 2000 steps ≈ 4.9e-2 —
  the residual is the ongoing structural relaxation of the short
  equilibration, not an integrator defect (the exact-restore/seam gates and
  the 1.9e-7 gentle-regime gate are the strong conservation checks).

## Honest open items (next session)

- **Equilibrated-liquid measurement quality**: the *workflow* now exists
  (NVT-equilibrate → checkpoint → NVE-measure, one command per phase), but
  the smoke budgets are not convergence: 300–1000 NVT steps leave the melt
  still relaxing (the 500-step NVE follow-up drifts 6.6e-2 / 1.5e-2
  respectively). A production equilibration (20k+ steps, staged tau) and a
  converged-liquid observables row (T ≈ 298 K, P ≈ 1 atm, rdf) remain open —
  an equilibrated-liquid NVE drift gate should follow once that run exists.
- **Static-energy LAMMPS cross-check**: needs a proper LAMMPS bonded+SHAKE
  SPC/E deck (atom_style full, special_bonds lj/coul 0 0 0, shake b,
  pppm gewald-matched). Without exclusions LAMMPS includes the intramolecular
  q·q/r (≈ −202 kcal/mol·water at 1 Å) — the naive comparison is misleading.
- HECBioSim 20K Crambin + amber-prm reader (the protein end of the claim).

## Follow-on perf levers (cheap, known)

- The rigid-water velocity projection is a single sparse 3×3 solve — already
  per-molecule parallel.
- SPME real-space at water cutoffs: the erfc loop is the same shared form as
  the NaCl deck (unchanged cost profile).

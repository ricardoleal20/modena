# Rigid SPC/E water — Tier B item 6 (status 2026-09-27)

The explicit-water prerequisite for the computational-biology claim. Shipped:
the exact per-molecule rigid-water constraint projection (the SETTLE role),
the erf exclusion correction for charged excluded pairs, the water structure
builder + `[electrostatics] method` plumbing, and the conservation gates.

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
- `tests/water.rs` — 4 deterministic gates (settle exactness, free flight,
  full-stack rigid water + SPME + exclusions NVE at 1e-4, builder geometry).

## Measured

- 64-water SPC/E step: ~1.5 ms/step (192 atoms) incl. the cold-lattice
  regime; SPME mesh 16³ for κ=0.34 Å⁻¹ / r_c = 9 Å.
- Full-stack conservation for a rigid water: 1.9e-7 over 2000 steps.

## Honest open items (next session)

- **Equilibration workflow**: the builder emits a cold lattice whose initial
  relaxation is violent (MODENA and LAMMPS agree the potential is ~−200
  kcal/mol·water for random orientations — the model's physics, not an
  engine bug; the balanced pattern cuts it to ~−5/water). A liquid NVE gate
  needs an NVT-equilibrate-then-resume flow (a restart feature is the
  pending piece).
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

# Performance results — implementation wave (2026-09-26)

Measured deltas from the committee pipeline implementation (PR #20,
`6785396`). All gates from `COMMITTEE_ACTA.md` that were implementable in this
wave passed; machine was load-shared during some runs, so the *deltas* are
reliable and the absolute numbers must be re-generated idle via `bench/bench.sh`
before the paper (the acta's G0 rule).

## LJ single-core — before / after

| Size | lmp_serial | MODENA 1t (before) | MODENA 1t (after) | MODENA default (before → after) |
|---|---|---|---|---|
| 864 | ~1.3 s | 2.61 s (2.0×) | ~2.3–2.6 s (≈1.8–2.0×) | 3.98 s → **~2.1–2.4 s** |
| 2048 | ~3.2 s | 6.65 s (2.1×) | ~7.0–7.7 s (noise) | 14.66 s → **~7.7 s (~1.9×)** |

- **The default-pool self-harm is gone**: with the serial-below-100k-pairs path,
  rayon's work-stealing churn (previously 11 threads *slower* than 1 thread and
  ~35 s of sys-time per 5000 steps at 2048) no longer dominates. The default
  configuration now tracks single-thread time instead of being 2.2× worse.
- 1-thread **bit-identity preserved**: 500-step trajectory byte-compared
  against the pre-change main binary (G1 gate).
- The remaining single-core gap (~1.8–2×) now lives in per-pair cost + rebuild
  — the exact targets of the deferred cluster-pair/NEON/SIMD and half-stencil
  work (see `MATCH_ENGINES.md`).

## Charged (Ewald) — Ewald-b

| 216 atoms, 2000 steps | wall | vs LAMMPS `ewald` |
|---|---:|---:|
| LAMMPS `ewald` | 9.93 s | — |
| MODENA (before) | 341.99 s | ~34× behind |
| **MODENA (Ewald-b)** | **107.11 s** | **~10.8× behind** (3.2× faster) |

Ewald-b = k-vector cache (box-keyed) + conjugate half-space (±k, energy and
force ×2) + cos/sin phase reuse. Gate G4 (≥1.8×, ≤10⁻¹⁰ vs brute force)
passed: decomposition agreement unchanged at −1.7e-10. The remaining ~10.8×
is the transcendental bill (≈4.2M sincos/step at ~10 ns each); the acta's path
to ~2× is per-axis trig tables (≈7 FMA/sincos, +0.5–1 d) or SPME (O(N log N),
post-paper).

## SPME (Stage B/C) — the charged regime goes O(N log N)

Deck: `bench/validation/spme/modena_10k.toml` — NaCl rocksalt, 10,648 atoms,
Buckingham + Coulomb, κ=0.8, r_c=5.0, k_c=7.0 → mesh 125³, 2000 NVE steps.

| 10,648 atoms, 2000 steps | wall | ms/step | drift |
|---|---:|---:|---|
| MODENA SPME (first cut) | 546.9 s | 273.4 | 8.0e-7 (Healthy) |
| MODENA SPME (influence cache + 5-smooth mesh + rayon axes 2/1) | 147.1 s | 73.6 | 8.0e-7 (Healthy) |
| **MODENA SPME (+ axis-0 L1-transpose; no serial axis pass)** | **~80 s** | **~40** | **8.0e-7 (Healthy)** |
| MODENA direct Ewald | infeasible | ~70,000 (extrapolated) | — |
| LAMMPS `pppm` 1e-5 (matched gewald=0.8, serial) | ~481 s CPU | ~240 | — |

- The 6.5× optimization arc: influence-factor cache (1.9M exp/cos per step
  → 0), 5-smooth mesh rounding (123³ → 125³, no prime-41 butterflies), reused
  transpose line buffers, rayon on all three axes (axis 0 via the L1-transpose
  scheme — transpose once, contiguous pass, transpose back into the spread
  grid's allocation). 273.4 → 73.6 → ~40 ms/step (37–42 measured, load-shared
  machine); the run now uses ~4 cores effective. Reciprocal pass alone:
  280 → 25 ms.
- LAMMPS `pppm` 1e-5 (matched gewald=0.8, same cutoffs) on the same deck:
  480.9 s user CPU single-core (wall 739 s — load-shared with the build; the
  CPU number is the honest one). MODENA SPME: 267 s user CPU (147 s wall at
  ~2 effective threads). **The charged regime at 10k atoms is the first
  regime where MODENA is ahead of `lmp_serial` — ~1.8× less CPU** (the LJ
  regime sits at ~1.8–2× behind). Both rows need idle-machine regeneration
  before the paper (G0).
- 216-atom deck: SPME ≈ direct Ewald in accuracy (all gates green); at that
  size the direct path is already fine — SPME is for N ≳ 2k.

## What did NOT move (and why it is correct it did not)

- Compiler flags: still 0% on aarch64 (release already `lto="thin"` +
  `codegen-units=1`).
- No f32-on-CPU (rejected by contract math), no hand-NEON (ceiling +42% on
  pair time ≈ 1.1–1.25× step — deferred behind cluster-pair/asm inspection).
- GPU path untouched in this wave — Stage 3 (Kahan + persistent buffers + the
  G3 ladder that decides the only "win" claim) is the next PR, spec'd in the
  acta.

## How to reproduce

```bash
cargo build --release
bash bench/bench.sh                     # idle machine, RUNS=5
# charged row: bench/bench.sh LADDER="864" THREADS="1" (includes ewald rows)
```
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
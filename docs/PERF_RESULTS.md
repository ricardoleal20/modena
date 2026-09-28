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

## Rigid SPC/E water — Tier B item 6 pathfinder

Deck: `bench/validation/water/modena_spc_e.toml` — 64 rigid SPC/E waters
(192 atoms), SPME (κ=0.34 Å⁻¹, r_c=9 Å, mesh 16³), intra-molecular
exclusions, exact per-molecule rigid projection.

| 64 waters (192 atoms) | per step |
|---|---:|
| MODENA SPME (cold-lattice regime) | ~1.5 ms (~7.8 µs/(atom·step)) |
| MODENA SPME (NVT-equilibrated, via restart) | **~1.4 ms (~7.3 µs/(atom·step))** |

The equilibrated row is measured through the restart workflow that shipped
with it: `modena_spce_nvt.toml` (1000 NVT steps, tau = 0.005 ps) writes a
checkpoint, `modena_spce_nve_restore.toml` restores it and runs 2000 NVE
steps at ~1.36–1.45 ms/step (load-shared machine, load avg ~10/11 cores —
regenerate idle before the paper). The equilibrated melt relaxes less
violently than the cold lattice, so the step time settles slightly below the
cold-regime row. Caveats from `docs/WATER_MODEL.md` still apply: 1000 NVT
steps is a smoke budget, not convergence — the 2000-step NVE follow-up
drifts |ΔE|/|E| ≈ 4.9e-2 from the ongoing structural relaxation. The
conservation gates that MATTER — rigid water + SPME + exclusions NVE at
1.9e-7 over 2000 steps, excluded-pair free flight at 1e-12, and the
restart seam gates (bit-exact restore; split L+M NVE ≡ one-go L+M) — are
committed tests and hold exactly.

## GPU-d decision ladder (gate G3) — the GPU claim, resolved

Protocol (per `docs/COMMITTEE_ACTA.md` G3 and the Hoefler&Belli design in
`MATCH_ENGINES.md`): matched LJ-NVE decks (FCC, ρ=0.8442, shifted-force LJ,
dt=0.005, seed 87287; 5000/10000/20000-atom decks use rectangular FCC blocks
— 5×10×25, 10×10×25, 10×20×25 cells — to hit the ladder sizes exactly),
**median (min) ms per integration step over 5 measured runs after warm-up
discard**, fence-synced on the GPU paths (every step blocks on the
staging-buffer map — `device.poll(PollType::wait)` — before the clock reads).
Engines: `lmp_serial` 22Jul2025 (serial, its own decks), MODENA CPU (release
CLI on the same decks, default threads — the shipped configuration), MODENA
GPU all-pairs (`GpuGather::AllPairs`, the pre-G3 kernel, kept as the
determinism oracle) and MODENA GPU cell-list (`GpuGather::CellList`, the G3
remedy). Driver: `examples/gpu_ladder.rs`; rig: `bench/gpu_ladder.sh`.

**Honesty caveats, same table**: the GPU kernels compute in **f32**
(Kahan-compensated in-kernel, f64 host reduction; WGSL has no f64 on Metal),
and their neighbor bins are **built on the CPU** each rebuild — the LAMMPS
GPU-package `neigh no` mode (lammps#3531), not a GPU-resident build. The
machine was **load-shared** during the runs (load average 4–11: a concurrent
worktree's test binary plus desktop composition on the shared GPU); medians
and minima are reported for that reason, and every number below must be
regenerated on an idle machine before the paper (the acta's G0 rule). GPU
rows used 1000/1000/500/250/100 steps per run at 864/2048/5000/10000/20000
(the O(N²) kernel is the wall-clock constraint); CPU/LAMMPS rows ran their
decks' 5000/5000/2000/2000/2000 steps.

| atoms | lmp_serial (f64) | MODENA CPU default (f64) | GPU all-pairs (f32) | GPU cell-list (f32) | gather probe (f32) |
|---:|---:|---:|---:|---:|---:|
| 864 | 0.292 (0.290) | 0.454 (0.438) | 1.73 (1.68) | 1.78 (1.62) | 1.65 (1.57) |
| 2048 | 0.800 (0.756) | 1.63 (1.53) | 2.66 (2.50) | 2.41 (2.39) | 1.95 (1.78) |
| 5000 | 2.37 (2.23) | 1.94 (1.85) | 2.75 (1.69) | 2.53 (2.52) | 2.02 (1.95) |
| 10000 | 3.56 (3.51) | 3.98 (3.94) | 4.38 (4.29) | 4.05 (3.95) | 4.57 (4.24) |
| 20000 | 7.94 (7.52) | 8.19 (7.96) | 15.39 (15.27) | **7.90 (7.89)** | 15.52 (15.41) |

At 20 000 atoms that is 0.40 µs/(atom·step) for `lmp_serial`, 0.41 for MODENA
CPU, 0.40 for the GPU cell-list — a three-way statistical parity (run maxima
overlap: 8.54 / 9.89 / 8.65), not a win.

### Why the ladder looks the way it does (phase decomposition)

A static-system probe (100 consecutive `compute` calls, same load regime)
separates the three costs the ladder mixes together:

| n | CPU neighbor rebuild | GPU cell-list steady | cell-list w/ forced rebuild | GPU all-pairs steady |
|---:|---:|---:|---:|---:|
| 256 | 0.13 ms | 1.61 ms | 2.18 ms | 1.82 ms |
| 2048 | 4.8 ms | 1.55 ms | 7.6 ms | 1.93 ms |
| 20000 | 42.4 ms | 1.83 ms | 47.3 ms | 13.6 ms |

- **Fence-sync round-trip floor ≈ 1.5–1.6 ms/step** (the 256-atom steady
  column): submit + copy + map + readback on a desktop-composited, shared
  GPU. This floor — not the kernel — dominates every GPU row below ~10k
  atoms, which is why all four GPU engines cluster at 864–10000 in the
  ladder.
- **The cell-list kernel works**: at 20k atoms the steady-state GPU step is
  1.83 ms vs 13.6 ms all-pairs — **7.5×** at the kernel level. The gather
  probe (all-pairs j-loop, LJ math skipped) matches the full all-pairs time
  at every size, i.e. the j-loop gather + minimum-image + distance work is
  ~90–95% of the all-pairs kernel wherever the kernel is resolvable above
  the floor — exactly the >30% G3 trigger that mandated the cell-list.
- **The remaining 20k bottleneck is CPU-side**: the serial neighbor rebuild
  (42 ms per rebuild at 20k, ~54 ns/pair — a shared MODENA cost the CPU row
  pays too, amortized ~4–5 ms/step) plus the CSR build + upload (~5 ms per
  rebuild). LAMMPS rebuilds the same-size list 224 times in 2000 steps and
  stays at 7.9 ms/step total, so its rebuild is sub-millisecond — the
  rebuild, not the GPU, is the next target (acta Tier A, `rebuild_count`
  hygiene).

### Drift classification (G3 drift gate)

2000-step NVE at seed 87287, same initial state, `ConservationReport` on the
final-vs-initial total energy (band = (E_max−E_min)/|E_0| over all steps):

| engine | 864 atoms | 2048 atoms |
|---|---|---|
| MODENA CPU (f64) | 1.24e-5, band 2.89e-4 — **Healthy** | 1.26e-5, band 2.83e-4 — **Healthy** |
| GPU all-pairs (f32+Kahan) | 1.05e-5, band 2.89e-4 — **Healthy** | 1.53e-5, band 2.83e-4 — **Healthy** |
| GPU cell-list (f32+Kahan) | 1.05e-5, band 2.89e-4 — **Healthy** | 1.17e-5, band 2.83e-4 — **Healthy** |

The f32+Kahan GPU trajectories conserve the NVE Hamiltonian to the same
order as the f64 CPU run at the same seed — the drift gate **passes** for
both kernels (Kahan-under-no-f64 is doing its job; WGSL has no f64 on
Metal).

### G3 verdict

**The GPU "win" claim does not survive G3 as a headline.** The pre-G3
all-pairs kernel loses to `lmp_serial` and to MODENA CPU at every ladder
size (1.7–15.4 ms/step vs 0.29–7.94). The gather share it wastes is the
measured ~90–95% of kernel time, so the cell-list gather (CPU Verlet bins →
contiguous per-atom j-ranges → one GPU thread per i, no atomics) was
implemented — it is 7.5× at the kernel level and lifts the 20 000-atom row
to statistical parity with `lmp_serial` (7.90 vs 7.94 ms/step, overlapping
spreads) — but parity at the top of the ladder is not a win, and below
~10k atoms the fence-synced round-trip floor (~1.5 ms/step on this
composited GPU) keeps every GPU configuration behind the CPU. The honest
paper wording, with caveats in the same table: *f32 (Kahan) kernel,
CPU-built bins (`neigh no` mode), fence-synced timing, load-shared machine
— GPU reaches parity with `lmp_serial` at 20k LJ atoms and is behind at all
smaller sizes; the binding constraints are the per-step host↔GPU round trip
and the serial CPU neighbor rebuild, not the pair kernel.* The path to an
actual win is identified and deferred: an async double-buffered force
pipeline to hide the round-trip (changes the per-step sync contract), and a
rebuilt neighbor search (which lifts the CPU row too). The all-pairs kernel
is retained as the determinism oracle (`tests/gpu_pair.rs` compares
cell-list vs all-pairs vs serial f64); both kernels hold Healthy drift.

## Parallel neighbor rebuild — the 20k bottleneck, retired

The GPU-d ladder identified the **serial CPU neighbor rebuild** (42 ms/rebuild
at 20k, ~54 ns/pair) as the shared bottleneck of the large rows. The rebuild
now parallelizes over `i`-ranges with rayon (`par_chunks().map().collect()`
preserves the per-chunk order and the chunk sequence, and the per-`i` scan is
unchanged, so the pair list is **byte-identical to the serial build** — the
force-accumulation order, and hence runnable determinism, is untouched; a
`MODENA_SERIAL_NEIGHBOR=1` env knob forces the old path for bench A/B).

Rebuild timing (release, load-shared machine, one core vs the rayon pool):

| n | serial | parallel | speedup |
|---:|---:|---:|---:|
| 10 000 | 12.0 ms | 3.06 ms | 3.9× |
| 20 000 | 42.5 ms | 8.58 ms | 4.9× |
| 40 000 | 204 ms | 40.9 ms | 5.0× |

Ladder 20 000-atom row, re-measured with the parallel rebuild (5 runs, median
ms/step, full VV steps, same decks/seed — load-shared machine, G0 still
mandatory before the paper):

| 20000 atoms | before (G3 ladder) | after parallel rebuild |
|---:|---:|---:|
| MODENA CPU default (f64) | 8.19 | **4.05** |
| MODENA GPU cell-list (f32) | 7.90 | **3.69** |

Both rows roughly halve, drift stays Healthy on both paths (~1.7e-5), and
the GPU cell-list now edges past the CPU row at 20k within spread. The
fence-sync floor (~1.5 ms) still caps the GPU's advantage, and below ~10k the
serial path (n < 4096) means small-system behavior is unchanged.

## What did NOT move (and why it is correct it did not)

- Compiler flags: still 0% on aarch64 (release already `lto="thin"` +
  `codegen-units=1`).
- No f32-on-CPU (rejected by contract math), no hand-NEON (ceiling +42% on
  pair time ≈ 1.1–1.25× step — deferred behind cluster-pair/asm inspection).
- No async GPU force pipeline (would hide the ~1.5 ms round-trip floor but
  breaks the fence-synced per-step contract G3 measures under; deferred with
  the G3 verdict above).

## How to reproduce

```bash
cargo build --release
bash bench/bench.sh                     # idle machine, RUNS=5
# charged row: bench/bench.sh LADDER="864" THREADS="1" (includes ewald rows)
# GPU-d ladder (G3): build the gpu example first
cargo build --release --features gpu --example gpu_ladder
bash bench/gpu_ladder.sh                # idle machine, RUNS=5
# drift rows (2000-step NVE, CPU vs both GPU kernels, same seed):
target/release/examples/gpu_ladder gpu 2048 2000 1
target/release/examples/gpu_ladder gpu-csr 2048 2000 1
target/release/examples/gpu_ladder cpu 2048 2000 1
```
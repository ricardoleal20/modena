# Performance roadmap — approaching LAMMPS (and where we can win)

Research date: 2026-09-26. Method: two-agent deep research (bibliography of
LAMMPS optimizations + adversarial verification), synthesized below. Baseline:
MODENA is ~2.3–4.5× behind `lmp_serial` single-core on LJ at 256–4000 atoms
(M3 Pro, measured 2026-09-26, `bench/run_comparison.sh`).

## Executive summary (honest envelope)

- **Single-core CPU: closing 4.5× → ~1.5–2× behind LAMMPS is achievable before
  the paper** (waste-removal × compiler flags × scalar/SIMD kernel work ≈
  1.8–2.3× multiplicative). **Single-core parity is not a realistic first-paper
  target**, and chasing it invites a reviewer to run one benchmark and reject
  the paper.
- **The defensible "we win" claim is the Apple GPU.** LAMMPS has **no Metal /
  Apple-GPU backend** (GPU package = CUDA/HIP/OpenCL [1]; OpenCL broken on Apple
  silicon — lammps#3531 [2]; Kokkos has no Apple backend — kokkos#8461 [3]).
  MODENA's wgpu/Metal path is the *only* head-to-head where LAMMPS forfeits by
  absence. Measure it at ≥10k atoms; publish that.
- **Framing for the paper**: "correct-by-construction Rust engine, measured
  performance gap closed from ~4.5× to ~1.5–2×, and the only engine with a
  working Apple-GPU path." That matches what the 2025–26 literature rewards
  (JAX-MD published while slower — differentiability; TorchSim published
  conceding classical-speed — MLIP ergonomics [4][5]).

## What makes LAMMPS fast (evidence)

| # | Technique | Source | Gain |
|---|---|---|---|
| 1 | Spatial-decomposition + half neighbor lists (Newton's 3rd law, each pair once) | Plimpton 1995, JCP 117, 1 [6]; Thompson et al. CPC 271, 108171 (2022) [7] | foundational — **MODENA already has this** |
| 2 | Atoms sorted by spatial locality (contiguous per-atom data → cache) | LAMMPS "Measuring performance" docs [8] | qualitative serial factor |
| 3 | Skin-vs-rebuild frequency tuning | LAMMPS docs [8] | rebuild less often |
| 4 | **Branch-free templated pair loop** | LAMMPS OPT package [9] | **5–20% on Pair time** |
| 5 | **Explicit SIMD inner loop** (INTEL package: mixed precision f32-pair/f64-accumulate, `newton off`, sorted bins) | LAMMPS INTEL docs [10]; Watanabe & Nakagawa arXiv:1806.05713 [11] | LN: **+42% Haswell vs optimized non-vectorized; +34% vs auto-vec on KNL**; AVX2 ≈ AVX-512 on Skylake |
| 6 | Single-pass energy+force with shared r⁻⁶ intermediates | dreid-kernel (Rust) benchmarks [12] | combined 1.19 ns vs ~0.94+0.97 ns separate (~40% cheaper) |
| 7 | GROMACS MxN cluster-pair kernel (branchless, SIMD-width clusters) | MD-Bench, FGCS 2023, arXiv:2302.14660 [13] | documented CPU state of the art (reference design) |

Reference target: LAMMPS INTEL test-suite rates (mix-precision vectorized LJ)
≈ **0.12–0.17 µs/(atom·step·core)** on 2015–17 Xeons [14], derived from
199.5–317.3 steps/s @ 512k atoms on 18–20 cores. **Aim within ~2× of that on
the M3 Pro before hand-SIMD** (LLVM may already vectorize the SoA loop).

## Prioritized roadmap

### P0 — Delete pure waste (½–1 day, zero contract risk, multiplicative)

1. **Ewald k-vectors: regenerate only on box change, reuse the buffer.** The
   k-space grid depends only on box + cutoffs; it is re-enumerated and
   re-allocated every step (audit-flagged). Hoist it into the `EwaldForceModel`
   state.
2. **One neighbor-list build per composite step.** The audit found
   `PairForceModel` and `EwaldForceModel` each call `neighbors.build()`; the
   second is a no-op only when the skin absorbed it — pass the same built list
   to both kernels.
3. **Single-pass U+F with shared intermediates for LJ** (and the analytic
   potentials): compute r⁻⁶ once, reuse for U, F, and the shift terms; drop
   redundant `powi`/`sqrt` (dreid-kernel shows ~40% on the combined pass [12]).
4. **Kill per-step allocations** in the hot loop (neighbor buffers, forces —
   audit-flagged).

Expected: **1.3–1.5×** on existing benchmarks.

### P1 — Compiler flags, measured (½ day)

A/B with the criterion benches in `crates/modena-bench`: `target-cpu=native`,
`lto=fat`, `codegen-units=1`. Expect **~10–15%** typical [15], but
`native`+fat-LTO has documented pathological regressions (rust#146497, #48371)
[16] — keep the flags only where the per-kernel numbers win. This also tells
us whether LLVM already auto-vectorizes the pair loop before deciding on
hand-SIMD.

### P2 — Cache + branch elimination (1 day)

1. **Sort atoms by spatial locality** at list build (advice #2).
2. **Per-type-pair sorted iteration**: iterate pairs grouped by `(ti, tj)` so
   the inner loop is branch-free and each table lookup amortizes (MODENA's
   `energy_force` pays the shift/type branches per pair — exactly LAMMPS OPT's
   5–20% [9]; a type-sorted inner kernel can do better).
3. Re-check the LoA: the only fetched layout study (Watanabe) found
   **AoS-with-padding beat SoA on AVX2** (minor on AVX-512) [11] — benchmark
   before assuming SoA stays optimal.

### P3 — The paper's performance centerpiece: Apple-GPU, measured (1–2 days)

1. Benchmark `GpuPairForceModel` vs CPU at N = 1k/4k/10k/50k on the M3 Pro.
   The O(N²) all-pairs kernel wins at moderate N (2.13× at 11k, measured in
   `bench/scaling/SCALING.md`); publish the crossover, then consider a GPU cell
   list **only if the measurement demands it**.
2. State it plainly: **LAMMPS has no Apple-GPU path** (GPU = CUDA/HIP/OpenCL
   [1]; OpenCL hangs on Apple silicon [2]; Kokkos no Apple backend [3]). This
   is MODENA's only head-to-head where LAMMPS cannot compete — and it is
   currently unmeasured.

### Not worth it before the paper (explicitly)

- **Hand-written NEON f64x2 LJ kernel**: the best published SIMD gain on LJ is
  +42% [11], i.e. ~1.1–1.25× on a full step; only after P1 shows LLVM failing.
- **f32 force accumulation on CPU**: NEON f32x4 vs f64x2 ceiling is 2×, eaten
  by memory-boundedness; and single-precision diverges trajectories fastest at
  small N — **a direct threat to the conservation contract**, our actual
  contribution. (On GPU, mixed mode f32-pair/f64-accumulate is standard — LAMMPS
  does it too.)
- **`std::simd` portable SIMD**: still nightly-only [17]; use `std::arch`
  (stable NEON on aarch64 [18]) or `wide`/`fearless_simd` (stable, maintained)
  [19][20] if hand-SIMD is ever justified.
- **GROMACS cluster-pair rewrite**: state of the art, but a rewrite — after.
- **Beating `lmp_serial` at small N**: fixed per-step overheads dominate; nobody
  publishes that number.

## Bottom line

1. Do P0+P1 now (≈1–2 days) → ~1.5–2× behind LAMMPS single-core, from ~4.5×.
2. Do P3 before writing the performance section → the Apple-GPU measurement is
   the paper's only "we win" and LAMMPS forfeits it by absence.
3. P2/structural rewrites are post-paper.

## Sources

[1] LAMMPS GPU package — https://docs.lammps.org/stable/Speed_gpu.html
[2] lammps#3531 (OpenCL hang on Apple Silicon) — https://github.com/lammps/lammps/issues/3531
[3] kokkos#8461 (no Apple GPU backend) — https://github.com/kokkos/kokkos/issues/8461
[4] JAX-MD — https://arxiv.org/abs/1912.04232 ; perf report — https://github.com/jax-md/jax-md/issues/377
[5] TorchSim (IOP 2025) — https://doi.org/10.1088/3050-287X/ae1799
[6] Plimpton 1995 — https://doi.org/10.1006/jcph.1995.1039
[7] Thompson et al. 2022 — https://doi.org/10.1016/j.cpc.2021.108171
[8] LAMMPS "Measuring performance" — https://docs.lammps.org/latest/Speed_bench.html
[9] LAMMPS OPT — https://docs.lammps.org/latest/Speed_opt.html
[10] LAMMPS INTEL — https://docs.lammps.org/stable/Speed_intel.html
[11] Watanabe & Nakagawa, LJ SIMD — https://arxiv.org/abs/1806.05713
[12] dreid-kernel (Rust) — https://github.com/caltechmsc/dreid-kernel/blob/51922aab51f5a064c826623279d3e0aff4800ad6/BENCHMARKS.md
[13] MD-Bench — https://arxiv.org/abs/2302.14660 ; repo https://github.com/RRZE-HPC/MD-Bench
[14] LAMMPS INTEL test-suite — https://raw.githubusercontent.com/lammps/lammps/develop/src/INTEL/TEST/README
[15] rust#93321 (LTO/codegen-units gains) — https://github.com/rust-lang/rust/issues/93321
[16] rust#146497 / #48371 (flag regressions) — https://github.com/rust-lang/rust/issues/146497 , https://github.com/rust-lang/rust/issues/48371
[17] std::simd nightly-only — https://doc.rust-lang.org/nightly/std/simd/index.html
[18] std::arch aarch64 NEON stable — https://doc.rust-lang.org/stable/core/arch/aarch64/index.html
[19] wide crate — https://docs.rs/wide/latest/wide/
[20] fearless_simd 1.0 — https://linebender.org/blog/fearless-simd-1-0/
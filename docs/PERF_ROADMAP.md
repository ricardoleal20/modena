# Performance roadmap — approaching LAMMPS (and where we can win)

> **IMPORTANT — superseded by measurement (2026-09-26, committee acta).**
> A four-expositor research committee plus an audit on the actual repo made
> these quantified corrections. See `docs/COMMITTEE_ACTA.md` for the full
> consensus record (4/4 vote), the doubt-resolution registry, and the final
> gated pipeline. Executive corrections:
> - **Baseline**: the honest single-core gap is **~2.1× flat** (2.06× at 864,
>   2.10× at 2048, pinned 1 thread), inflated to 2.5–4.6× by default-thread
>   self-harm (rayon 11t is 31% *slower* than 1t at N=864). The "3–4.5×" and
>   committed CSV numbers below are stale/pre-thread-fix. There is **no cache
>   superlinearity** at 1t + skin — the earlier "17→38 ns/pair" was a rebuild
>   artifact of the force-kernel bench (52–68% rebuild).
> - **Charged systems**: MODENA direct-Ewald is **~34× behind LAMMPS** on the
>   matched 216-atom deck (341.99 s vs 9.93 s for 2000 steps; 171 ms/step,
>   ~16.9M sincos/step at ~10 ns each = the entire step). Ewald-b (±k symmetry
>   + phase reuse) leaves ~8.6×; only per-axis trig tables (~7 FMA/sincos,
>   +0.5–1 d) or SPME (O(N log N), post-paper) reach ~2×.
> - **Compiler flags**: measured **0%** on aarch64 (release already
>   `lto="thin"` + `codegen-units=1`); P1 is an appendix row, arm64-scoped.
> - **Bit-determinism**: byte-identical at 1 thread; **not** at default threads.
> - Pipeline order (gates G0–G5): 0 rig → 1 thread plumbing (serial-below-N,
>   `with_min_len`, scratch reuse; the biggest default-pool lever, 1.3–2.1×) →
>   2 kernel PR (B monomorphize — bit-identical, unblocks everything; C fused,
>   gated; K rebuild half-stencil for dims ≤ 3) → 3 GPU-e Kahan + persistent
>   buffers & GPU-d ladder (G3 decides the only optional "win" claim) →
>   4 Ewald-b exact hardening → 5 freeze.

Research date: 2026-09-26. Method: two-agent deep research (bibliography of
LAMMPS optimizations + adversarial verification), synthesized below; the
committee (4 expositors, 2 cross-review rounds, audit, 4/4 convergence vote)
refined it. Baseline (measured this page's date, corrected): MODENA ~2.1×
behind `lmp_serial` single-core at 1t on LJ (864/2048), inflated to 2.5–4.6×
by default threading.

## Executive summary (honest envelope, committee-corrected)

- **Single-core CPU: closing focused at 1t from ~2.1× to ~1.4–1.7× is the
  defensible pre-paper target** — via K (rebuild half-stencil, rebuild share
  32%→<15%) × B (monomorphize, 5–20% pair time, bit-identical) × C (fused,
  gated). The thread-plumbing lever (A) fixes the *default* pool (1.3–2.1×),
  not the 1t baseline. **Single-core parity is not a first-paper target**;
  chasing it invites a reviewer benchmark rejection. (Flags contribute ≈ 0%
  on arm64.)
- **The only optional "we win" claim is the Apple GPU**, gated on a measured
  win (G3): wording = "the only supported/tested/validated Apple-GPU path",
  LAMMPS has no Metal backend (GPU = CUDA/HIP/OpenCL [1]); OpenCL degraded
  `neigh no` on Apple silicon (lammps#3531 [2]); Kokkos none (kokkos#8461
  [3]). The current all-pairs GPU kernel likely loses to `lmp_serial`
  everywhere (estimates: 0.87–1.43 vs ~0.25–0.35 µs/(atom·step)); the
  cell-list kernel's ~25–100× is **unmeasured** (Metal gather-throughput spike
  is G3).
- **Charged systems**: direct Ewald is measured ~34× behind LAMMPS at N=216;
  correctness-only for the paper unless Ewald-b + trig tables land (→ ~2×),
  else SPME post-paper.
- **Framing**: "correct-by-construction Rust engine, measured LJ gap ~2.1×
  closed to ~1.4–1.7× (1t), charged path characterized, bit-deterministic at
  1t, and the only supported Apple-GPU acceleration path." Matches the 2025–26
  literature reward structure (JAX-MD published slower — differentiability;
  TorchSim conceded classical speed — MLIP ergonomics [4][5]).

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

### P0 — Measure first, then delete (committee-corrected ordering)

> Note: the original P0 items were **re-scoped by measurement**. The Ewald
> k-vector regeneration is <1% of Ewald time (and Ewald is absent from every
> LJ deck) → hygiene inside Ewald-b; the "second neighbor build" is already a
> no-op in production (skin) → docs note. The real first moves are A (thread
> plumbing — the default-pool lever) and K (rebuild half-stencil), because the
> rig proved orchestration/rebuild, not FP codegen, dominate.

1. **Thread/orchestration (Stage 1, G1)**: serial path below a size threshold;
   `with_min_len` so ~17k pairs don't fan into 11+ subtasks; reuse per-thread
   scratch instead of per-fold `vec![0.0; 3n]`. Measured ceiling: 1.31× (864)
   to 2.33× (2048) on the default pool; 1t trajectory stays byte-identical.
2. **Rebuild half-stencil (Stage 2, K, G2)**: at dims ≤ 3 per axis the ±1
   stencil covers the whole box → ~11× overscan of distance checks. dims-aware
   rebuild (plain contiguous j>i scan when dims ≤ 3, or finer cells with a
   computed stencil radius). Rebuild share 32%→<15% of 1t step.
3. **Monomorphize the potential (Stage 2, B, G2)**: kill the per-pair vtable
   call and four table indirections; bit-identical, unblocks inlining.
4. **Fused single-pass U+F with shared intermediates (Stage 2, C, gated)**:
   r⁻⁶ once for U, F, virial (dreid-kernel ~40% on the combined pass [12]).

Expected (1t path): **~1.4–1.7×** cumulative (K × B × C); default-pool path
additionally **1.3–2.1×** from A.

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
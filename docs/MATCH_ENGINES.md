# Matching the representative MD / computational-biology engines

Five-agent deep-research synthesis (2026-09-26): what MODENA can adopt to go
from ~2× behind `lmp_serial` single-core to *matching* — and being accepted
alongside — the engines that define the field (GROMACS, OpenMM, AMBER, NAMD,
HOOMD-blue). Sources and per-engine detail are in the agent reports; this is
the consolidated, prioritized, code-mapped recommendation.

## What each reference engine actually does differently

| Engine | The core speed idea | Evidence (cited in the agent reports) |
|---|---|---|
| **GROMACS** | Cluster-pair (MxN) kernels: i-clusters × j-clusters, SIMD-width blocks, no per-pair gathers; `r^-6 = (1/r²)³` with `rsqrt`+Newton (no sqrt, no pow); kernel **modules** code-generated per (electrostatics × VDW × flags); 14-instruction-set SIMD layer; dual list + cheap rolling prune | Päll et al. JCP 2020 / arXiv:2006.09167; source `nbnxm/`; MD-Bench arXiv:2302.14660 |
| **OpenMM** | Two-tier **ReferencePlatform / CpuPlatform** (correct oracle + fast SIMD path — MODENA's split already mirrors it); f32 `fvec4` kernels, **per-thread force buffers**, voxel list with skin = 0.25·cutoff, deterministic-forces mode | OpenMM devguide ch.4–5; CpuKernels.cpp; OpenMM 7 paper |
| **AMBER (GPU)** | W×W **warp-tile** pair kernels (no atomics), **SPFP** fixed-point accumulation ("energy conservation equivalent to double"), GPU-resident streaming | Götz et al. 2012; Amber24 manual §22.5–6 |
| **NAMD** | Charm++ compute-objects + **GPU-offload then GPU-resident**; PME charge-spread on GPU, FFT on CPU; tile non-bonded kernels in CUDA and AVX-512 | Phillips 2005/2020; NAMD 3.0 release |
| **HOOMD-blue** | GPU-native: cell-list on GPU → particle sort → per-atom j-ranges; **per-block register accumulation, no atomics**; runtime autotuner over block sizes | Glaser 2015 (12.5×/node **on 2014 discrete GPUs — not transferable to an iGPU**); HOOMD docs |
| **Rust/new-gen + literature** | dreid-kernel shows combined U+F = 1.19 ns/eval (shared r⁻⁶); MD-Bench: half-list vs full-list crossover is SIMD-width dependent; Watanabe: hand-SIMD LJ = +34–46%, AVX2 ≈ AVX-512; **splines preserve dU/dr=−F only if force = derivative of the same interpolant** (never LAMMPS-style separate U/F splines under a strict contract) | dreid-kernel BENCHMARKS.md; arXiv:1806.05713; arXiv:2302.14660; GROMACS manual (tabulated) |

## The consolidated matching roadmap (code-mapped)

### Tier A — close the LJ single-core gap (2.0× → ~1.1–1.3× of `lmp_serial`)

1. **Cluster-pair (MxN) inner loop** — `force_model.rs` + `neighbor.rs`. The
   flat half-list of `(i,j)` tuples (16 B/pair, two gathers/pair) becomes
   i-clusters (4) × j-cluster j-lists; 128-bit NEON loads one cluster
   contiguously; Newton's 3rd law per i-cluster. This is GROMACS's core
   invention and the single change most likely to retire the remaining gap.
   Effort 3–5 d. Keep the scalar kernel as the oracle (GROMACS keeps
   `kernels_reference/` for the same reason).
2. **The LJ fast path without sqrt/pow** — `potential/lennard_jones.rs`:
   `r⁻⁶ = (1/r²)³` with `vrsqrteq_f64` + 1–2 Newton steps (stable,
   `std::arch::aarch64`), keeping the FD `dU/dr = −F` gate to choose the
   refinement count (MD-Bench: each NR iteration ≈ up to 11% DP). Effort 1–2 d.
3. **Half-stencil + locality** — `neighbor.rs`: dims≤3 fast path is in
   (`6785396`); extend to per-axis cell-sorted atom order (Morton-lite) when
   the cluster kernel lands. Effort 1 d.
4. **Per-thread force buffers with deterministic reduction** — OpenMM's
   pattern; compose with rayon, add a `deterministic=true` mode (fixed chunk
   order) for the "reproducible down to summation order" story. Effort 1–2 d.

Projection: cluster-pair × rsqrt-LJ × remaining fixing ≈ **1.2–1.5× more**,
i.e. ~1.0–1.3× of `lmp_serial` — matching, per the committee's fix + the
engines' own published technique program (GROMACS's 2015→2020 kernel work was
worth ~2× on identical hardware).

### Tier B — protein-scale (explicit water + PME), the "computational biology" claim

5. **SPME** (B-spline spread → FFT → ik-solve → interpolate), `ewald.rs`,
   with `rustfft` (CPU; NEON path) — direct Ewald is O(N²) with a measured
   10.8× residual even after Ewald-b; the reciprocal FFT path is the only way
   a 20k-atom explicit-water benchmark exists at all. Preconditions recorded:
   strain-form reciprocal virial on the grid (kills the documented +0.00632
   pressure defect), matched-κ/order/accuracy protocol vs LAMMPS `pppm`.
   Effort 1.5–3 wks; post-paper unless the charged claim goes in.
6. **Exclusions + SETTLE** for water; force-field input path (amber-prm/CHARMM
   subset); HECBioSim 20K Crambin as the community-benchmarked protein target
   (1-2/1-3/1-4 exclusions are the actual missing piece for any protein).

### Tier C — the GPU "win" (only if G3 says so)

7. **GPU-e** (Kahan compensation + persistent buffers, fixed overhead ~2 ms →
   ~0.3 ms) then **GPU-d** ladder (modena CPU/GPU vs lmp_serial, fence-synced,
   drift-classified). **GPU cell-list gather** (CPU bins → contiguous
   per-atom j-ranges → one thread per i, no atomics — the HOOMD/AMBER design)
   only if the ladder demands it. WGSL has **no f64** on Metal: Kahan in f32 +
   f64 host reduction is the only route; the "only supported/tested Apple-GPU
   path" claim needs the `neigh no` caveat (lammps#3531) + same-table
   precision disclosure.

### Benchmark methodology (what makes every claim citable)

Hoefler & Belli (SC'15) protocol: median (+spread) over ≥5 runs after warm-up
discard, runs ≥30 s, report the design; pinned threads for the headline
(`RAYON_NUM_THREADS=1` both engines); full stamp (SHA, rustc, lmp version,
flags); `µs/(atom·step)` — LAMMPS's own metric. GPU: fence-synced timing,
min AND median, precision disclosed in the same table. `bench/bench.sh`
(`6785396`) encodes this; the paper's numbers must come from an idle-machine
regeneration.

## What the research explicitly says NOT to do

- Not spline tables for LJ (GROMACS's own manual: lookup is slower than analytic
  LJ); tables only for tabGAP-style/user potentials, with the
  exact-derivative spline pattern.
- Not f32-on-CPU for MODENA's headline path (contract math; 2× SIMD-width is
  only 11–38% on these gather-bound kernels anyway — MD-Bench measured SP
  vs DP).
- Not FFTW-FFI over rustfft (NEON path mature; <30% at power-of-2 grids).
- Not chasing GROMACS/AMBER absolute GPU throughput on an iGPU (M3 shares
  ~150 GB/s with the CPU; a 1.5–2 TB/s discrete-GPU number is a different
  regime — cite and scope honestly).
- Not beating `lmp_serial` at N ≤ 256 as a claim.

## Realistic targets (all honest, per the research)

- LJ single-core parity-adjacent (1.0–1.3×) via Tier A.
- Small protein with PME: "not embarrassed" = within ~5–10× of GROMACS
  reference class (0.3–0.5 µs/(atom·step) all-cores on the M3) — conditional
  on SPME existing; the direct-Ewald 10.8× row means no protein number at
  all without it.
- GPU: match-class OpenMM-OpenCL on M4/M5 Max (DHFR-PME ~600-750 ns/day),
  never AMBER-class (1,600+).

Primary sources: GROMACS JCP 2020 (arXiv:2006.09167) + `nbnxm` source; MD-Bench
(arXiv:2302.14660, arXiv:2207.13094); Watanabe & Nakagawa (arXiv:1806.05713);
OpenMM devguide + CpuPlatform/CpuKernels + 8.0 notes; AMBER (Götz 2012,
Amber24 manual, GPUPerformance.php); NAMD (Phillips 2005/2020, NAMD 3.0); HOOMD
(Glaser 2015) + docs; dreid-kernel BENCHMARKS.md; HECBioSim suite; Hoefler &
Belli SC'15; LAMMPS Speed_* docs.
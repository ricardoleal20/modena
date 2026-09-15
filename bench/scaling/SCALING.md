# MODENA M9 Parallel-Performance Scaling Benchmarks

Measured on 2026-09-15. All numbers are wall-clock timings of actual runs —
none are extrapolated.

## Hardware

| Component | Details |
|-----------|---------|
| CPU       | Apple M3 Pro (11 cores: 5 performance + 6 efficiency) |
| RAM       | 18 GB unified memory |
| GPU       | Apple M3 Pro built-in (14-core, Metal 4) |
| MPI       | Open MPI 5.0.10 (`mpiexec --oversubscribe --allow-run-as-root`) |
| Toolchain | Rust 1.97.1 (c980f4866, 2026-06-30) |
| OS        | macOS (Darwin 25.6.0) |

## Methodology

All benchmarks use an FCC Lennard-Jones system in reduced units (sigma = eps =
1, cutoff = 2.5, rho ~ 0.8442, lattice constant a = 1.6796). Positions are
static throughout the timed loop — this measures steady-state force-kernel
cost (ghost exchange + pair computation), not integration overhead.

- **MPI strong scaling**: fixed problem size (~10 976 atoms), 100 force
  evaluations, rank counts 1/2/4/8. Grid chosen by `RankGrid::balanced`
  (subdomain edge >= cutoff). Timed region bracketed by `world.barrier()` on
  all ranks. One untimed warm-up step precedes the timed loop. Speedup =
  T(1) / T(N).
- **MPI weak scaling**: atoms scaled with ranks (~1 372 atoms/rank), 100
  force evaluations. Parallel efficiency = T(1, N0) / T(N, N*N0).
- **GPU vs CPU**: `GpuPairForceModel` (f32 WGSL O(N^2) all-pairs shader) vs
  `PairForceModel` (f64 rayon O(N) cell-list half-list), 100 force
  evaluations per size. Each step zeroes forces + virial before the compute
  call (trait contract is accumulate). One untimed warm-up step precedes the
  timed loop. Speedup = T_cpu / T_gpu.

Binaries built with `cargo build --release` (optimized). Source:
`crates/modena-core/examples/mpi_scaling.rs` and `gpu_scaling.rs`.

Raw CSV data: `bench/scaling/results/`.

## Results

### Strong Scaling (10 976 atoms, 100 steps)

| Ranks | Grid      | Wall (s) | us/step | Speedup |
|-------|-----------|----------|---------|---------|
| 1     | [1,1,1]   | 2.092    | 20 917  | 1.00    |
| 2     | [2,1,1]   | 3.308    | 33 079  | 0.63    |
| 4     | [2,2,1]   | 1.895    | 18 949  | 1.10    |
| 8     | [2,2,2]   | 3.717    | 37 171  | 0.56    |

### Weak Scaling (~1 372 atoms/rank, 100 steps)

| Ranks | Atoms  | Atoms/rank | Wall (s) | us/step | Parallel efficiency |
|-------|--------|------------|----------|---------|---------------------|
| 1     | 1 372  | 1 372      | 0.266    | 2 662   | 1.00                |
| 2     | 2 916  | 1 458      | 0.386    | 3 858   | 0.69                |
| 4     | 5 324  | 1 331      | 0.958    | 9 585   | 0.28                |
| 8     | 10 976 | 1 372      | 3.660    | 36 604  | 0.07                |

### GPU vs CPU Speedup (100 steps)

| Atoms  | CPU wall (s) | CPU us/step | GPU wall (s) | GPU us/step | Speedup |
|--------|-------------|-------------|-------------|-------------|---------|
| 864    | 0.132       | 1 322       | 0.219       | 2 189       | 0.60x   |
| 10 976 | 2.045       | 20 450      | 0.959       | 9 588       | 2.13x   |
| 48 668 | 12.609      | 126 087     | 6.963       | 69 626      | 1.81x   |

## Interpretation

The MPI strong-scaling results are poor at this problem size (~11K atoms):
only the 4-rank configuration achieves a speedup (1.10x), while 2 and 8 ranks
are slower than serial. The root cause is the blocking point-to-point ghost
exchange: at ~11K atoms each rank owns only ~1 400 atoms, and the 26-neighbor
halo exchange (with blocking send/receive ordering) dominates the force
computation. The 2-rank [2,1,1] grid is worst because each rank sends its
entire owned set to a single neighbor, and the blocking exchange serialises
communication. Weak scaling confirms this: efficiency drops to 7% at 8 ranks —
communication cost grows superlinearly with rank count because each rank
exchanges with up to 26 neighbors while the per-rank compute shrinks. These
results are honest for a blocking-exchange prototype on a single shared-memory
node; production MD engines achieve good scaling with non-blocking
communication (compute/communication overlap) and much larger problem sizes
(10^5-10^6 atoms) where the compute-to-communication ratio is favourable.
The GPU results tell a different story: the f32 WGSL shader is 2.1x faster
than the f64 rayon CPU path at 11K atoms despite using an O(N^2) all-pairs
algorithm (vs the CPU's O(N) cell-list). At 864 atoms the GPU is slower (0.6x)
due to buffer allocation and kernel launch overhead, and at 49K atoms the
speedup drops to 1.8x as the O(N^2) scaling erodes the throughput advantage.
The f32 vs f64 precision difference (~1e-5 relative) is the standard GPU MD
trade-off and does not affect the speedup comparison.

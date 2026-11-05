# Native CPU pair scratch reuse (B02)

`PairForceModel` owns a reusable dense parallel workspace. The serial path still
applies below 100,000 half pairs and in a one-thread Rayon pool. Its arithmetic
and `compute_deterministic` pair order are unchanged. LJ specialization and the
generic `PairPotential` path share the workspace policy. B09's immutable,
prebuilt bounded snapshot list reaches the force pass directly, even with zero
Verlet skin.

The model attempts a nonblocking `try_lock` lease for a parallel pass. On
contention (or poisoning), it computes serially. This permits shared models and
nested same-pool evaluations without waiting on a lock held by fork/join work.
The contention path can be slower; sharing one model among concurrent
simulations is supported for correctness, not promised to scale in throughput.

Each useful partition owns a `3N` f64 force vector plus energy and scalar
virial. Every pass resets all used force/scalar elements. Partitions cover
contiguous half-pair chunks, with at most the current pool's worker count.
Local scalar sums are written back once per partition. Disjoint force output
blocks reduce the dense partitions in parallel; this remains a dense O(TN)
workspace, not a sparse/touched-atom accumulator.
Parallel pool shrinkage drops excess partitions. Remaining vectors retain
their largest atom-count capacity; serial/empty calls leave the idle cache
untouched. Drop releases the whole workspace. Thus retained memory is O(TN),
not sparse or bounded independently of system size. Stable warmed passes
allocate no N-sized force buffers; initial use, N growth and partition growth
can allocate. Allocator calls alone do not establish a speed improvement.

## Rust migration

External `PairForceModel { potential }` literals no longer compile because the
workspace is private. Replace them with `PairForceModel::from_boxed(potential)`
or `PairForceModel::new(concrete_potential)`. Both constructor signatures and
the public `potential` field are unchanged. `ForceModel` remains `Send + Sync`.
No TOML or checkpoint migration is needed.

## Reproduction and numerical policy

`pair_scratch_profile ATOMS THREADS PASSES` is an isolated example with an
instrumented allocator; its timings are not used as primary performance
evidence. It records actual pair counts, force-sized allocation/reallocation
counts and bytes, live allocated bytes/high water, released model bytes and
fixed static numerical gates. The allocator only exists in the test/example.
Large static E/W gates are `atol=1e-8, rtol=1e-12`; every force component uses
`atol=1e-10, rtol=1e-12`. One-thread E/F/W uses `to_bits` equality, including
the sign of zero; CLI checkpoint/thermo comparisons use actual file hashes.

Build ordinary baseline and candidate release CLIs with the same compiler,
flags and checkout, and preserve both before rebuilding either:

```sh
rtk cargo build --locked --offline --release -p modena-cli -p modena-core --bin modena --example pair_scratch_profile
rtk proxy "$PY" bench/pair_scratch.py --before "$BASELINE" --after "$CANDIDATE" \
  --sizes 864 2048 5000 10000 --threads 1 4 8 --steps 1000 \
  --repetitions 5 --output "$FRESH_EVIDENCE"
```

The stdlib-only POSIX Python >=3.11 driver invokes `git`, `rustc` and macOS
`sysctl` directly for provenance; it does not require RTK. The RTK prefixes in
these examples are the agent shell convention. It alternates AB/BA pairs after separate
discarded fresh-process warmups. It retains every input, stdout/stderr, full
checkpoint, thermo log and wall/user/system CPU/RSS/load row, plus source and
binary hashes. It rejects existing output directories, failed subprocesses,
incomplete steps, nonfinite states, wrong physical signatures and one-thread
checkpoint/thermo differences. Multi-thread chaotic trajectories are validated
for complete/finite matched physics, with static numerical and scientific
gates handled separately.

All runs use reduced LJ sigma=epsilon=mass=1, cutoff=2.5, skin=0.3,
FCC a=1.6796, velocity Verlet dt=0.005, T=1 and seed=87287. The sizes use
replicas [6,6,6], [8,8,8], [10,5,25], [10,10,25], respectively. Perturbed
static profiling yields 33,696 / 79,872 / 195,000 / 390,000 half pairs;
N=864/2048 are serial controls even in larger pools.

## Measured acceptance and limits

The retained block variant has a reproducible scoped full-CLI benefit at
N=5000, four threads on this shared Apple M3 Pro (11 logical CPUs). Two separate
1000-step series gave before/after median ratios 1.108× (three pairs) and
1.086× (five pairs): 1.356→1.248 seconds in the complete matrix. All eight paired
ratios favored the candidate. This supports retaining the variant for that
case; it does not establish a general speedup for all systems, pools or hosts.

The initial dense variant used a serial force reduction. It allocated no warm
force buffers but regressed the four large parallel cases by about 20–30%,
with higher median RSS. It was rejected. The final variant combines local
energy/virial scalars and parallel output blocks; their individual performance
contributions were not isolated. Both variants and all negative rows remain.

Ordinary-release full-CLI results, 1000 steps and five AB/BA pairs per case:

| N | threads | before median s | block median s | before/block |
|---:|---:|---:|---:|---:|
| 864 | 1 | 0.626840 | 0.638524 | 0.982× |
| 864 | 4 | 0.619774 | 0.615893 | 1.006× |
| 864 | 8 | 0.611273 | 0.611487 | 1.000× |
| 2048 | 1 | 1.784583 | 1.820879 | 0.980× |
| 2048 | 4 | 1.846600 | 1.854466 | 0.996× |
| 2048 | 8 | 1.774144 | 1.761625 | 1.007× |
| 5000 | 1 | 4.094979 | 4.072741 | 1.005× |
| 5000 | 4 | 1.355683 | 1.248346 | 1.086× |
| 5000 | 8 | 1.060009 | 1.024522 | 1.035× |
| 10000 | 1 | 8.487229 | 8.619399 | 0.985× |
| 10000 | 4 | 2.739052 | 2.695475 | 1.016× |
| 10000 | 8 | 2.032944 | 1.944128 | 1.046× |

Serial controls varied by roughly 2%; no serial performance benefit is claimed.
N5000/8T's separate pilot ratio was 1.033×. Larger N results were noisier:
separate 3000-step/five-pair N10000 repeats gave ratios 1.052× at 4T and 1.027×
at 8T, with regressing individual pairs. These durations/series are kept
separate; their samples are not pooled into the table above.

The normal host pool count is 11. A separate 11T/1000-step series gave 1.069×
at N5000 (all five paired ratios positive), but **0.791× at N10000**, as both
binaries slowed sharply during rising host load. A separate 11T/3000-step
N10000 repeat gave 1.029×, with one regressing pair. Default-pool large-N wall
performance therefore remains uncertain. No uniform default/large-N speedup
is claimed, and the negative series is retained in full.

Child CPU and wall time are distinct. For N5000/4T in the complete matrix,
median child user time fell 4.367→4.100 s while system time rose 0.124→0.143 s.
At N5000/11T, wall time fell while user/system CPU rose 5.591/1.115→5.854/1.585 s.
N10000/11T's long repeat also used more child CPU. Thus a wall improvement does
not imply lower CPU consumption or energy use. Recorded one-minute loads
ranged 3.45–8.63 in the complete block matrix, 6.22–22.07 in the default-count
series, and 12.43–24.51 in the long repeats. Our other build/test jobs were idle;
unrelated host work remained active. These are local shared-host experiments.

Allocation/static profiles confirm zero warmed force-sized allocations at
parallel N5000/10000 and T4/8/11; baseline allocated 8 buffers per pass at 4T,
16 at 8T and 16 at 11T. All static gates passed, with maximum force-component
error 5.67e-15. All ordinary CLI samples passed complete finite matched-physics
checks; all one-thread full checkpoint/thermo comparisons were exact.

Fresh-model retained scratch bytes (including partition metadata) are:

| N | 4T | 8T | 11T |
|---:|---:|---:|---:|
| 5000 | 480,224 | 960,384 | 1,320,504 |
| 10000 | 960,224 | 1,920,384 | 2,640,504 |

This increases idle live memory. At N10000/8T, the isolated Rust-requested
allocated-byte peak fell 14.538→12.618 MB, while full-CLI median OS RSS in the
complete matrix rose 29.098→33.915 MB. RSS was mixed across cases; neither a
general RSS reduction nor bounded memory independent of T/N is claimed.

Baseline source `7c7fa34` and measured production source `0cea614` were clean
builds in the same worktree using rustc 1.97.1/aarch64-apple-darwin and identical
unset Rust flags. At pre-review follow-up `a03dab0`, all 48 measured core/CLI
production source/manifest hashes and retained CLI bytes were identical.
Review fix `3ec0e69` corrects four documentation lines in `force_model.rs`:
47 source/manifest hashes remain identical; that file's hash changes solely
because of the same-line-count comment correction. The ordinary CLI rebuilt
from this clean review-fix commit is **byte-identical** to the immutable
`0cea614` CLI (SHA256
`1b39924b2ab9141748c49a1936df9a57cda8a70d3b879c5008da33996e0113f8`).
The comment now distinguishes fixed chunk/reduction order at a fixed leased
pool from order changes with pool size or serial contention fallback.
Raw untracked local artifacts under `target/b02-evidence/` are unavailable in a
fresh clone; rebuild/rerun to reproduce:

- `baseline/`, `dense/`, `block/`: immutable CLI/profile binaries, compiler,
  clean-source epochs, binary hashes and all allocation/static profile rows.
- `dense-full-cli/`: rejected 144-run series and retrospective box/time guard
  audit. The loaded driver hash/initial source are distinct from the later
  guard committed as `0248372`; no timing rows were replaced.
- `block-pilot/`, `block-full-cli/`, `block-default11/`,
  `block-long-n10000/`: each separate raw series, every input/log/checkpoint,
  wall/user/system/RSS/load row, paired ratios and source/binary provenance.
- `default11-profiles.json`, `production-equivalence.json`: non-power-of-two
  static/allocation checks and historical pre-review test-only equivalence.
- `fix-round-1/`: RTK-absent/signed-zero RED/GREEN, rebuilt CLI equivalence,
  fresh 14-case bit-gated profiles, and retrospective validation of all 380
  retained rows and 48 one-thread byte-exact file pairs. Original loaded
  driver/source hashes and timing files remain unchanged; the new validator
  hash is recorded separately. No performance retiming was required.

The [fresh B02 verification campaign](B02_VERIFICATION.md) passed 357 Rust
tests (one existing ignored doctest), 35 Python tests with actual ASE, five
worker tests, 34 benchmark-tool tests and the 67-observation scientific smoke.
The supported base interpreter executed 26 tests with nine explicit optional
skips. Default/GPU lint and Rust 1.90 workspace/harness checks passed. The fresh
release CLI remained byte-identical to the accepted measured block binary.
GPU compilation/lint provides no GPU runtime evidence; MPI and other-engine
performance are outside this result.

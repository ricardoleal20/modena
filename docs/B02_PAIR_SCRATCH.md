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
`atol=1e-10, rtol=1e-12`. One-thread results require exact equality.

Build ordinary baseline and candidate release CLIs with the same compiler,
flags and checkout, and preserve both before rebuilding either:

```sh
rtk cargo build --locked --offline --release -p modena-cli -p modena-core --bin modena --example pair_scratch_profile
rtk proxy "$PY" bench/pair_scratch.py --before "$BASELINE" --after "$CANDIDATE" \
  --sizes 864 2048 5000 10000 --threads 1 4 8 --steps 1000 \
  --repetitions 5 --output "$FRESH_EVIDENCE"
```

The stdlib-only Python >=3.11 driver alternates AB/BA pairs after separate
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

Measured acceptance and limitations will be recorded after paired ordinary
release measurements; this description alone makes no runtime claim.

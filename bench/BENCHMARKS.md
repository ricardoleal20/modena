# Reproducible benchmark tools

`benchmark.py` requires Python >=3.11 and the prebuilt release CLI. It generates
matched LJ decks from one set of parameters, uses the actual requested steps,
retains input hashes, stdout/stderr, return codes, wall/user/system CPU and RSS,
seed, hardware, compiler and binary hashes in `benchmark.json`. Missing metrics
are null. A failed process stops the campaign with exit 1 and retains evidence.
The renderer uses numeric atoms/steps. Output directories must be fresh.

```sh
cargo build --release --bin modena
python3 bench/benchmark.py --lammps lmp_serial --ladder 864 2048 20000 \
  --threads 1 --steps 1000 --warmup-steps 200 --runs 3 --output target/bench-new
python3 -m unittest discover -s bench/tests -v
```

`bench.sh`/`run_comparison.sh` delegate to this runner. Set PYTHON to a supported
interpreter. Default threads=1 matches serial LAMMPS; additional thread counts
are scaling observations, never a same-resource speed win. Whole-process timing
includes startup, initialization and integration. No trajectory/thermo files are
written; stdout and endpoint LAMMPS thermo are captured. Warmup is a separate
fresh process. MODENA/LAMMPS Gaussian initializations use different RNGs, the
same T and DOF, and the same Hamiltonian, not identical microstates.

Criterion `force_cached_neighbors` includes force clearing, compute and neighbor
displacement checks with cached lists; `force_with_rebuild` additionally forces
a rebuild each iteration. `neighbor_build` isolates list creation/build.
`integrator_cached_neighbors` measures one Velocity-Verlet step with setup and
neighbor construction outside timing; kicks, clearing, compute and displacement
checks are inside. The owned cloned system and neighbor buffers are also
destroyed before the measured closure returns, so buffer teardown/deallocation
is included. No trajectory startup is included. Workload ladders contain
LJ N256/864/2048, neutral ionic SPME N64/216/512 and molecular SPME N81/192/375.
Molecular fixtures use SPC/E geometry/charges with intramolecular exclusions;
they benchmark electrostatics, not a full water simulation or chemically complete
force field. Default Rayon resources are controlled with RAYON_NUM_THREADS.

```sh
RAYON_NUM_THREADS=1 cargo bench -p modena-bench --bench force_kernel
RAYON_NUM_THREADS=1 cargo bench -p modena-bench --bench integrator_step
```

`gpu_ladder.py`/`gpu_ladder.sh` retain logs, process resources, DOF, precision,
steps, seed and integration rows. Their in-process step timings exclude startup
and initial force; outer resource timings include these. One fresh trajectory
warms the same process. CPU uses f64; GPU is f32 host offload with synchronous
readback. The `gather` probe uses rc=.01 and a different Hamiltonian; its timings
cannot be sold as a production LJ simulation. These scopes must not be combined
with the whole-process CLI/LAMMPS benchmark as a strict speedup.

Missing or unreadable GPU binary hashes are null. OS launch errors retain
return code 127, stdout/stderr and `gpu.json`, then the tool exits 1. The optional
`--timer PATH` selects the resource timer (default `/usr/bin/time`); if that path
is absent, the engine runs directly and unavailable CPU/RSS measurements are null.

```sh
cargo build --release -p modena-core --features gpu --example gpu_ladder
python3 bench/gpu_ladder.py --ladder 864 --steps 100 --runs 3 --threads 1 --output target/gpu-new
```

[Measured results and historical caveats](../docs/scientific-hardening/RESULTS.md).
Original reports and raw baseline values are preserved in the evidence archive.

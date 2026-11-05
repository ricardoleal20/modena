# B09: ASE sessions and calculator verification

B09 implements two separate adapter directions:

- The ordinary MODENA real-unit runner can use pinned official ASE EMT as a
  total energy/forces/stress provider. `transport="one_shot"` remains the
  default; `transport="persistent"` explicitly reuses one Python process and
  EMT calculator per provider. Both use the same physical identity, units and
  compatible checkpoint signature. See [provider usage](B01_EXTERNAL_PROVIDER.md).
- Optional `modena.ase.ModenaCalculator` exposes native shifted-force LJ to ASE
  through `modena snapshot`. It supports energy and forces for one species,
  full PBC and positive diagonal orthorhombic cells. Constraints, nonzero
  charge/spin, mixed species and unsupported cells fail. Tensor stress is not
  implemented. Base `import modena` needs neither ASE nor NumPy. See the
  [calculator guide](../python/README.md#modena-as-an-optional-ase-calculator).

These are bounded interface/numerical checks. They do not establish chemical
or material accuracy, converged ensembles, DFT/ML support, independent-engine
parity or GPU hardware speed.

## Source and retained evidence

The driver and process-boundary tests were committed **before** measurement as
`146c097bb79bf5d82b6e581bc6f2cad970d95e54`. The source tree was clean during
both release builds, the transport campaign and scientific smoke. Native and
adapter production source is unchanged from `e547758`; this later revision
adds the benchmark driver/tests and prior progress documentation.

Local raw evidence is retained under the ignored shared target directory:
`target/b09-evidence/task3-146c097/`. It includes preserved executable copies,
build commands/logs, `transport/ase_transport.json`, every input deck, stdout,
stderr, thermo CSV and final checkpoint, and `scientific-smoke/validation.*`.
Generated binaries and large raw series are not committed. The timing artifact
contains all tracked source-file SHA-256 hashes, input hashes, commands,
thread environment, return codes, completed steps, per-run loads and timings.

| Artifact | SHA-256 |
| --- | --- |
| Sorted source-file hash manifest | `75cc4d8d06e7a41d0476c7d9c441055b2c09e2e9a032d8c2592be992f22bc2fc` |
| `bench/ase_transport.py` | `520da3856900a783f249948bc740ec451f62a4d75dd9caf68c77146640e08a18` |
| Release CLI | `b8d812fa4db93b17d2cda6a642fdc2c72fed57b1206c1e3164a5900ae04c3a19` |
| Release scientific example | `47fb3cb131b089cd77212d41a050f2ce4ea498c0b61a2d62386ebda75d2614d4` |

Host: Apple M3 Pro, 11 logical CPUs, macOS 26.6.2 arm64; rustc 1.97.1
(`8bab26f4f`, LLVM 22.1.6), Python 3.14.5, ASE 3.29.0, NumPy 2.5.4. EMT source
SHA-256 was `8ce8bb2f90e8afad8f091b8c3f600ea260534e5c270a6b2c82cae1817cdf3881`,
with `asap_cutoff=false`. The initial load averages were 3.1924/3.1631/3.1582;
recorded one-minute loads ranged from 2.7393 to 3.3247. This was a shared host,
without concurrent builds/tests during samples, not an isolated machine.

## Full-run timing experiment

Measured 2026-11-05 16:33:16–16:35:15 UTC. Each mode used equivalent FCC Cu
configurations with replicas `[1,1,1]` (4 atoms) and `[3,3,3]` (108 atoms),
lattice constant 3.61 Å, mass 63.546 amu, seed 42, initial temperature 10 K,
Velocity-Verlet and dt 0.0001 ps for **40 integration steps**. Both modes used
thermo every step, no trajectory file, and a final checkpoint at step 40.
`OMP_NUM_THREADS=OPENBLAS_NUM_THREADS=RAYON_NUM_THREADS=1` in both modes.

For each atom count, separate 10-step warmups in fresh processes were discarded.
Five measured pairs alternated one-shot/persistent, persistent/one-shot,
one-shot/persistent, persistent/one-shot, one-shot/persistent. All 24 CLI runs
exited zero; all measured runs completed 40 steps. Each had a new CLI process
and provider session. `time.perf_counter()` measured the whole child process,
including Python startup/imports, discovery, setup, integration, output and
cleanup. It excludes the driver's separate initial pin-discovery call.

| Cu atoms | One-shot median (s) | Persistent median (s) | One-shot / persistent |
| --- | ---: | ---: | ---: |
| 4 | 10.296735208 | 0.303316750 | 33.9471 |
| 108 | 10.852868125 | 0.744721583 | 14.5731 |

Measured wall seconds, in repetition order:

| Repetition | Cu4 one-shot | Cu4 persistent | Cu108 one-shot | Cu108 persistent |
| --- | ---: | ---: | ---: | ---: |
| 1 | 10.382357167 | 0.305705375 | 10.894942250 | 0.748168917 |
| 2 | 10.237236833 | 0.298168917 | 10.852868125 | 0.746060209 |
| 3 | 10.306946500 | 0.303316750 | 10.818342291 | 0.744721583 |
| 4 | 10.240277209 | 0.304100375 | 10.829097292 | 0.738641167 |
| 5 | 10.296735208 | 0.299172042 | 10.964610917 | 0.736879792 |

The driver required exact full checkpoint structural/numeric equality and
exact thermo CSV values across transports and repetitions at each step count.
This includes final full-precision energy, forces, virial, positions,
velocities and physical signature, plus all 41 thermo rows, 41 diagnostic
samples and complete runtime metadata/identity. Warmup comparisons use a
separate 10-step baseline. All comparisons passed. The timings are therefore
not explained by skipped steps, stale cached snapshots or changed physics.

These local ratios measure savings from avoiding repeated Python/EMT startup
in short full runs. They are not isolated steady-evaluation measurements,
pristine cold-file-cache measurements or a claim that EMT arithmetic became
faster. OS file caches may be warm. No controller pilot timings are included;
the results cannot be generalized to DFT, ML, other engines or GPU hardware.

## Reproduction

Use a clean checkout of the measured revision, Python 3.11+ with ASE installed,
and a **fresh** output path. Preserve the venv executable symlink. From the
repository root (set `CARGO_TARGET_DIR` if using a shared build directory):

```bash
rtk proxy cargo build --release --bin modena --locked --offline
rtk proxy cargo build --release -p modena-core --example scientific_validation --locked --offline
rtk proxy /absolute/ase-venv/bin/python bench/ase_transport.py \
  --modena /absolute/target/release/modena \
  --python /absolute/ase-venv/bin/python \
  --output /absolute/target/b09-evidence/new-transport-run
rtk proxy /absolute/ase-venv/bin/python bench/validate.py \
  --binary /absolute/target/release/examples/scientific_validation \
  --output /absolute/target/b09-evidence/new-scientific-smoke
```

`--offline` requires dependencies already cached. The driver discovers exact
installed pins, never overwrites an existing output directory, and exits
nonzero while retaining evidence on command, metadata, completed-step or
parity errors. macOS hardware/resource queries require host permission in a
restricted sandbox. See `build-commands.json` and the timing artifact for the
exact local paths and invocation.

## Regression and numerical checks

The frozen release scientific smoke passed **67/67 observations**, retaining
six raw observable series of 13,667 values each and their sampling analysis.
The raw series remain in `validation.stdout`; the finite-series estimates do
not establish ensemble convergence.

After measurement, the affected benchmark-tool suite passed **27/27** tests
(`python -m unittest discover -s bench/tests -v`), including five new driver
tests using actual child processes for success, alternating order, retained
failure evidence, incomplete steps, identity/thermo/checkpoint mismatch and
overwrite refusal. Rust formatting and whitespace checks passed.

Earlier Task 1/2 validation on unchanged functional source covered default/GPU
lint and actual Rust 1.90 locked checks. Combined non-overlapping workspace
coverage was **351 passing tests and one intentionally ignored doctest**:
the original workspace command had 349 passes and one stale MCP fixture
failure, followed by the corrected complete 47-test MCP suite and doctests.
It was not a zero-exit original workspace command. Python checks passed all
33 with actual ASE, and 26 with seven explicit skips without optional ASE.
These heavy checks were not repeated for driver/documentation changes.
The minimal standalone validation-harness lock repair preserved all prior
pins. GPU runtime/hardware and MPI runtime were not exercised. The upstream
ASE/NumPy shape-assignment warning remains visible in direct reference checks.

Task 1/2 focused read-only reviews and fixes preceded this measurement.
Whole-branch review and remote CI/shipping outcomes belong to the eventual PR;
this record does not claim them. No B09/B02 merge is authorized here.

Persistent session failure/drop kills and reaps the direct worker child;
protocol errors and deadlines leave the session terminal. This is not a
process-tree sandbox. Cancellation is observed at safe setup/step boundaries;
an in-flight call is bounded by its configured deadline plus bookkeeping,
without immediate preemption or a hard global wall-time guarantee.

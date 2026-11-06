# B03 packed layout experiment

This experiment compares scalar SoA, four-lane AoSoA and spatially sorted
four-lane cluster masks with the unchanged deterministic interleaved pair
implementation. Production storage, dispatch and defaults remain unchanged.
The production-auto interleaved implementation is a separate threading control;
its multi-thread results are not causal layout comparisons with serial candidates.
No SIMD or parallel packed implementation is evaluated.

## Reproduce

The driver uses Python 3 standard-library subprocess argv directly, without RTK.
RTK prefixes below are the repository agent-shell convention only.

```sh
rtk cargo build --locked --offline --release -p modena-bench --example layout_step
rtk proxy python3 -m unittest discover -s bench/tests -p test_layout_step.py
rtk proxy python3 bench/layout_step.py pilot --out target/b03-evidence/new-pilot --binary target/release/examples/layout_step
rtk proxy python3 bench/layout_step.py campaign --out target/b03-evidence/new-campaign --binary target/release/examples/layout_step --config /absolute/path/config.json
rtk proxy python3 bench/layout_step.py audit --out target/b03-evidence/new-campaign
rtk proxy python3 bench/layout_step.py secondary --out target/b03-evidence/new-no-snapshot --primary target/b03-evidence/new-campaign --binary target/b03-evidence/new-campaign/layout_step
```

For the 2026-11-06 campaign, save this configuration as a JSON file and pass its
absolute path to `--config`. The seed is fixed at 17 in the campaign driver.
The collection entry point currently uses Unix `wait4` and macOS `sysctl` host
metadata; the recorded campaign ran on macOS. The transparent snapshot reader
and numerical audit use the Python standard library.

```json
{
  "dt": 0.000001,
  "steps": {
    "864": {"no-rebuild": 700, "forced": 190},
    "2048": {"no-rebuild": 200, "forced": 40},
    "5000": {"no-rebuild": 90, "forced": 20},
    "10000": {"no-rebuild": 45, "forced": 20}
  }
}
```

The output directory must be new; existing raw subprocess outputs and snapshots
are never reused. The script copies the ordinary release executable and its own
loaded source, and records SHA-256 hashes, compiler, host CPU, source hashes and
Git HEAD. Launching from a dirty checkout requires preserving untracked source,
Git status and diff separately. A later commit does not change the recorded
launch source epoch. No source or binary changes are allowed during timing.

The frozen configuration declares timestep and steps per N and policy, selected
from an untimed exploratory pilot before the campaign. Current matrix: N864,
2048, 5000 and 10000; requested pools 1, 4, 8 and 0 (actual default reported);
explicit no-rebuild and forced-rebuild policies. Each scalar candidate receives
its own serial baseline: a separate A/B warmup, then five AB/BA alternating
pairs. Every primary timing-positive case receives another warmup and a separate
five-pair reproduction with reversed initial ordering. Production-auto has a
separate warmup and five samples per case and is excluded from scalar acceptance.

All elapsed values are ordinary release measurements. The internally fenced
full velocity-Verlet step includes kick/drift, validation, packing, force
clearing, neighbor maintenance, pair evaluation, scatter and final kick. Cluster
sorting and mask reconstruction are included whenever required. Setup and pool
creation have separate timers. Snapshot serialization, final active counting and
JSON output are outside the full-step interval. External child wall, CPU, maximum
RSS and host load are recorded separately and include process/setup/snapshot
costs. Initial snapshot Value allocations are created during setup and held
through the timed steps. Consequently process RSS and cache state include this
shared instrumentation; absolute RSS is not the footprint of production runs
without snapshots. The current example writes snapshot JSON directly to a File without a
buffered writer. Source inspection and process wall/CPU telemetry suggest
snapshot I/O contributes strongly to collection overhead; no syscall trace was
collected, so this is an inference rather than a measured syscall count. It is
not a kernel speedup or a reason to mix binary epochs. RSS is actual wait4 child high water (Darwin bytes; Linux KiB converted
to bytes), not a subtraction of cumulative maxima. Host load is context, not a
claim that unrelated user processes were controlled.

## Physical and identity gates

The workload is deterministic truncated FCC with a=1.6796, seed 17 SplitMix64,
position jitter 0.002, velocity scale 0.1 and COM velocity removal. Non-4*cube N
truncates the final FCC cells; this is a synthetic workload with vacancies, not
an equilibrated ensemble. All species and atomic numbers are zero, masses one,
charges zero. Orthorhombic periodic cells use reduced shifted-force LJ with
sigma=epsilon=1, cutoff 2.5 and skin 0.3. Unsupported inputs reject explicitly.

The driver checks the actual stdout mode, N, requested/actual threads, timestep,
seed, policy, requested/completed steps, execution class, finite E/F/W, timing
and rebuild/counter fields. Each snapshot includes initial input, force-ready
input and final state. Actual shape, logical IDs, cell/PBC, masses/species/charges/
atomic numbers, reconstructed RNG/FCC coordinates and velocities, forces,
energies, pair identities/scales and completed steps are validated. Two children
sharing the same invalid state cannot pass by parity alone. FNV checksums are
reproducibility metadata, not an integrity or physics proof.

Schema1 `snapshot.row` is exactly stdout **minus `pool_setup_seconds`**. All
shared fields must match; any other missing, extra or mismatched field fails.
The stdout-only pool time must independently be finite and nonnegative. The
snapshot row and stdout are deliberately not treated as interchangeable.

After all timing, the frozen driver independently enumerates linked-cell half
pairs from actual coordinates. It recomputes minimum images and analytical
shifted-force LJ forces, potential, virial, kinetic energy, temperature and
pressure without production helpers. Initial force-ready and forced-final
neighbor sets must equal the independently enumerated skin set. Reused neighbor
lists must cover every final active pair; their original skin set need not equal
the final-coordinate skin set. Oracle arrays and pair sets are cached by an SHA256 key over exact actual
coordinates and validated physical context (cell/PBC, N, masses, species,
charges, atomic numbers and potential). Every force component and E/W value in
every saved force-ready/final state is checked directly against that oracle;
no tolerance is accumulated through a canonical snapshot. A separate trajectory
comparison covers every saved field and unknown mismatches. Both large caches
use explicit four-entry LRUs. Half-pair identities and multiplicity remain
exact. Fixed gates are absolute 1e-8 plus relative 1e-12 for energy/virial and
absolute 1e-10 plus relative 1e-12 for each force component. An audit failure
invalidates the affected measurement, irrespective of a favorable timing.

Task1 Rust tests additionally verify exact mask/permutation mappings, padding,
scaled/excluded half pairs, cutoff and periodic boundaries, current-coordinate
packing, additive forces, cache invalidation, resize/reorder/box changes and B06
convergence. Those correctness tests are not primary performance data.

## Raw evidence and decisions

`rows.jsonl` retains every attempted invocation, including warmups, failed and
unfavorable rows; `index.jsonl` adds phase/pair/candidate/role metadata. The
original stdout and stderr bytes are retained separately with SHA-256. Large
snapshot JSON is compressed after child timing with gzip level 1, mtime 0. Before
removing the duplicate uncompressed file, the complete decompressed bytes are
verified equal. Both original JSON and gzip SHA-256 are recorded. No state or
original byte is discarded. `read_json` transparently reads JSON or gzip; the
audit verifies both hashes. These gzip files can also be inspected using ordinary
`gzip -dc`. They are local ignored evidence, not committed generated artifacts.

`summary-preaudit.json` is only a timing screen. `audit.json` and
`summary-audited.json` include independent correctness invalidation. Improvement
is 1 minus candidate median full-step seconds divided by matched-baseline median
seconds; all five individual paired ratios and their median are also preserved.
The declared promotion bar is at least 5% median improvement, favorable direction
in at least 4/5 pairs, repeated in a separate five-pair campaign, with all
numerical/invariant gates and an explicit memory/regression assessment passing.
Retained packed capacity is partial allocation telemetry, distinct from process
RSS; it excludes core neighbor/workspace and transient cluster reconstruction.

A rejected campaign keeps production unchanged and rejects only the scalar
variants and declared workloads measured here. A positive case requires a
controller-approved concrete promotion spec before core integration, with
explicit experimental selection and unchanged defaults. No engine or chemistry
parity claim follows from this reduced LJ experiment.

## Primary timing screen, 2026-11-06

Launch HEAD was `f13aa727a507f543717af3872202324859adcb78` with the new
driver and tests untracked. They were frozen by full source copies and hashes;
a later commit preserves those same bytes and does not relabel the launch as
clean. Host: Apple M3 Pro, macOS, actual pools 1/4/8/default 11. Release compiler
was rustc 1.97.1, aarch64-apple-darwin, LLVM 22.1.6. Binary SHA256:
`5001bb6b8bdfddc60e74e9fb9bee391e177b1a9992851f47e9ae1a3c8f53c2dc`. Primary loaded driver SHA256:
`5764f6b110841d181f6ce90f670ecca020441aeae5e9d334451fe08df965ef82`.

All 1,344 invocations completed; all 96 timing comparisons failed the 5% bar.
There were no positive cases requiring primary reproduction. The best median
time saving was -1.5416% (SoA, N2048, pool 4, forced), with 0/5 favorable pairs.
The later direct-oracle numerical certification passed all 1,344 invocations.
The separately frozen no-snapshot instrumentation check also rejected every case; its complete results follow below.

The table retains all cases. A is the matched interleaved-serial baseline; B
is the named scalar candidate. Saved%=100*(1-median(B)/median(A)). Setup and
RSS values are medians of each pair side; RSS is each child process high water
including snapshots. Full raw rows retain CPU/load/counters and failed/warmup
records. Pool11d means requested0, actual11.

| N | Pool | Policy | B | A VV ms | B VV ms | Saved % | Favorable | Setup A/B ms | RSS A/B MiB | Packed B KiB |
|---:|---:|:---|:---|---:|---:|---:|---:|---:|---:|---:|
| 864 | 1 | no-rebuild | soa | 255.059 | 370.779 | -45.37 | 0/5 | 2.22/2.35 | 17.94/17.98 | 40.5 |
| 864 | 1 | no-rebuild | aosoa4 | 254.374 | 374.995 | -47.42 | 0/5 | 2.17/2.36 | 17.94/17.97 | 40.5 |
| 864 | 1 | no-rebuild | cluster4 | 255.010 | 436.621 | -71.22 | 0/5 | 2.17/4.53 | 17.95/21.97 | 1994.1 |
| 864 | 1 | forced | soa | 242.291 | 272.797 | -12.59 | 0/5 | 2.13/2.30 | 17.89/17.98 | 40.5 |
| 864 | 1 | forced | aosoa4 | 249.360 | 286.470 | -14.88 | 0/5 | 2.25/2.40 | 17.94/17.98 | 40.5 |
| 864 | 1 | forced | cluster4 | 248.287 | 637.920 | -156.93 | 0/5 | 2.24/4.68 | 17.91/22.05 | 1991.1 |
| 864 | 4 | no-rebuild | soa | 268.096 | 388.780 | -45.02 | 0/5 | 2.30/2.48 | 18.00/18.06 | 40.5 |
| 864 | 4 | no-rebuild | aosoa4 | 265.184 | 388.227 | -46.40 | 0/5 | 2.39/2.49 | 17.98/18.05 | 40.5 |
| 864 | 4 | no-rebuild | cluster4 | 254.720 | 433.506 | -70.19 | 0/5 | 2.21/4.49 | 18.00/22.00 | 1994.1 |
| 864 | 4 | forced | soa | 244.057 | 277.117 | -13.55 | 0/5 | 2.21/2.34 | 18.03/18.03 | 40.5 |
| 864 | 4 | forced | aosoa4 | 246.030 | 280.854 | -14.15 | 0/5 | 2.24/2.43 | 17.98/18.02 | 40.5 |
| 864 | 4 | forced | cluster4 | 245.339 | 628.363 | -156.12 | 0/5 | 2.19/4.55 | 18.00/22.12 | 1991.1 |
| 864 | 8 | no-rebuild | soa | 255.302 | 368.125 | -44.19 | 0/5 | 2.23/2.34 | 18.14/18.16 | 40.5 |
| 864 | 8 | no-rebuild | aosoa4 | 255.698 | 372.565 | -45.71 | 0/5 | 2.22/2.37 | 18.14/18.16 | 40.5 |
| 864 | 8 | no-rebuild | cluster4 | 256.971 | 427.595 | -66.40 | 0/5 | 2.21/4.54 | 18.14/22.19 | 1994.1 |
| 864 | 8 | forced | soa | 246.257 | 279.449 | -13.48 | 0/5 | 2.24/2.39 | 18.12/18.16 | 40.5 |
| 864 | 8 | forced | aosoa4 | 243.804 | 276.517 | -13.42 | 0/5 | 2.24/2.36 | 18.14/18.19 | 40.5 |
| 864 | 8 | forced | cluster4 | 243.695 | 627.880 | -157.65 | 0/5 | 2.22/4.51 | 18.16/22.23 | 1991.1 |
| 864 | 11d | no-rebuild | soa | 255.806 | 382.813 | -49.65 | 0/5 | 2.29/2.38 | 18.22/18.27 | 40.5 |
| 864 | 11d | no-rebuild | aosoa4 | 254.850 | 376.293 | -47.65 | 0/5 | 2.24/2.36 | 18.20/18.25 | 40.5 |
| 864 | 11d | no-rebuild | cluster4 | 254.306 | 426.678 | -67.78 | 0/5 | 2.24/4.52 | 18.22/22.28 | 1994.1 |
| 864 | 11d | forced | soa | 243.765 | 275.299 | -12.94 | 0/5 | 2.23/2.34 | 18.23/18.27 | 40.5 |
| 864 | 11d | forced | aosoa4 | 243.083 | 278.788 | -14.69 | 0/5 | 2.22/2.37 | 18.19/18.30 | 40.5 |
| 864 | 11d | forced | cluster4 | 244.159 | 628.155 | -157.27 | 0/5 | 2.22/4.53 | 18.20/22.34 | 1991.1 |
| 2048 | 1 | no-rebuild | soa | 240.053 | 282.481 | -17.67 | 0/5 | 8.31/8.58 | 38.55/38.64 | 96.0 |
| 2048 | 1 | no-rebuild | aosoa4 | 238.928 | 275.027 | -15.11 | 0/5 | 8.34/8.54 | 38.56/38.64 | 96.0 |
| 2048 | 1 | no-rebuild | cluster4 | 238.418 | 329.822 | -38.34 | 0/5 | 8.40/13.82 | 38.55/47.78 | 4622.7 |
| 2048 | 1 | forced | soa | 249.738 | 254.511 | -1.91 | 0/5 | 8.31/8.47 | 38.55/38.64 | 96.0 |
| 2048 | 1 | forced | aosoa4 | 249.668 | 254.754 | -2.04 | 0/5 | 8.26/8.47 | 38.56/38.64 | 96.0 |
| 2048 | 1 | forced | cluster4 | 255.099 | 445.546 | -74.66 | 0/5 | 8.56/13.94 | 38.56/47.94 | 4620.9 |
| 2048 | 4 | no-rebuild | soa | 245.268 | 281.771 | -14.88 | 0/5 | 8.57/8.63 | 38.62/38.69 | 96.0 |
| 2048 | 4 | no-rebuild | aosoa4 | 240.375 | 274.757 | -14.30 | 0/5 | 8.33/8.43 | 38.62/38.70 | 96.0 |
| 2048 | 4 | no-rebuild | cluster4 | 241.128 | 326.091 | -35.24 | 0/5 | 8.36/13.47 | 38.61/47.86 | 4622.7 |
| 2048 | 4 | forced | soa | 250.139 | 253.995 | -1.54 | 0/5 | 8.37/8.44 | 38.62/38.69 | 96.0 |
| 2048 | 4 | forced | aosoa4 | 261.380 | 266.793 | -2.07 | 0/5 | 8.64/8.82 | 38.62/38.69 | 96.0 |
| 2048 | 4 | forced | cluster4 | 261.098 | 453.218 | -73.58 | 0/5 | 8.76/14.56 | 38.64/48.00 | 4620.9 |
| 2048 | 8 | no-rebuild | soa | 253.757 | 294.522 | -16.06 | 0/5 | 8.85/8.85 | 38.75/38.81 | 96.0 |
| 2048 | 8 | no-rebuild | aosoa4 | 252.785 | 291.964 | -15.50 | 0/5 | 8.77/8.84 | 38.77/38.81 | 96.0 |
| 2048 | 8 | no-rebuild | cluster4 | 249.685 | 336.882 | -34.92 | 0/5 | 8.74/13.99 | 38.75/47.98 | 4622.7 |
| 2048 | 8 | forced | soa | 261.406 | 266.043 | -1.77 | 1/5 | 8.81/8.90 | 38.75/38.83 | 96.0 |
| 2048 | 8 | forced | aosoa4 | 257.047 | 263.112 | -2.36 | 1/5 | 8.66/8.72 | 38.72/38.81 | 96.0 |
| 2048 | 8 | forced | cluster4 | 262.416 | 451.166 | -71.93 | 0/5 | 8.71/13.83 | 38.77/48.11 | 4620.9 |
| 2048 | 11d | no-rebuild | soa | 248.977 | 287.379 | -15.42 | 0/5 | 8.72/8.80 | 38.84/38.89 | 96.0 |
| 2048 | 11d | no-rebuild | aosoa4 | 248.490 | 286.601 | -15.34 | 0/5 | 8.63/8.71 | 38.81/38.91 | 96.0 |
| 2048 | 11d | no-rebuild | cluster4 | 254.537 | 335.984 | -32.00 | 0/5 | 8.78/14.12 | 38.83/48.09 | 4622.7 |
| 2048 | 11d | forced | soa | 261.452 | 265.753 | -1.65 | 2/5 | 8.75/8.93 | 38.84/38.89 | 96.0 |
| 2048 | 11d | forced | aosoa4 | 260.238 | 265.001 | -1.83 | 1/5 | 8.74/8.85 | 38.84/38.91 | 96.0 |
| 2048 | 11d | forced | cluster4 | 264.863 | 450.278 | -70.00 | 0/5 | 8.81/14.03 | 38.81/48.25 | 4620.9 |
| 5000 | 1 | no-rebuild | soa | 265.022 | 297.214 | -12.15 | 0/5 | 18.40/18.56 | 88.30/88.56 | 234.4 |
| 5000 | 1 | no-rebuild | aosoa4 | 261.928 | 301.041 | -14.93 | 0/5 | 18.50/18.89 | 88.31/88.58 | 234.4 |
| 5000 | 1 | no-rebuild | cluster4 | 264.788 | 382.106 | -44.31 | 0/5 | 18.46/32.34 | 88.33/110.56 | 11021.0 |
| 5000 | 1 | forced | soa | 262.811 | 272.073 | -3.52 | 0/5 | 18.44/18.85 | 88.48/88.66 | 234.4 |
| 5000 | 1 | forced | aosoa4 | 265.347 | 271.693 | -2.39 | 1/5 | 18.64/18.87 | 88.48/88.66 | 234.4 |
| 5000 | 1 | forced | cluster4 | 263.968 | 504.165 | -90.99 | 0/5 | 18.48/31.51 | 88.50/111.02 | 11021.4 |
| 5000 | 4 | no-rebuild | soa | 254.405 | 288.269 | -13.31 | 0/5 | 10.61/10.97 | 88.58/88.80 | 234.4 |
| 5000 | 4 | no-rebuild | aosoa4 | 254.412 | 289.502 | -13.79 | 0/5 | 10.60/11.25 | 88.58/88.78 | 234.4 |
| 5000 | 4 | no-rebuild | cluster4 | 253.641 | 363.948 | -43.49 | 0/5 | 10.65/23.50 | 88.59/110.78 | 11021.0 |
| 5000 | 4 | forced | soa | 115.767 | 121.164 | -4.66 | 0/5 | 10.61/10.93 | 88.70/88.86 | 234.4 |
| 5000 | 4 | forced | aosoa4 | 115.171 | 122.031 | -5.96 | 0/5 | 10.65/11.00 | 88.70/88.88 | 234.4 |
| 5000 | 4 | forced | cluster4 | 119.652 | 359.462 | -200.42 | 0/5 | 11.04/24.53 | 88.73/111.30 | 11021.4 |
| 5000 | 8 | no-rebuild | soa | 265.244 | 300.998 | -13.48 | 0/5 | 10.53/10.48 | 90.27/90.39 | 234.4 |
| 5000 | 8 | no-rebuild | aosoa4 | 265.132 | 297.860 | -12.34 | 0/5 | 10.40/10.55 | 90.02/90.23 | 234.4 |
| 5000 | 8 | no-rebuild | cluster4 | 260.932 | 369.679 | -41.68 | 0/5 | 10.05/23.82 | 90.08/112.52 | 11021.0 |
| 5000 | 8 | forced | soa | 107.159 | 112.399 | -4.89 | 1/5 | 10.39/10.84 | 90.61/90.95 | 234.4 |
| 5000 | 8 | forced | aosoa4 | 106.222 | 111.554 | -5.02 | 0/5 | 10.25/10.70 | 90.45/90.69 | 234.4 |
| 5000 | 8 | forced | cluster4 | 107.018 | 348.254 | -225.42 | 0/5 | 10.52/23.81 | 90.58/113.09 | 11021.4 |
| 5000 | 11d | no-rebuild | soa | 266.086 | 298.203 | -12.07 | 0/5 | 10.15/10.78 | 90.62/91.12 | 234.4 |
| 5000 | 11d | no-rebuild | aosoa4 | 261.200 | 291.515 | -11.61 | 0/5 | 10.08/10.39 | 90.61/90.89 | 234.4 |
| 5000 | 11d | no-rebuild | cluster4 | 253.803 | 373.344 | -47.10 | 0/5 | 9.61/22.56 | 90.61/112.47 | 11021.0 |
| 5000 | 11d | forced | soa | 97.162 | 101.850 | -4.82 | 0/5 | 9.74/9.93 | 91.08/91.23 | 234.4 |
| 5000 | 11d | forced | aosoa4 | 97.797 | 104.026 | -6.37 | 0/5 | 9.76/10.15 | 90.95/91.17 | 234.4 |
| 5000 | 11d | forced | cluster4 | 95.738 | 324.828 | -239.29 | 0/5 | 9.64/22.72 | 90.98/113.69 | 11021.4 |
| 10000 | 1 | no-rebuild | soa | 257.471 | 293.425 | -13.96 | 1/5 | 34.71/34.90 | 178.48/178.91 | 468.8 |
| 10000 | 1 | no-rebuild | aosoa4 | 256.746 | 288.047 | -12.19 | 0/5 | 33.86/34.42 | 178.48/178.91 | 468.8 |
| 10000 | 1 | no-rebuild | cluster4 | 257.223 | 355.446 | -38.19 | 0/5 | 34.00/60.97 | 178.48/215.64 | 21544.8 |
| 10000 | 1 | forced | soa | 476.872 | 490.989 | -2.96 | 0/5 | 34.13/34.85 | 179.09/179.41 | 468.8 |
| 10000 | 1 | forced | aosoa4 | 477.517 | 490.882 | -2.80 | 0/5 | 34.12/34.61 | 179.08/179.42 | 468.8 |
| 10000 | 1 | forced | cluster4 | 476.517 | 954.023 | -100.21 | 0/5 | 33.83/60.85 | 179.08/222.52 | 21537.4 |
| 10000 | 4 | no-rebuild | soa | 259.029 | 285.186 | -10.10 | 0/5 | 20.41/20.93 | 173.09/173.45 | 468.8 |
| 10000 | 4 | no-rebuild | aosoa4 | 257.298 | 286.874 | -11.49 | 0/5 | 20.37/20.91 | 173.06/173.45 | 468.8 |
| 10000 | 4 | no-rebuild | cluster4 | 259.444 | 355.946 | -37.20 | 0/5 | 20.30/47.28 | 173.08/215.91 | 21544.8 |
| 10000 | 4 | forced | soa | 218.299 | 229.884 | -5.31 | 0/5 | 20.12/21.00 | 174.16/174.11 | 468.8 |
| 10000 | 4 | forced | aosoa4 | 218.949 | 231.159 | -5.58 | 0/5 | 20.13/21.38 | 174.12/173.92 | 468.8 |
| 10000 | 4 | forced | cluster4 | 219.448 | 691.602 | -215.15 | 0/5 | 20.55/46.82 | 174.12/217.81 | 21537.4 |
| 10000 | 8 | no-rebuild | soa | 257.222 | 287.567 | -11.80 | 0/5 | 19.07/19.73 | 175.52/175.95 | 468.8 |
| 10000 | 8 | no-rebuild | aosoa4 | 257.030 | 287.263 | -11.76 | 0/5 | 19.23/19.87 | 175.36/176.02 | 468.8 |
| 10000 | 8 | no-rebuild | cluster4 | 256.762 | 358.188 | -39.50 | 0/5 | 19.15/46.73 | 175.50/218.31 | 21544.8 |
| 10000 | 8 | forced | soa | 193.315 | 206.538 | -6.84 | 0/5 | 19.19/19.73 | 176.06/176.53 | 468.8 |
| 10000 | 8 | forced | aosoa4 | 194.388 | 206.236 | -6.09 | 1/5 | 19.43/19.80 | 176.31/176.53 | 468.8 |
| 10000 | 8 | forced | cluster4 | 193.974 | 665.784 | -243.23 | 0/5 | 19.33/46.81 | 176.30/220.09 | 21537.4 |
| 10000 | 11d | no-rebuild | soa | 257.097 | 283.950 | -10.44 | 0/5 | 18.54/19.05 | 176.91/177.28 | 468.8 |
| 10000 | 11d | no-rebuild | aosoa4 | 256.726 | 289.982 | -12.95 | 0/5 | 18.53/19.11 | 176.98/177.09 | 468.8 |
| 10000 | 11d | no-rebuild | cluster4 | 256.729 | 355.883 | -38.62 | 0/5 | 18.75/45.97 | 176.55/219.62 | 21544.8 |
| 10000 | 11d | forced | soa | 182.463 | 193.904 | -6.27 | 0/5 | 18.47/19.55 | 178.59/178.58 | 468.8 |
| 10000 | 11d | forced | aosoa4 | 183.254 | 195.342 | -6.60 | 0/5 | 18.40/19.13 | 178.23/178.53 | 468.8 |
| 10000 | 11d | forced | cluster4 | 182.520 | 653.959 | -258.29 | 0/5 | 18.88/46.01 | 177.88/222.17 | 21537.4 |

## Direct audit and snapshot-free protocol

The exact primary source was preserved after timing in commit
`4b63e9b12a6bc7ea3c462c9a5db64c0928b4d86b`; its SHA matches the launch copy.
This commit correspondence does not change the actual launch HEAD or dirty
status. The original audit used a tolerant canonical-state comparison after its
first oracle. Before acceptance, a real subprocess regression demonstrated that
a 0.75-tolerance canonical force bias followed by a 1.5-tolerance candidate bias
could pass both hops. The fixture updates Fmax and recomputes consistent FNV, so
its failure is physical. The later direct-oracle validator rejects that case
without changing any tolerance. Its full source and hash are preserved as a
separate post-timing audit epoch.

The direct audit also independently recomputes FNV from every actual saved
state and SHA256 from every complete initial input. It verifies raw stdout,
stderr, compressed snapshot and original uncompressed snapshot SHA256. It emits
`secondary-references.json`, tied to the audit and binary hashes. The audit
process exits before the secondary process starts, releasing its arrays and
neighbor sets. The secondary loads only small reference rows and fingerprints.

The separately authorized `secondary` action runs the same complete matrix,
steps, timestep, AB/BA pairs, warmups and production controls, omitting only
`--snapshot`. It requires the same executable SHA256. Before timing it freezes
its own complete orchestrator/test source, compiler/source identity, copied
configuration and audited references with explicit input hashes. Every available
stdout physics field, counter, type, domain/identity field and FNV is compared
against the same-mode/parameter audited references. Timing values are separately
checked for finite/nonnegative values; physical values use the fixed tolerances,
and every other available field/type/counter must match exactly. No vectors or full-state SHA are
emitted by the secondary child. Its recomputed-reference FNV is a consistency
check, not a cryptographic full-state proof or a replacement for primary arrays.

The primary and secondary medians, durations and resource rows remain separate.
Conditional five-pair reproductions are retained for any secondary case that
meets the original timing screen; no candidate is promoted by either script.

The authoritative direct audit source SHA256 is
`c1a312eef14a74aa5769b09c374c459d198b88397923ccd0f64db18293f03d78`.
It passed all 1,344 saved invocations using 13 exact-coordinate oracle contexts;
maximum force/energy/virial errors were 1.6362e-14, 1.6953e-9 and 5.7917e-9.
All 4,032 saved state fingerprints were independently recomputed; each of the
160 mode/parameter reference groups had one stable final fingerprint. Actual
initial-input SHA256 values, derived during this post-timing audit, were:

| N | Actual complete initial-input SHA256 |
|---:|:---|
| 864 | e30da513aaeee6f11cb85ad3282dd83aa57fa1577316c891f83d25599220e9fa |
| 2048 | fde06a3b57b9d6895153447b1743806709e019ce60be2492ed8d492007b4837a |
| 5000 | 87785b1b54f158eecb02af5cbd4dd3c79ce6caaef7713a349df41cfcb8231784 |
| 10000 | 3f97a1b9ae562df9cdc2816944249b2123f0fa312eee261bb5eacac662e93852 |

The primary and already-running audit did not record Python identity at startup.
A timestamped later observation resolved the worktree python3 to Python 3.14.5;
that observation is not retroactively attributed to their launch metadata. The
secondary/future audit startup explicitly captures sys.executable/sys.version.
The seven driver tests passed under both Python 3.14.5 and CI-compatible
Python 3.12.8 before secondary timing. The biased-force regression had a genuine
RED with the original audit accepting the bad outputs, then GREEN with direct
oracle checking. No numerical tolerance or Rust source changed.

## Snapshot-free results

The separate protocol retained 1,344 invocations. All available-field and independently recomputed-reference FNV gates passed. None of its 96 measured comparisons passed the original timing bar; no secondary positive-case reproduction was triggered.

Its best median time saving was -0.9117% (aosoa4, N5000, requested pool1, forced), 0/5 favorable pairs. Primary and secondary observations are separate; neither provides a 5% positive case.

| N | Pool | Policy | B | A VV ms | B VV ms | Saved % | Favorable | Setup A/B ms | RSS A/B MiB | Packed B KiB |
|---:|---:|:---|:---|---:|---:|---:|---:|---:|---:|---:|
| 864 | 1 | no-rebuild | soa | 264.502 | 383.356 | -44.94 | 0/5 | 1.53/1.66 | 2.84/2.88 | 40.5 |
| 864 | 1 | no-rebuild | aosoa4 | 266.031 | 389.069 | -46.25 | 0/5 | 1.55/1.67 | 2.83/2.86 | 40.5 |
| 864 | 1 | no-rebuild | cluster4 | 264.458 | 435.838 | -64.80 | 0/5 | 1.54/3.91 | 2.84/6.88 | 1994.1 |
| 864 | 1 | forced | soa | 288.085 | 307.906 | -6.88 | 0/5 | 1.71/1.85 | 2.84/2.86 | 40.5 |
| 864 | 1 | forced | aosoa4 | 250.804 | 286.067 | -14.06 | 0/5 | 1.52/1.66 | 2.83/2.86 | 40.5 |
| 864 | 1 | forced | cluster4 | 253.444 | 645.450 | -154.67 | 0/5 | 1.52/3.92 | 2.84/8.08 | 1991.1 |
| 864 | 4 | no-rebuild | soa | 267.141 | 383.956 | -43.73 | 0/5 | 1.56/1.75 | 2.91/2.94 | 40.5 |
| 864 | 4 | no-rebuild | aosoa4 | 264.909 | 390.821 | -47.53 | 0/5 | 1.56/1.70 | 2.92/2.95 | 40.5 |
| 864 | 4 | no-rebuild | cluster4 | 264.535 | 431.718 | -63.20 | 0/5 | 1.56/3.95 | 2.91/6.97 | 1994.1 |
| 864 | 4 | forced | soa | 251.102 | 282.791 | -12.62 | 0/5 | 1.50/1.68 | 2.91/2.94 | 40.5 |
| 864 | 4 | forced | aosoa4 | 250.824 | 285.096 | -13.66 | 0/5 | 1.52/1.67 | 2.91/2.95 | 40.5 |
| 864 | 4 | forced | cluster4 | 253.161 | 647.913 | -155.93 | 0/5 | 1.51/3.88 | 2.92/8.19 | 1991.1 |
| 864 | 8 | no-rebuild | soa | 264.488 | 383.344 | -44.94 | 0/5 | 1.53/1.71 | 3.05/3.14 | 40.5 |
| 864 | 8 | no-rebuild | aosoa4 | 264.391 | 389.262 | -47.23 | 0/5 | 1.55/1.73 | 3.03/3.11 | 40.5 |
| 864 | 8 | no-rebuild | cluster4 | 266.589 | 440.998 | -65.42 | 0/5 | 1.57/3.97 | 3.05/7.08 | 1994.1 |
| 864 | 8 | forced | soa | 252.381 | 286.819 | -13.65 | 0/5 | 1.55/1.71 | 3.06/3.09 | 40.5 |
| 864 | 8 | forced | aosoa4 | 271.988 | 289.770 | -6.54 | 0/5 | 1.77/1.81 | 3.05/3.06 | 40.5 |
| 864 | 8 | forced | cluster4 | 250.598 | 648.007 | -158.58 | 0/5 | 1.52/3.89 | 3.06/8.30 | 1991.1 |
| 864 | 11d | no-rebuild | soa | 261.941 | 380.387 | -45.22 | 0/5 | 1.51/1.71 | 3.22/3.23 | 40.5 |
| 864 | 11d | no-rebuild | aosoa4 | 266.897 | 390.408 | -46.28 | 0/5 | 1.56/1.73 | 3.16/3.20 | 40.5 |
| 864 | 11d | no-rebuild | cluster4 | 265.549 | 448.172 | -68.77 | 0/5 | 1.54/4.00 | 3.14/7.22 | 1994.1 |
| 864 | 11d | forced | soa | 252.409 | 286.772 | -13.61 | 0/5 | 1.55/1.73 | 3.12/3.22 | 40.5 |
| 864 | 11d | forced | aosoa4 | 253.008 | 285.809 | -12.96 | 0/5 | 1.53/1.78 | 3.17/3.22 | 40.5 |
| 864 | 11d | forced | cluster4 | 254.584 | 659.569 | -159.08 | 0/5 | 1.58/4.15 | 3.14/8.48 | 1991.1 |
| 2048 | 1 | no-rebuild | soa | 246.765 | 283.560 | -14.91 | 0/5 | 6.87/6.99 | 3.77/3.86 | 96.0 |
| 2048 | 1 | no-rebuild | aosoa4 | 249.460 | 284.140 | -13.90 | 0/5 | 6.86/6.94 | 3.78/3.86 | 96.0 |
| 2048 | 1 | no-rebuild | cluster4 | 245.502 | 331.246 | -34.93 | 0/5 | 6.80/12.16 | 3.78/13.00 | 4622.7 |
| 2048 | 1 | forced | soa | 258.758 | 266.251 | -2.90 | 0/5 | 6.83/7.05 | 3.77/3.86 | 96.0 |
| 2048 | 1 | forced | aosoa4 | 263.839 | 266.802 | -1.12 | 0/5 | 6.97/7.14 | 3.78/3.86 | 96.0 |
| 2048 | 1 | forced | cluster4 | 264.197 | 452.537 | -71.29 | 0/5 | 6.93/12.30 | 3.77/15.64 | 4620.9 |
| 2048 | 4 | no-rebuild | soa | 248.324 | 287.899 | -15.94 | 0/5 | 6.87/6.98 | 3.86/3.94 | 96.0 |
| 2048 | 4 | no-rebuild | aosoa4 | 248.299 | 286.906 | -15.55 | 0/5 | 7.04/7.07 | 3.84/3.95 | 96.0 |
| 2048 | 4 | no-rebuild | cluster4 | 248.574 | 332.024 | -33.57 | 0/5 | 6.86/12.13 | 3.84/13.06 | 4622.7 |
| 2048 | 4 | forced | soa | 258.740 | 266.150 | -2.86 | 0/5 | 6.87/7.01 | 3.84/3.94 | 96.0 |
| 2048 | 4 | forced | aosoa4 | 255.209 | 260.682 | -2.14 | 0/5 | 6.74/6.99 | 3.84/3.94 | 96.0 |
| 2048 | 4 | forced | cluster4 | 256.137 | 445.034 | -73.75 | 0/5 | 6.83/12.12 | 3.86/15.72 | 4620.9 |
| 2048 | 8 | no-rebuild | soa | 247.749 | 287.429 | -16.02 | 0/5 | 6.93/7.11 | 4.03/4.09 | 96.0 |
| 2048 | 8 | no-rebuild | aosoa4 | 248.213 | 286.581 | -15.46 | 1/5 | 6.92/7.00 | 4.00/4.08 | 96.0 |
| 2048 | 8 | no-rebuild | cluster4 | 249.182 | 335.826 | -34.77 | 1/5 | 6.90/12.29 | 3.98/13.28 | 4622.7 |
| 2048 | 8 | forced | soa | 287.194 | 297.617 | -3.63 | 1/5 | 7.77/7.91 | 4.00/4.06 | 96.0 |
| 2048 | 8 | forced | aosoa4 | 255.099 | 260.036 | -1.94 | 0/5 | 6.81/6.93 | 4.02/4.06 | 96.0 |
| 2048 | 8 | forced | cluster4 | 256.121 | 444.055 | -73.38 | 0/5 | 6.82/12.19 | 4.00/15.84 | 4620.9 |
| 2048 | 11d | no-rebuild | soa | 248.184 | 287.097 | -15.68 | 0/5 | 6.84/6.98 | 4.08/4.22 | 96.0 |
| 2048 | 11d | no-rebuild | aosoa4 | 249.751 | 286.051 | -14.53 | 0/5 | 6.86/7.10 | 4.11/4.22 | 96.0 |
| 2048 | 11d | no-rebuild | cluster4 | 244.546 | 327.152 | -33.78 | 0/5 | 6.85/12.11 | 4.17/13.36 | 4622.7 |
| 2048 | 11d | forced | soa | 260.163 | 267.766 | -2.92 | 1/5 | 6.90/7.03 | 4.12/4.25 | 96.0 |
| 2048 | 11d | forced | aosoa4 | 259.245 | 264.406 | -1.99 | 0/5 | 6.93/7.06 | 4.08/4.17 | 96.0 |
| 2048 | 11d | forced | cluster4 | 259.081 | 453.283 | -74.96 | 0/5 | 6.88/12.21 | 4.12/16.02 | 4620.9 |
| 5000 | 1 | no-rebuild | soa | 264.940 | 296.066 | -11.75 | 0/5 | 14.47/14.74 | 11.44/11.44 | 234.4 |
| 5000 | 1 | no-rebuild | aosoa4 | 269.916 | 304.870 | -12.95 | 0/5 | 14.59/14.88 | 11.44/11.44 | 234.4 |
| 5000 | 1 | no-rebuild | cluster4 | 264.861 | 377.999 | -42.72 | 0/5 | 14.58/27.96 | 11.44/33.22 | 11021.0 |
| 5000 | 1 | forced | soa | 259.299 | 267.663 | -3.23 | 0/5 | 14.29/14.68 | 11.75/11.89 | 234.4 |
| 5000 | 1 | forced | aosoa4 | 265.229 | 267.647 | -0.91 | 0/5 | 14.33/14.53 | 11.78/11.89 | 234.4 |
| 5000 | 1 | forced | cluster4 | 259.643 | 498.046 | -91.82 | 0/5 | 14.27/27.55 | 11.77/39.98 | 11021.4 |
| 5000 | 4 | no-rebuild | soa | 266.640 | 296.921 | -11.36 | 0/5 | 7.02/7.32 | 12.22/12.03 | 234.4 |
| 5000 | 4 | no-rebuild | aosoa4 | 264.329 | 299.774 | -13.41 | 0/5 | 6.95/7.40 | 12.20/12.08 | 234.4 |
| 5000 | 4 | no-rebuild | cluster4 | 266.186 | 377.898 | -41.97 | 0/5 | 7.08/20.41 | 11.84/33.78 | 11021.0 |
| 5000 | 4 | forced | soa | 121.392 | 126.180 | -3.94 | 0/5 | 7.03/7.45 | 12.55/12.80 | 234.4 |
| 5000 | 4 | forced | aosoa4 | 120.026 | 127.033 | -5.84 | 0/5 | 7.06/7.32 | 12.58/12.67 | 234.4 |
| 5000 | 4 | forced | cluster4 | 119.938 | 360.070 | -200.21 | 0/5 | 6.96/20.26 | 12.56/41.08 | 11021.4 |
| 5000 | 8 | no-rebuild | soa | 266.123 | 327.231 | -22.96 | 0/5 | 6.33/6.70 | 12.55/12.44 | 234.4 |
| 5000 | 8 | no-rebuild | aosoa4 | 261.674 | 296.407 | -13.27 | 0/5 | 6.18/6.61 | 12.56/12.36 | 234.4 |
| 5000 | 8 | no-rebuild | cluster4 | 261.367 | 371.565 | -42.16 | 0/5 | 6.27/19.37 | 12.66/34.19 | 11021.0 |
| 5000 | 8 | forced | soa | 106.505 | 111.465 | -4.66 | 0/5 | 6.29/6.53 | 13.41/13.23 | 234.4 |
| 5000 | 8 | forced | aosoa4 | 104.733 | 111.283 | -6.25 | 0/5 | 6.24/6.57 | 13.64/13.52 | 234.4 |
| 5000 | 8 | forced | cluster4 | 104.264 | 341.766 | -227.79 | 0/5 | 6.21/19.31 | 13.69/41.53 | 11021.4 |
| 5000 | 11d | no-rebuild | soa | 265.410 | 295.967 | -11.51 | 0/5 | 6.09/6.37 | 12.30/12.75 | 234.4 |
| 5000 | 11d | no-rebuild | aosoa4 | 264.981 | 301.293 | -13.70 | 0/5 | 5.94/6.46 | 12.62/12.31 | 234.4 |
| 5000 | 11d | no-rebuild | cluster4 | 261.501 | 366.943 | -40.32 | 0/5 | 5.96/19.37 | 12.64/34.30 | 11021.0 |
| 5000 | 11d | forced | soa | 101.901 | 108.055 | -6.04 | 0/5 | 5.83/6.24 | 13.78/14.42 | 234.4 |
| 5000 | 11d | forced | aosoa4 | 103.382 | 110.198 | -6.59 | 0/5 | 6.10/6.42 | 14.62/14.86 | 234.4 |
| 5000 | 11d | forced | cluster4 | 102.504 | 344.286 | -235.88 | 0/5 | 6.15/19.43 | 14.19/43.05 | 11021.4 |
| 10000 | 1 | no-rebuild | soa | 267.400 | 296.873 | -11.02 | 0/5 | 26.93/27.50 | 17.16/17.17 | 468.8 |
| 10000 | 1 | no-rebuild | aosoa4 | 268.546 | 301.024 | -12.09 | 0/5 | 26.95/27.59 | 17.17/17.16 | 468.8 |
| 10000 | 1 | no-rebuild | cluster4 | 267.960 | 368.249 | -37.43 | 0/5 | 27.04/55.35 | 17.17/53.64 | 21544.8 |
| 10000 | 1 | forced | soa | 509.204 | 519.401 | -2.00 | 1/5 | 27.12/31.00 | 17.81/18.30 | 468.8 |
| 10000 | 1 | forced | aosoa4 | 488.480 | 501.172 | -2.60 | 0/5 | 26.75/27.43 | 17.81/18.28 | 468.8 |
| 10000 | 1 | forced | cluster4 | 498.250 | 993.073 | -99.31 | 0/5 | 27.04/55.29 | 17.81/73.89 | 21537.4 |
| 10000 | 4 | no-rebuild | soa | 269.093 | 297.332 | -10.49 | 0/5 | 13.42/13.74 | 18.61/18.59 | 468.8 |
| 10000 | 4 | no-rebuild | aosoa4 | 268.488 | 299.137 | -11.42 | 0/5 | 13.16/13.74 | 20.98/18.34 | 468.8 |
| 10000 | 4 | no-rebuild | cluster4 | 268.567 | 373.095 | -38.92 | 0/5 | 13.08/41.81 | 18.58/60.59 | 21544.8 |
| 10000 | 4 | forced | soa | 225.425 | 235.723 | -4.57 | 0/5 | 13.09/13.69 | 23.42/24.34 | 468.8 |
| 10000 | 4 | forced | aosoa4 | 225.291 | 238.178 | -5.72 | 0/5 | 13.31/13.67 | 23.47/24.36 | 468.8 |
| 10000 | 4 | forced | cluster4 | 224.388 | 715.453 | -218.85 | 0/5 | 13.00/41.30 | 23.00/78.58 | 21537.4 |
| 10000 | 8 | no-rebuild | soa | 269.370 | 296.963 | -10.24 | 0/5 | 11.95/12.54 | 20.23/20.22 | 468.8 |
| 10000 | 8 | no-rebuild | aosoa4 | 268.947 | 301.441 | -12.08 | 0/5 | 11.77/12.60 | 19.16/19.69 | 468.8 |
| 10000 | 8 | no-rebuild | cluster4 | 267.819 | 372.455 | -39.07 | 0/5 | 11.81/39.55 | 19.61/61.89 | 21544.8 |
| 10000 | 8 | forced | soa | 201.319 | 212.037 | -5.32 | 0/5 | 11.84/12.48 | 23.33/23.38 | 468.8 |
| 10000 | 8 | forced | aosoa4 | 202.534 | 217.147 | -7.22 | 0/5 | 12.04/12.61 | 22.97/23.02 | 468.8 |
| 10000 | 8 | forced | cluster4 | 201.156 | 697.405 | -246.70 | 0/5 | 11.66/40.81 | 23.42/77.45 | 21537.4 |
| 10000 | 11d | no-rebuild | soa | 263.935 | 291.279 | -10.36 | 0/5 | 11.14/11.79 | 20.77/20.67 | 468.8 |
| 10000 | 11d | no-rebuild | aosoa4 | 268.374 | 300.041 | -11.80 | 0/5 | 11.67/12.07 | 20.45/21.41 | 468.8 |
| 10000 | 11d | no-rebuild | cluster4 | 272.404 | 384.548 | -41.17 | 0/5 | 11.51/39.98 | 20.09/63.14 | 21544.8 |
| 10000 | 11d | forced | soa | 197.118 | 208.934 | -5.99 | 0/5 | 11.57/12.28 | 27.78/28.78 | 468.8 |
| 10000 | 11d | forced | aosoa4 | 196.647 | 208.934 | -6.25 | 0/5 | 11.53/12.09 | 25.44/28.86 | 468.8 |
| 10000 | 11d | forced | cluster4 | 196.303 | 693.030 | -253.04 | 0/5 | 11.57/39.63 | 25.81/83.12 | 21537.4 |

## Memory and regression assessment

No candidate has a case meeting the fixed time-saving bar in either protocol. The process high-water ranges below include every warmup, measured and control row. Primary includes retained snapshot Values and output allocation; secondary omits those allocations. Packed capacity is partial telemetry and excludes core neighbor/workspace and transient reconstruction maps. No new memory acceptance threshold is invented.

| N | Mode | Primary RSS min/max MiB | Secondary RSS min/max MiB | Packed capacity min/max KiB |
|---:|:---|---:|---:|---:|
| 864 | interleaved-serial | 17.88/18.33 | 2.83/3.28 | 0.0/0.0 |
| 864 | interleaved | 17.89/18.36 | 2.83/3.25 | 0.0/0.0 |
| 864 | soa | 17.94/18.33 | 2.86/3.28 | 40.5/40.5 |
| 864 | aosoa4 | 17.94/18.34 | 2.86/3.27 | 40.5/40.5 |
| 864 | cluster4 | 21.94/24.56 | 6.88/10.67 | 1991.1/1994.1 |
| 2048 | interleaved-serial | 38.52/38.94 | 3.77/4.22 | 0.0/0.0 |
| 2048 | interleaved | 38.53/38.89 | 3.77/4.14 | 0.0/0.0 |
| 2048 | soa | 38.59/38.97 | 3.84/4.31 | 96.0/96.0 |
| 2048 | aosoa4 | 38.61/38.95 | 3.84/4.30 | 96.0/96.0 |
| 2048 | cluster4 | 47.73/48.27 | 13.00/16.05 | 4620.9/4622.7 |
| 5000 | interleaved-serial | 88.28/92.42 | 11.03/14.95 | 0.0/0.0 |
| 5000 | interleaved | 88.28/94.69 | 11.30/16.38 | 0.0/0.0 |
| 5000 | soa | 88.53/92.50 | 11.02/14.89 | 234.4/234.4 |
| 5000 | aosoa4 | 88.56/93.23 | 11.44/15.78 | 234.4/234.4 |
| 5000 | cluster4 | 110.50/114.98 | 32.83/43.69 | 11021.0/11021.4 |
| 10000 | interleaved-serial | 172.98/179.20 | 17.16/31.47 | 0.0/0.0 |
| 10000 | interleaved | 173.98/184.03 | 17.16/30.81 | 0.0/0.0 |
| 10000 | soa | 173.39/179.45 | 17.16/32.56 | 468.8/468.8 |
| 10000 | aosoa4 | 173.42/179.44 | 17.16/31.14 | 468.8/468.8 |
| 10000 | cluster4 | 215.64/225.72 | 53.64/86.52 | 21537.4/21544.8 |

| Protocol | Layout | Time-saved range across 32 cases |
|:---|:---|---:|
| primary | soa | -49.65% to -1.54% |
| primary | aosoa4 | -47.65% to -1.83% |
| primary | cluster4 | -258.29% to -32.00% |
| secondary | soa | -45.22% to -2.00% |
| secondary | aosoa4 | -47.53% to -0.91% |
| secondary | cluster4 | -253.04% to -33.57% |

Packing and scatter run on every force evaluation. Cluster sorting/mask construction runs once during setup and on every forced rebuild; its exact timed reconstruction count is steps for forced and zero for no-rebuild. Setup and full-step policy times above charge these costs. There are no exclusive packing/sort timers, so differences are not presented as isolated component timings. The cluster retained buffers and RSS are materially larger; full-step timing already includes that implementation cost.

## Separate production-auto controls

These five-sample medians are descriptive controls and never enter scalar layout acceptance. Each cell also has a separately retained warmup.

| N | Pool | Policy | Primary VV ms | Secondary VV ms | Primary/secondary child CPU ms | Primary/secondary child wall ms | Primary/secondary peak RSS MiB |
|---:|---:|:---|---:|---:|---:|---:|---:|
| 864 | 1 | no-rebuild | 254.084 | 264.371 | 980.37/268.14 | 983.27/269.65 | 17.97/2.86 |
| 864 | 1 | forced | 243.723 | 252.478 | 953.17/256.24 | 956.53/257.61 | 17.97/2.86 |
| 864 | 4 | no-rebuild | 256.275 | 297.435 | 854.00/301.78 | 858.42/303.34 | 18.05/2.92 |
| 864 | 4 | forced | 245.583 | 255.218 | 911.13/259.02 | 914.29/260.31 | 18.03/2.94 |
| 864 | 8 | no-rebuild | 254.936 | 264.966 | 894.21/269.56 | 897.44/270.44 | 18.19/3.08 |
| 864 | 8 | forced | 244.003 | 249.367 | 900.93/254.41 | 903.79/254.92 | 18.22/3.09 |
| 864 | 11d | no-rebuild | 253.687 | 264.691 | 958.17/270.44 | 961.15/270.19 | 18.36/3.19 |
| 864 | 11d | forced | 244.608 | 289.152 | 962.94/295.01 | 966.23/295.11 | 18.28/3.25 |
| 2048 | 1 | no-rebuild | 239.954 | 244.654 | 1818.28/254.40 | 1824.29/255.77 | 38.58/3.78 |
| 2048 | 1 | forced | 253.976 | 257.523 | 1752.80/267.42 | 1761.23/268.82 | 38.56/3.80 |
| 2048 | 4 | no-rebuild | 240.614 | 247.716 | 1665.18/257.85 | 1671.65/259.29 | 38.64/3.89 |
| 2048 | 4 | forced | 260.672 | 255.819 | 1918.53/265.71 | 1924.45/266.86 | 38.64/3.88 |
| 2048 | 8 | no-rebuild | 250.014 | 248.411 | 1841.89/259.34 | 1847.14/260.07 | 38.77/4.06 |
| 2048 | 8 | forced | 259.381 | 258.872 | 1847.04/269.66 | 1851.98/270.38 | 38.78/4.08 |
| 2048 | 11d | no-rebuild | 247.996 | 245.146 | 1858.62/256.60 | 1861.63/256.95 | 38.89/4.14 |
| 2048 | 11d | forced | 263.279 | 259.239 | 1788.97/270.97 | 1807.01/271.02 | 38.86/4.14 |
| 5000 | 1 | no-rebuild | 264.305 | 298.494 | 3972.92/319.64 | 3981.33/321.40 | 88.34/11.45 |
| 5000 | 1 | forced | 265.126 | 262.089 | 4145.19/281.26 | 4157.14/283.26 | 88.52/11.78 |
| 5000 | 4 | no-rebuild | 70.761 | 74.791 | 3857.54/280.26 | 3672.39/85.42 | 89.11/12.22 |
| 5000 | 4 | forced | 77.659 | 77.416 | 4236.34/293.66 | 4040.81/88.04 | 89.98/13.52 |
| 5000 | 8 | no-rebuild | 64.706 | 64.273 | 4157.42/413.33 | 3796.17/74.35 | 94.69/13.20 |
| 5000 | 8 | forced | 67.262 | 60.696 | 4165.97/392.42 | 3822.65/70.59 | 93.84/14.97 |
| 5000 | 11d | no-rebuild | 63.188 | 67.269 | 4240.61/523.93 | 3827.06/76.86 | 91.84/12.72 |
| 5000 | 11d | forced | 56.812 | 57.446 | 4147.76/453.46 | 3784.13/67.19 | 93.70/16.38 |
| 10000 | 1 | no-rebuild | 256.595 | 267.621 | 7988.78/301.19 | 8015.65/302.93 | 178.52/17.17 |
| 10000 | 1 | forced | 479.834 | 496.455 | 7888.48/529.80 | 7923.30/531.86 | 179.08/17.81 |
| 10000 | 4 | no-rebuild | 69.901 | 81.366 | 8266.33/340.69 | 8081.09/100.44 | 174.27/21.00 |
| 10000 | 4 | forced | 136.044 | 139.300 | 8554.36/542.99 | 8195.81/156.48 | 175.28/27.06 |
| 10000 | 8 | no-rebuild | 58.835 | 61.673 | 8272.69/417.56 | 7994.74/77.29 | 177.56/19.94 |
| 10000 | 8 | forced | 104.582 | 108.625 | 8504.83/700.00 | 7988.36/124.22 | 178.77/23.55 |
| 10000 | 11d | no-rebuild | 56.847 | 61.400 | 8694.56/500.54 | 8326.64/76.65 | 180.31/21.95 |
| 10000 | 11d | forced | 90.662 | 106.479 | 8842.86/780.07 | 8265.12/122.48 | 184.03/30.81 |

## Observed counters and scoped decision

| N | Policy | Steps | Initial/final neighbor pairs | Final active pairs | Timed neighbor visits | Packed evaluated pairs | Timed rebuilds / cluster reconstructions |
|---:|:---|---:|---:|---:|---:|---:|---:|
| 864 | no-rebuild | 700 | 33696/33696 | 23328 | 23587200 | 16329600 | 0/0 |
| 864 | forced | 190 | 33696/33696 | 23328 | 6402240 | 4432320 | 190/190 |
| 2048 | no-rebuild | 200 | 79872/79872 | 55296 | 15974400 | 11059200 | 0/0 |
| 2048 | forced | 40 | 79872/79872 | 55296 | 3194880 | 2211840 | 40/40 |
| 5000 | no-rebuild | 90 | 187600/187600 | 129798 | 16884000 | 11681820 | 0/0 |
| 5000 | forced | 20 | 187600/187600 | 129798 | 3752000 | 2595960 | 20/20 |
| 10000 | no-rebuild | 45 | 372400/372400 | 257808 | 16758000 | 11601360 | 0/0 |
| 10000 | forced | 20 | 372400/372400 | 257808 | 7448000 | 5156160 | 20/20 |

Cluster reconstruction is zero in both interleaved controls and the non-cluster layouts; the last table reconstruction column describes Cluster4 only. Setup rebuilds are one in every run, and setup cluster reconstruction is one only for Cluster4. Interleaved evaluated-pair counts remain null/unavailable, never zero.

Reject these measured serial-scalar prototypes for the declared reduced-LJ workloads and keep production unchanged. No explicit SIMD or parallel packed implementation was measured or rejected. The second matrix addresses the material snapshot-retention concern, but its output lacks vectors/full-state SHA and is not a replacement for the complete-state primary. Host user processes were not stopped; observed load is retained and cross-session absolute timings are not a causal measurement of instrumentation alone. No engine/chemistry parity or universal storage-layout result is claimed.

The snapshot-free launch was clean commit `0c54db355b830c49e687e97e2e9aa8502bbaa5d1`, loaded driver SHA256 `fc3baeb1cc309667fadd6b0160bf6de4b682ed3a35f89c5380c40ad1853df58b`, using `/Users/ricardoleal20/.pyenv/versions/3.14.5/bin/python3` (Python 3.14.5). It completed in 6.779 minutes, with binary/source/config/reference hashes unchanged. The final documentation commit does not relabel either measurement epoch. Seven targeted driver tests pass under Python 3.12.8 and 3.14.5; Rust/core source is unchanged from the reviewed Task1 implementation.

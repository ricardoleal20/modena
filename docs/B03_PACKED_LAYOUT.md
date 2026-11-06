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
rtk proxy python3 target/b03-evidence/new-campaign/layout_step.py audit --out target/b03-evidence/new-campaign
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
the final-coordinate skin set. Oracle work is cached by N, steps, timestep and
policy; every other measured snapshot is compared on every physics field with
that independently audited state. Half-pair identities and multiplicity remain
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
clean. Host: Apple M3 Pro, macOS, actual pools 1/4/8/default11. Release compiler
was rustc1.97.1, aarch64-apple-darwin, LLVM22.1.6. Binary SHA256:
`5001bb6b8bdfddc60e74e9fb9bee391e177b1a9992851f47e9ae1a3c8f53c2dc`. Primary loaded driver SHA256:
`5764f6b110841d181f6ce90f670ecca020441aeae5e9d334451fe08df965ef82`.

All 1,344 invocations completed; all 96 timing comparisons failed the 5% bar.
There were no positive cases requiring primary reproduction. The best median
time saving was -1.5416% (SoA, N2048, pool4, forced), with 0/5 favorable pairs.
This is a timing screen: final independent numerical certification and the
authorized separate no-snapshot instrumentation check are still pending.

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

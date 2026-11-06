# B04 implementation and verification record

B04 numerical work is starting in the isolated `ricardo/b04-periodic-neighbors`
branch based on reviewed B03 source `0c54db355b830c49e687e97e2e9aa8502bbaa5d1`.
The core tree is `39b46c7496c6820219056fe38a821378a32633e8`, identical to reviewed
B03 prototype66af0e4. No core/CLI files changed between those revisions.
B03 primary/secondary timing and direct audit finished before this work; its
CPU lease was explicitly released. B03 reports/review remain separate closure
work. Before B04 publication, synchronize its unpublished base with final
reviewed B03 without force-push or an unauthorized merge.

Controller preflight confirmed native managed worktree creation, branch prefix
`ricardo/`, source-tree equivalence and the actual Vec3 dot/wrap operation
orders used by the proof. Subsequent G1 implementation evidence is recorded below.

| Gate | Status | Required evidence |
| --- | --- | --- |
| G1 | executable checks passed; independent Task 3 review pending | Fixed-cube/resource tests and explicitly executed official ASE exact tuples |
| G2 | pending | Actual initialized official MACE graph/E/F parity |
| G3 | pending | Scalar/stencil pairs, E/F/W, skin/bounded regressions |
| G4 | pending | Paired complete-step benefit and fallback overhead |

Source epochs, commands, numerical tolerances, package/model hashes and measured
results will be appended when actually executed. Prepared source contracts do
not count as acceptance evidence.

## G1 — independent periodic image geometry evidence (2026-11-06)

Task 3 adds 30 fixed numeric snapshots and an independent direct Cartesian
`[-16,16]^3` oracle. It imports no production geometry, inverse search or
minimum-image helper. Exact sorted integer-image tuples agree with the actual
borrowed Rust graph for empty/self/multi-image, unreduced and oblique cells,
all eight explicit PBC masks, gauge, deformed and displaced snapshots.
Adjacent-f64/equality boundaries are separate and excluded from ASE. Analytic
triangular numeric certificates, including a search-only +1e-8 cushion, are
in `crates/modena-core/tests/fixtures/periodic_graph/README.md`; their integer
limits are asserted for every fixture. Physical membership retains the
specified strict f64 squared-radius expression without an epsilon.

The first expanded Rust test run passed 17/17 on the already approved graph;
this is an honest test addition, not a manufactured production RED/fix.
Final tests passed 18/18 after adding the independent certificate assertion.
Negative controls detect a deliberately incomplete +/-1 image set and a
changed integer shift without changing production source. Existing Task 2
resource/error/atomic-publication tests remain in the focused suite; previous
Task 1/2 proof and allocation-injection review evidence remains separate.

The explicitly executed official gate passed 1/1 (zero ignored), comparing
all 27 official fixture tuples exactly and independent ASE vectors with
atol=rtol=2e-12. Every official fixture passed fixed-cube tuple and cutoff-margin
preflight before any external helper was invoked. Minimum finite margin was
`0.01922590896378784` Angstrom, above 1e-8; an empty eligible image set has
vacuously infinite margin (serialized as null in the manifest). Maximum
absolute vector error was `1.1102230246251565e-15` Angstrom. The helper consumed
the actual Rust edge lists and used ASE `neighbor_list('ijSD', ...)` as the
independent reference. No packages, interpreter or model weights were installed
or changed.

Commands executed (all shell invocations used `rtk`; output was retained):

```sh
rtk cargo test --locked -p modena-core --test periodic_graph
rtk proxy env MODENA_ASE_PYTHON=/Users/ricardoleal20/Documents/Code/Misc/modena/target/b01-ase-venv/bin/python cargo test --locked -p modena-core --test periodic_graph_official official_ase_images -- --ignored --nocapture
rtk proxy cargo +1.90.0 test --locked -p modena-core --test periodic_graph --test periodic_graph_official
rtk cargo clippy --locked -p modena-core --test periodic_graph --test periodic_graph_official -- -D clippy::all
rtk proxy rustfmt --edition 2021 --check crates/modena-core/tests/support/periodic_oracle.rs crates/modena-core/tests/periodic_graph.rs crates/modena-core/tests/periodic_graph_official.rs
```

The MSRV run passed 18 tests and compiled one explicitly ignored external gate;
that ignored result is not the official acceptance evidence. Focused Clippy
and rustfmt checks passed. The two unchanged Task 1 dead_code warnings
(`Interval::hull`, `scalar_min_image_search_radius`) remain for later consumer
work/final triage. A blanket `-D warnings` result is not claimed.

Source epoch: production graph base `bd95cd2796de16a3ccb4844c26924bd4c3c967d3`;
Task 3 changed only fixtures/test helpers/integration tests and this G1 record.
Approved private geometry SHA256 remains
`b715350a8544a969ed7fb4b84a92a769814a951124fb9f04e4e8e4dc2147342e`.
Fixture bytes SHA256:
`c3f1bf1ec691fb7131b70bbd67605d1eb9e9f8c6b5d6ea6b1421093f41fdf9e6`.
Cargo.lock SHA256:
`f5e02137c735dfff96179f4daa525abf859cca698b0999029f0978a31d1591ac`.

Actual official environment: ASE=3.29.0, NumPy=2.5.4, Python=3.14.5;
macOS-26.6.2-arm64-arm-64bit-Mach-O. Imported neighborlist.py SHA256:
`488c12c35e31d2406db01cae84de6f7fcda0cfd0288b45436d32b786a3eca379`.
Helper SHA256:
`1922883ecf90b69d2c87dc5da042e29067caec679971faaa8f2371d7483d256a`.
Full per-fixture responses, source/compiler/package/platform identity and log
hashes are retained in ignored
`target/b04-evidence/task3-oracle-20261106/manifest.json` (SHA256
`66d60903c9282f3c0622730a9b58faccb78ad6869ca9d36538c214f7636756e5`),
with separate first/final official logs. Logs were created in a new directory;
prior artifacts and immutable campaign ELFs were preserved.

These results cover the declared G1 finite fixtures and the existing supported
numerical/resource contract, subject to independent Task 3 review. They do not
establish official MACE E/F parity, classical stencil correctness or complete-step
performance. G2/G3/G4 remain pending; no performance measurements were run.
Task 1's recorded scalar-fixture Minor and final reviewed B03-base synchronization
remain with the controller and were not expanded into this task.

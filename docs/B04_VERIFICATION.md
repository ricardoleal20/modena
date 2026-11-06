# B04 implementation and verification record

B04 numerical work started in isolated `ricardo/b04-periodic-neighbors`,
initially based on reviewed B03 source `0c54db355b830c49e687e97e2e9aa8502bbaa5d1`.
Its initial core tree was `39b46c7496c6820219056fe38a821378a32633e8`, identical
to reviewed B03 prototype66af0e4. No core/CLI files changed between those bases.
B03 timing/audits finished and its CPU lease was released before this work.
After Tasks1–4 review, unpublished B04 was synchronized onto delivered B03
`12e4ab13efd2aff2d3693b7fc7fb4f50af12668f`. Every B04 source byte and commit
patch matched its pre-sync identity; the mapping and proof are recorded below.
B04 remains unpublished. G1–G3 have independent approval; G4 and final
whole-branch review/delivery remain pending.

Controller preflight confirmed native managed worktree creation, branch prefix
`ricardo/`, source-tree equivalence and the actual Vec3 dot/wrap operation
orders used by the proof. Subsequent G1 implementation evidence is recorded below.

| Gate | Status | Required evidence |
| --- | --- | --- |
| G1 | passed; Task 3 independently approved through ce04b4f | Fixed-cube/resource tests and explicitly executed official ASE exact tuples |
| G2 | passed; Task 4 independently approved at f45be2d, mapped to2aa3960 | Actual initialized official MACE graph/E/F parity |
| G3 | passed; Task 6 independently approved at04343c1 | Certified stencil, scalar/cube pairs, E/F/W, skin/bounded and cap checks |
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

### G1 review fix round 1 — mandatory tuple identity under optimization

Independent Task 3 review found that Python optimization removed the helper's
`assert actual == reference`, allowing count/vector agreement to stand in for
mandatory integer tuples. The controller authorized unconditional mismatch
exceptions in place of the plan snippet. The helper now raises a tuple-specific
ValueError before vector comparison/publication, and enforces ASE/NumPy pins
with unconditional RuntimeError checks. These guards survive PYTHONOPTIMIZE.

The focused ignored test
`official_ase_rejects_wrong_tuple_with_python_optimization` was explicitly run
under PYTHONOPTIMIZE=1. It builds a real Rust graph for two coincident distinct
indices: correct tuples (0,1,0) and (1,0,0), both vectors exactly zero. Replacing
the second tuple with a duplicate first tuple preserves edge count and every
Cartesian vector, so vector allclose cannot detect it. The original optimized
helper accepted this request (behavioral RED, Rust exit 101). After the fix,
the valid request succeeds and the wrong tuple fails specifically at the tuple
guard before vector comparison, publishing no success response (GREEN). The
focused control also passed when executed with Rust 1.90.0.

The actual positive `official_ase_images` gate was explicitly rerun with
PYTHONOPTIMIZE=1 against all 27 official snapshots and the fresh helper source;
1 test passed, zero ignored, 1 intentionally filtered negative-control test.
Minimum finite margin and maximum vector error remain respectively
0.01922590896378784 and 1.1102230246251565e-15 Angstrom. Actual ASE/NumPy/Python,
neighborlist.py hash, fixture bytes and production code remain unchanged.
Focused official-target Clippy and Rust formatting checks passed, with the two
previously disclosed production unused warnings. No broad suite or timing was run.

Current helper SHA256:
`6ea0dd6eb1072c95a5c14b0013009c798baa36ba564c82c0ec26dd9143ca6b53`.
Current official Rust driver SHA256:
`7c1cbaead1be5079591cd372df6f882e74e8e8581e6925be8bd3c8f19a108405`.
The earlier helper hash/manifests above identify the original historical epoch.
New RED/GREEN/MSRV/fmt/lint and all 27 official responses are retained separately
in `target/b04-evidence/task3-fix1-20261106/manifest.json`, SHA256
`7bfab955f16efd34e17e9f50b56d4070bcd1accb69d995d8138a27707a198dcd`.
Its base is `4f0f6f8aa333b908f8b8dfe7b62b7ae8980cbe14`; independent scoped
re-review remains required for this fix epoch. G2/G3/G4 remain pending.

### G1 review fix round 2 — synthetic control margin preflight

The optimized coincident-index control now checks its own fixed-cube cutoff
margin >=1e-8 immediately after changing coordinates, before any ASE invocation.
The focused actual external control passed on stable and Rust 1.90.0 (each:
1 passed, zero ignored, 1 other test intentionally filtered); focused rustfmt
passed. The helper and all 27 positive fixture inputs remain unchanged; their
fix-round-1 evidence above was neither rerun nor rehashed. Current official
Rust test SHA256 is
`25f2424471d1436ea48a057903dd030e630c1ec33027e91105c88a8ac3bf9d30`.
New focused evidence: `target/b04-evidence/task3-fix2-20261106/manifest.json`,
SHA256 `50cedf37f29690ec71499ece66ec77f80217c2d75a7b750a93885d3aaa84a99b`.
Base: `b84c36bedbc39c8b177a845305c629c8eaf5874a`. Scoped re-review remains
required; previously disclosed warnings/Minor and G2/G3/G4 status are unchanged.


### G1 independent review closure

The controller accepted Task 3 and both fixes through
`ce04b4f659af9654adef0f77ccd56d92d0d053e9`. This closes the review-pending
qualifiers of the historical G1 epochs above. Task 4 leaves all G1 fixture,
oracle, ASE helper and production geometry bytes unchanged.

## G2 — actual initialized official MACE software parity (2026-11-06)

The named external gate actually passed on stable Rust 1.97.1 and the Rust
1.90.0 floor: each final run executed one test, zero ignored, with the two ASE
tests intentionally filtered. Independent Task 4 review approved the complete
source/environment/gate record; see the closure below.
The helper uses official MACE code, official Configuration/AtomicData/Batch,
the same initialized model object, and actual serialized native PeriodicGraph
edges. Fresh position and node tensors enter every gradient-enabled forward.
Only edge_index, exact-integral unit_shifts and Cartesian shifts change in the
alternate batch; all other fields are checked for exact identity.

The original required NumPy 2.2.6 / matscipy 1.1.1 combination failed resolution:
[matscipy 1.1.1 metadata](https://pypi.org/pypi/matscipy/1.1.1/json) requires
NumPy >=1.16.0,<2.0.0. Controller contract revision
`0a01e55ea9f65ad609600d5a29df955d4c36a866` changed only NumPy to 1.26.4
before model evaluation. Exact installed direct versions are torch 2.6.0,
mace-torch 0.3.13, e3nn 0.4.4, numpy 1.26.4, ase 3.26.0 and matscipy 1.1.1.
Installed Python remains 3.12.8, macOS-26.6.2-arm64-arm-64bit. A pyenv shim
initially failed discovery; PYENV_VERSION=3.12.8 selected the already installed
exact runtime. No alternate Python/runtime/framework was substituted.

`bench/requirements-b04-mace-py312.lock` freezes all 44 runtime packages with
hashes. `uv pip sync --require-hashes` installed only into target/b04-mace-venv;
`uv pip check` passes. Installed RECORD hashes verified 21,938 files. Retained
uv archive receipts identify actual wheels; installed WHEEL tags disambiguate
multiple cached PyYAML archives. python-hostlist 2.3.0 has only a source release:
its hash-locked sdist built a retained wheel using setuptools 84.2.1 (recorded
in WHEEL). The source and built-wheel hashes are respectively
`e1a0b18e525a5fca573cb9862799f11b3f2bd3ba7aec70c4ecd8b95341bb71ea` and
`5f81f639eb73a350c84143877eb46b534387028751f2afca111d75d5e30b2c1c`.
The B01 environment and global Python packages were not modified.

Torch 2.6's weights-only default initially rejected builtins.slice in e3nn's
packaged constants. The controller approved a narrowly scoped
`torch.serialization.safe_globals([slice])` only around importing e3nn.o3.
The helper checks the packaged constants SHA256 before importing, verifies the
prior allowlist is restored, and rejects the unsafe environment override.
No weights_only=False, global loader relaxation or pretrained checkpoint is used.
Constants SHA256: `279fb1f7844b2b85df092550bf7854a70fe3b427b312672f6a1e1fef0a786876`.
Independent imports and the exact constructor succeeded before meaningful RED.

The model is the exact design constructor, seed 20261106, CPU float64, one thread,
r_max 3.1, two interactions, elements [6,8], hidden 8x0e+8x1o, correlation 2,
and 1,968 initialized parameters. SHA256 over sorted state_dict key/dtype/shape
headers and contiguous CPU tensor bytes (length-delimited) is
`1d1a81cf8cfaa7885469757c26b82d49944cd188034ba7f5c481ad020974066e`.
MACE release v0.3.13 resolves to official source commit
`b5faaa076c49778fc17493edfecebcabeb960155`; installed models.py, utils.py,
atomic_data.py and neighborhood.py are byte-identical to that immutable source.
The actual MACE wheel SHA256 is
`e0aeef5fdb493f94436e4a806f89e051240d79202bccd6bacf5c6007d9adf75a`.

All four frozen full-PBC fixtures (mace3/gauge/deform/displaced) have 76 edges.
Exact tuple parity and fixed-cube cutoff margins pass before evaluation.
Official original order, canonical official order and native graph forwards
agree within atol=rtol=1e-10. Maximum energy error is 5.551115123125783e-17;
maximum force-component error is 1.717376241217039e-16. Same-process repeated
forwards are exactly equal; final stable/MSRV model digests and all E/F arrays
also agree exactly. Gauge-translated E/F agreement passes.

For the displaced fixture, Rust rebuilds 18 actual native graphs for all nine
Cartesian central differences at h=1e-5. Their neighbor tuples stay unchanged,
all cutoff margins exceed 1e-8, and each displaced input also passes official
tuple and E/F comparison. Maximum autograd-versus-central-difference force
error is 1.693140072589472e-11, below atol=rtol=1e-6.

The deterministic controls drop the first nonzero image edge, drop all periodic
self images, and change the first nonzero component of the first nonzero edge
away from zero by one integer. Across the fixed suite their maximum E/F effects
are 0.0017405161288699955, 0.01902269751579285 and 0.0017405161288699816,
respectively, all above 1e-10. The original +1-to-x control accidentally created
a zero-self edge and failed with nonfinite model forces. Its preflight/log is
retained; the controller approved the stated deterministic repair before the
next evaluation. The failed attempt is not acceptance evidence, and no seed,
model, fixture, tolerance or mutation-sensitivity threshold was tuned.

All mandatory checks use unconditional exceptions or always-on checking
functions. Every model subprocess runs Python -O. Positive runs additionally
exercise non-graph and zero-sensitivity guard failures; external negative
controls reject duplicated tuples, fractional shifts and wrong schema without
publishing success. Separate optimized probes reject wrong dependency pins and
unsafe loader overrides before importing frameworks. The missing-interpreter
control errors clearly, as required. Numerical/source/model/input manifests
are captured before evaluation; failed attempts and successful epochs are kept
separately under the new ignored task4-mace-20261106 directory.

Commands (full outputs and setup/failed-attempt commands are in the task report):

```sh
rtk proxy env PYENV_VERSION=3.12.8 uv venv --python 3.12.8 --no-python-downloads target/b04-mace-venv
rtk proxy uv pip compile --python target/b04-mace-venv/bin/python --generate-hashes bench/requirements-b04-mace.in -o bench/requirements-b04-mace-py312.lock
rtk proxy uv pip sync --python target/b04-mace-venv/bin/python --require-hashes bench/requirements-b04-mace-py312.lock
rtk proxy env MODENA_MACE_PYTHON=target/b04-mace-venv/bin/python cargo test --locked -p modena-core --test periodic_graph_official official_mace_initialized_graph_energy_forces -- --ignored --nocapture
rtk proxy env MODENA_MACE_PYTHON=target/b04-mace-venv/bin/python cargo +1.90.0 test --locked -p modena-core --test periodic_graph_official official_mace_initialized_graph_energy_forces -- --ignored --nocapture
rtk cargo clippy --locked -p modena-core --test periodic_graph_official -- -D clippy::all
rtk proxy rustfmt --edition 2021 --check crates/modena-core/tests/periodic_graph_official.rs
rtk proxy uv pip check --python target/b04-mace-venv/bin/python
```

Evidence runs additionally set MODENA_MACE_EVIDENCE to the absolute ignored
artifact directory to retain preflight request/model/package manifests. The
first meaningful RED is log 13 (missing helper, after successful constructor
probe); log 12 retains the earlier broken-pipe diagnostic failure. Final
stable/MSRV logs are 23/21; focused final Clippy and rustfmt pass. The two
previously disclosed production unused warnings remain deferred. No broad
suite, production changes or performance campaign was run.

Final evidence manifest:
`target/b04-evidence/task4-mace-20261106/manifest.json`, SHA256 `8e6a85963863491d00f85386e6d330d9f8c2d5f3bfe4451259425b8d750a9ccf`.
Lock SHA256: `a94b437c2d9faf76fd0c404aa7a6b9702d6402d3b6cf098b88ba145c5125dcef`.
Helper SHA256: `4ddafd56ed5bee74fe81170ec3a16469eacdb798e9faf8ffded9028bbe0d3cff`.
Fixture bytes remain `c3f1bf1ec691fb7131b70bbd67605d1eb9e9f8c6b5d6ea6b1421093f41fdf9e6`.

This is software parity for this initialized MACE instance, snapshots, locked
packages and arithmetic only. It establishes no trained chemical accuracy,
FAIR-Chem task/checkpoint equivalence or B10 acceptance. G3/G4 remain pending.

## Task 4 independent review closure

Task 4 approved at original source `f45be2d363e1050173962b9ba8c21bf4e0082453`,
with no Critical or Important finding. The reviewer verified all121 artifact
hashes,44 lock entries and selected package artifacts, official MACE source
identity, official autograd implementation and identical stable/MSRV numerical
results. No model evaluation or suite was rerun during review.

The pinned dependencies also emit upstream TorchScript annotation warnings
during constructor compatibility/first-gate runs. These are retained diagnostic
noise, in addition to the two existing unused Rust helpers. Their disposition
is explicit disclosure without dependency changes or broad warning suppression;
final-review triage remains. G2 software parity does not close G3/G4 or B10.

## Unpublished base synchronization, 2026-11-06

Before synchronization, remote B04 was verified absent and the worktree clean.
Backup `ricardo/b04-pre-b03-sync-20261106` retains the complete old history.
Rebase onto B03 `12e4ab1` completed without conflicts, force-push or a merge.
At mapped head `2aa39602de5d8fb10cac20cbcadca7fbb746574d`, all15 B04 changed
files and every one of eight individual commit patches are byte-identical.
The complete B04 delta SHA256 before/after is
`06a15ec8cbc956dd0f4a31389d8add6981a93894c21d5ba152a5b986258ab076`;
core tree remains `743ce83c62ef42e4021fa6670beadc394aa5468e`.

| Historical commit | Mapped commit |
| --- | --- |
| `72be6d56571d2877d54e3fbf927cc8eeeb9d1e94` | `c679437b344e4eb71ac2f3406552339553d07a34` |
| `bac0a55a91c5ca5a161a89c91adb3747e19db392` | `97454a34efbe325f41f5fee1bb953f3044a13b4e` |
| `bd95cd2796de16a3ccb4844c26924bd4c3c967d3` | `ab3df0bcb03270363a49a37d09690da08b96e26e` |
| `4f0f6f8aa333b908f8b8dfe7b62b7ae8980cbe14` | `d31296a5989c987a3c1708975d027be526767174` |
| `b84c36bedbc39c8b177a845305c629c8eaf5874a` | `3cc4e2a560d707823a642d0fcb3f6d6cd82f7460` |
| `ce04b4f659af9654adef0f77ccd56d92d0d053e9` | `a91c20d6d86307d5eacb123ed905d8bb0dc01930` |
| `0a01e55ea9f65ad609600d5a29df955d4c36a866` | `d7148f42721136eafba03f6f096b54eb755378bd` |
| `f45be2d363e1050173962b9ba8c21bf4e0082453` | `2aa39602de5d8fb10cac20cbcadca7fbb746574d` |

The full receipt and before/after patches are retained under
`target/b04-evidence/base-sync-20261106/`. Earlier commands, failures, numerical
results and source hashes retain their actual historical epochs; they are not
relabeled as new evaluations after rebase. Source identity is the evidence
for carrying those reviews forward. Later delivery metadata changes only docs.


## G3 — certified fractional triclinic stencil (2026-11-06)

Task 6 implements forced FractionalStencil execution after reviewed Task 5 base
`76143e76d6d190f4f2bbd74810fa1651ebf5c824`. Auto still executes Scalar without
running stencil preparation. G3 functional/numerical checks and independent
Task 6 review passed. G4 is pending, with no timing or default promotion.

The existing Task 1 scalar radius proof is consumed, not replaced. Its
subtraction/dot error and inner-index calculations were extracted into private
helpers for direct numeric tests without changing their arithmetic. The scalar
operation audit still matches `SimBox::min_image` and Vec3: one coordinate
subtraction, five dot operations, exact `t-round(t)` below 2^52, coefficient
subtraction, and five Cartesian reconstruction operations. The norm envelope
covers three products/two sums and gradual underflow; J must be <=4096. The
stored list cutoff squared is converted to an upward-rounded search radius,
so its reconstructed square covers the original threshold. Final pruning uses
the original `min_image(...).norm2() < list_cutoff2`, unchanged.

For each atom, interval evaluation of A*r plus the certified inverse residual
`eta*||r||inf` encloses the true U*r. Equal endpoint floors certify periodic
range reduction with an exactly representable integer. Outward subtraction and
multiplication by m must place the entire closed interval in one bin; seam or
face ambiguity falls back. For a scalar-accepted image, Q encloses its real
Cartesian displacement. Hence each fractional image difference is bounded by
w=Q*||U_row||, and bin differences necessarily satisfy
`dist(q/m, integers) <= w + 1/m`. Both the offset range and filter are outward;
a rejection requires a lower distance bound to exceed the upper threshold.
The reciprocal's lower endpoint determines a possibly coarser m at an integer
boundary, which preserves this proof.

For m>1, m*w<=1; one outward endpoint ulp gives at most seven enumerated
pre-wrap offsets per axis. For w>=1, `m=max(1,floor(1/w))=1` and all integer
offsets are the same residue, collapsed before any range enumeration. Modulo
residues are deduplicated before their Cartesian product, then translated to
unique target bins. Thus the fixed 4096-offset cap is defensive under this
grid rule (at most 7^3 offsets), not a physically triggered failure claimed by
the fixtures. The fixed cap was not retuned.

A fallibly allocated, sorted (bin, atom) list inserts each atom once. Visiting
only source-bin<=target-bin counts each unordered occupied-bin pair once;
same-bin counts use C(n,2), cross-bin counts use n*m. Checked u128 arithmetic
implements the specified candidate_estimate < 0.8*all_pairs rule before pair
emission. The grid, offset, N*offset, i32-index and workspace caps precede heap
requests. Workspace includes the requested atom/offset buffers and fixed
plan/certificate/axis storage, excludes the existing output buffer, and is
checked again using actual Vec capacities. There is no dense global pair set,
unbounded candidate staging, or accepted-pair deduplication. Output is sorted
lexicographically. The shared bounded emitter propagates max_pairs failures;
these cannot become successful scalar retries.

Actual optimized fixtures all compare exact ordered pairs with retained Scalar
and independent fixed-cube membership. The original 30 fixed snapshots also
compare against Scalar; their full-PBC subset compares to the cube. Partial-PBC
fixture masks remain graph metadata, since classical neighbors are full PBC.
Empty builds retain their existing no-diagnostics behavior.

| Actual selected shape (8 atoms, 4 accepted pairs) | Predicate checks | Requested workspace |
| --- | ---: | ---: |
| 32-cell, tilts -8/0/8/24 | 4 vs 28 Scalar | 1216 bytes |
| Unreduced tilts -128/+128, individual atom lattice translations | 12 vs 28 Scalar | 1216 bytes |
| First-axis dimension 1 (length 4) | 4 vs 28 Scalar | 784 bytes |
| First-axis dimension 2 (length 6) | 4 vs 28 Scalar | 1000 bytes |
| First-axis dimension 3 (length 8) | 4 vs 28 Scalar | 1216 bytes |

Real fallback diagnostics cover uncertified geometry, fractional boundary,
image domain, grid cap, all bins and candidate cost, with exact Scalar parity.
WorkCap and MemoryCap are exercised by count-only preflight at the fixed
production limits, without allocating huge atom vectors. StencilCap is likewise
a defensive count test; overflow cases are checked directly. These are separate
from admitted-geometry diagnostics. Memory/allocation uncertainty always falls
back before candidate emission. No allocator-failure injection was added.

The meaningful RED run executed 13 tests: eight existing tests passed and five
new tests failed specifically because forced requests still selected Scalar
or lacked the boundary reason. Final targeted tests cover exact ordered pairs,
zero/independent/amplified pair policies, additive force/virial seeds, E/F/W
at atol=rtol=1e-10, list-radius skin membership, the half-skin threshold, a
periodic seam crossing, full-box cache invalidation, and actual private bounded
stencil overflow/recovery. Numeric tests directly cover t_error, the 2^52
wrapping domain, J<=4096, certified closed bins, and the actual scalar-selected
real-image counterpart on every fixed cube fixture. The earlier Task 1 Q test
now derives all tested displacement components from its supplied snapshot.

Verification commands (all run through RTK with complete logs):

```sh
rtk cargo test --locked -p modena-core --features neighbor-diagnostics --test triclinic_neighbors --test force_hardening --test snapshot --test pair_scales
rtk cargo test --locked -p modena-core --features neighbor-diagnostics --lib periodic_geometry::tests
rtk cargo test --locked -p modena-core --features neighbor-diagnostics --lib neighbor::
rtk proxy cargo +1.90.0 test --locked -p modena-core --features neighbor-diagnostics --test triclinic_neighbors --test force_hardening --test snapshot --test pair_scales
rtk proxy cargo +1.90.0 test --locked -p modena-core --features neighbor-diagnostics --lib periodic_geometry::tests
rtk proxy cargo +1.90.0 test --locked -p modena-core --features neighbor-diagnostics --lib neighbor::triclinic::tests
rtk proxy cargo +1.90.0 test --locked -p modena-core --features neighbor-diagnostics --lib bounded_stencil_overflow_never_succeeds
rtk cargo test --locked -p modena-core --no-default-features --test force_hardening --test snapshot --test pair_scales
rtk cargo test --locked -p modena-core --no-default-features --lib bounded_stencil_overflow_never_succeeds
rtk proxy cargo +1.90.0 check --locked -p modena-core --no-default-features
rtk cargo clippy --locked -p modena-core --features neighbor-diagnostics --lib --test triclinic_neighbors --test force_hardening --test snapshot --test pair_scales -- -D warnings
rtk cargo clippy --locked -p modena-core --no-default-features --lib -- -D warnings
rtk cargo fmt --all -- --check
```

Final feature integration suite: 36 pass (16+11+3+6), including Rust 1.90.
Geometry: 18 pass; neighbor unit scope: 18 pass. Focused MSRV cap/bounded units:
3+1 pass. Default integration: 20 pass; default actual bounded stencil: 1 pass.
Feature and default Clippy are warning-free with `-D warnings`; previously
staged Rust warnings are resolved by real helper/variant use and feature-gating
only the diagnostic report fields. Existing TorchScript warnings belong to the
unchanged G2 environment and were not rerun or suppressed.

Fresh RED/GREEN/failed-attempt/MSRV/default/lint logs, source identities and full
commands are retained under `target/b04-evidence/task6-stencil-20261106` and in
`.superpowers/sdd/2026-11-06-b04-periodic-neighbors/task-6-report.md`. The first
implementation compile needed an explicit generic error type; a later private
test was initially placed outside its module. Both failed logs are preserved.
No graph/model/helper semantics, System/Vec3 arithmetic, exclusions, frozen
ELFs or historical evidence changed. These results do not establish G4,
trained chemical accuracy or B10 acceptance.

### G3 independent review closure

Task6 independently approved `04343c16f9f2149fa7f5472b5c932902f2568fde`,
with no Critical, Important or Minor finding. The reviewer verified the finite-
precision/search-only bounds, certified bins/offset coverage and uniqueness,
actual selected-path fixtures, exact unordered occupancy accounting, allocation
caps, bounded failure propagation, pair ordering and original final predicate.
Ten source and21 log hashes match the retained manifest; actual stable/MSRV,
feature/default, numeric/neighbor and warning-clean lint evidence was checked
without rerunning it. Task1's Q-fixture mismatch and staged native warnings are
closed. Historical G1/G2 approvals remain scoped to their actual execution
epochs. No G4 timing/default promotion or chemical efficacy is established.

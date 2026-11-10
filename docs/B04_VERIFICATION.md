# B04 implementation and verification record

B04 numerical work started in isolated `ricardo/b04-periodic-neighbors`,
initially based on reviewed B03 source `0c54db355b830c49e687e97e2e9aa8502bbaa5d1`.
Its initial core tree was `39b46c7496c6820219056fe38a821378a32633e8`, identical
to reviewed B03 prototype66af0e4. No core/CLI files changed between those bases.
B03 timing/audits finished and its CPU lease was released before this work.
After Tasks1–4 review, unpublished B04 was synchronized onto delivered B03
`12e4ab13efd2aff2d3693b7fc7fb4f50af12668f`. Every B04 source byte and commit
patch matched its pre-sync identity; the mapping and proof are recorded below.
B04 remains unpublished. G1–G3 have independent approval. G4 is not passed
(incomplete frozen campaign); Auto remains Scalar. Final whole-branch
review/delivery remain pending.

Controller preflight confirmed native managed worktree creation, branch prefix
`ricardo/`, source-tree equivalence and the actual Vec3 dot/wrap operation
orders used by the proof. Subsequent G1 implementation evidence is recorded below.

| Gate | Status | Required evidence |
| --- | --- | --- |
| G1 | passed; Task 3 independently approved through ce04b4f | Fixed-cube/resource tests and explicitly executed official ASE exact tuples |
| G2 | passed; Task 4 independently approved at f45be2d, mapped to2aa3960 | Actual initialized official MACE graph/E/F parity |
| G3 | passed; Task 6 independently approved at04343c1 | Certified stencil, scalar/cube pairs, E/F/W, skin/bounded and cap checks |
| G4 | not passed / incomplete; Auto=Scalar | Natural-rebuild/overhead requirements unresolved; one 600s timeout, 35 later unstarted cases |

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

## G4 — Task 7 harness preparation, 2026-11-06

Status: prepared and functionally checked; **G4 is pending**. No campaign,
performance pilot, accepted throughput result or selector promotion was run.
Task 8 requires the controller's exclusive team CPU/write reservation, with B05
paused first. The runner requires `MODENA_B04_CPU_LEASE=exclusive`; that marker
records the controller's reservation and does not acquire a machine lock.
Auto remains Scalar. B03 campaign files and production core are unchanged.

The opt-in bench feature `neighbor-diagnostics` forwards the core feature;
the named `triclinic_neighbors` binary requires it. Runtime needs this source
checkout, Git, Python 3 and `shasum`, with no RTK dependency. Public commands:

```sh
rtk proxy python3 -m unittest discover -s bench/tests -p test_triclinic_neighbors.py
rtk cargo test --locked -p modena-bench --features neighbor-diagnostics --bin triclinic_neighbors
rtk proxy cargo +1.90.0 test --locked -p modena-bench --features neighbor-diagnostics --bin triclinic_neighbors
rtk cargo build --locked --release -p modena-bench --features neighbor-diagnostics --bin triclinic_neighbors
rtk proxy python3 bench/triclinic_neighbors.py --binary PATH --output-dir PATH --generate-only
rtk proxy python3 bench/triclinic_neighbors.py --binary PATH --output-dir PATH --run
rtk proxy python3 bench/triclinic_neighbors.py --output-dir PATH --analyze-only
```

Compilation must additionally set `B04_BUILD_REV`, `B04_BUILD_SOURCE_SHA256`
and `B04_RUSTC_VERSION` in the compiler environment, calculated from the actual
checkout's `source_identity()` and actual `rustc --version` before building.
The exact executed command/environment wrapper is retained with the evidence.
Source SHA256 covers sorted tracked Rust files, Cargo manifests/lock and the
new Rust/Python harness, with NUL-delimited paths and contents. Records bind
embedded build revision/source/compiler, runtime revision/source, actual
executable/input SHA256, feature and profile. Before kernels, a mismatched
runtime source/revision rejects a stale executable. Task 7's precommit
dirty-source build is a preparation epoch; Task 8 must rebuild at its clean,
reachable frozen HEAD and retain a fresh immutable executable/manifest.
Committing deliberately invalidates the precommit binary's runtime revision.

The manifest fixes N=256/2048/8192, density .8, the four exact design matrices
(including elongated's moderate off-diagonals), skins 0/.3, threads 1/4, both
modes, LJ sigma=epsilon=mass=1, cutoff 2.5, no tail correction, and the actual
`LennardJones::new` shifted-force treatment. Grid coordinates, u64-safe hash
products/modulo 65536, seed 20261106, velocity coefficient .05 and atom-order
mean removal follow the design. Resulting f64 positions/velocities are
serialized once and loaded identically by both strategies. Separate controls
use an eight-atom grid-cap cell, an origin-on-boundary input and a 256-atom
b_x=.2L/.21L full-basis change. The grid-cap control records its actual density.
There are 15 snapshots and 116 planned paired cases, with basis changes in
rebuild mode only. Generation is not numerical acceptance: execution checks
all initial pair distances finite/positive and force/list-cutoff margin >1e-8
before timing; failed fixtures stay invalid without input changes. Both
basis-change states receive preflight and scalar/candidate physical comparison.

Each process owns independent scalar/candidate System, NeighborList,
PairForceModel and VelocityVerlet instances. Actual sorted half pairs and all
E/F/W, position and velocity components are compared before timing, after
each paired block and finally. Checks use atol=rtol=1e-10; recorded maxima are
normalized errors `abs(a-b)/(1+max(abs(a),abs(b)))`. Validation never mutates
live trajectories/counters. Five warmup and 20 measured paired blocks each
advance ten steps at dt=.0005; AB/BA alternates by phase-local sample index.
Both replicas continue from their preceding block. Step mode times complete
VV/neighbor/force work plus per-step diagnostic accounting. Rebuild mode times
whole forced NeighborList rebuilds; VV advances outside those timers to keep
continuous trajectories. Those force/advance costs are explicitly excluded
from the rebuild metric. Basis-change cases alternate the two full cells
before forced builds and refresh forces before outside-timer VV advance.

Actual rebuild deltas gate diagnostic accounting: cached reports contribute
no repeated candidate checks/workspace. Records retain the last timed build's
selection/fallback/pair count, and counts of every actual selected/fallback
path. Warmups and measured records remain in chronological order. Errors
publish invalid JSON/nonzero exit. Child failures, missing/zero durations,
bool-as-int metadata, mismatched identities, absent counterparts, wrong order,
nonfinite physics and incomplete rosters cannot become valid analysis.
The case cap is 600 seconds including preflight; the parent also enforces a
timeout and 7200-second campaign cap. Aborts retain partial files/child receipts
as incomplete evidence. Exclusive file creation prevents historical overwrite.

Analysis uses matched candidate/scalar ratios, their median and 10,000
bootstrap resamples seeded 20261106, with bounds at indices 249/9749. Fixed
preparation-time screens mark uncertainty if even/odd order median ratios
differ by >10%, first/last-half medians drift by >10%, or bootstrap width
exceeds 10% of the median. These screens do not identify a thermal cause.
Natural skin .3 step cases need >=2 measured rebuilds in each replica; the
protocol is never extended after observing too few. Rebuild gain is reported
separately and cannot satisfy full-step promotion. Moderate N>=2048 at every
skin/thread regime needs median<=.95 and CI upper<1; small/fallback cases need
actual scalar selection and <=5% median overhead. Noisy fallback overhead
stays inconclusive. The offline report leaves G4 pending for gate review.

Memory scope is explicit because core neighbor/model capacities are private.
`requested_bytes` is maximum reported temporary neighbor workspace across
actual timed new builds, excluding pairs/model/process storage. Scalar/cache
zero means zero new neighbor temporary workspace, never zero process memory.
`allocation_count` counts actual process-wide Rust alloc/alloc_zeroed/realloc
calls during timers, including failed requests. `retained_bytes` is live
Rust-requested heap immediately after a timed block or its final timed rebuild:
it includes BOTH replicas, model scratch, pool, counters and prior output
records. It is not per-strategy memory, allocator usable size or OS RSS.
Live bytes change only on successful allocation. Result-vector storage and
reason maps are reserved outside timers. Pool/model setup is untimed; deferred
worker allocations count only if they occur during a timer, with no model
attribution. Both controls use identical counting allocator/relaxed atomic
instrumentation. Its overhead is unquantified; instrumented cost alone does
not establish uninstrumented production throughput. Any default promotion
requires final review of this limitation as well as the full-step criteria.

Preparation checks: 16 Python parser/statistical tests; three actual Rust tests
on stable and Rust 1.90; named release build; focused Clippy `-D warnings` and
rustfmt. Missing-counterpart RED/GREEN and a fixed nonconstant seeded-bootstrap
mutation control fail for the intended reasons. Synthetic durations test
arithmetic only. Five eight-atom G3-derived correctness smokes cover actual
forced-stencil step/rebuild, one/four threads, both skins, paired Auto and
alternating-basis rebuilds.
They use the fixed 5/20/10 protocol; every duration is excluded from G4.
Pairs and E/F/W/state agree exactly. Auto selects Scalar. The natural-skin
smoke sees zero measured rebuilds and is inconclusive. A bool-as-N snapshot
publishes retained invalid JSON/nonzero exit. The large ladder is generated
only and has not been numerically evaluated. Running without the CPU lease
marker fails before launching a child. Analysis of the unrun ladder retains
116 unavailable cases as invalid and exits nonzero.

Exact raw commands, source/binary/input hashes, failed attempts and immutable
smoke epochs are retained under `target/b04-evidence/task7-harness-20261106/`.
The detailed handoff is
`.superpowers/sdd/2026-11-06-b04-periodic-neighbors/task-7-report.md`.
These checks establish preparation only, with no model/tool-engine parity,
chemical efficacy, G4 pass or B10 claim.

### Task 7 I1 correction — explicit fallback coverage

Review I1 identified that the original small boundary cell returned AllBins
before reaching fractional bin certification. That original input and all
prior smoke/prepared epochs remain preserved and are superseded for boundary
coverage. Before new numerical evaluation, the controller fixed the separate
boundary recipe below and recorded exact coordinates/velocities in the Task7
handoff. The density-.8 ladder, grid-cap input, basis-change input, physics,
velocity formula and fixed timing protocol are unchanged; fourteen snapshot
hashes agree exactly with the previous prepared collection.

Boundary control: N=8, basis columns (32,0,0), (6.4,32,0), (3.2,4.8,32), actual
density 0.000244140625. Atom0 is isolated at (0,0,0) with v=(0,0,0). Atoms1..7
are the seven nonzero-bit corners of a Cartesian cube with spacing 1.4 anchored
at H*(.5,.5,.5)=(20.8,18.4,16); x/y/z correspond to bits0/1/2. Initial cluster
velocities use the same .05 formula for i=1..7, with their own seven-atom
component means removed in atom order. Exact f64 coordinates and velocities
are retained in the new snapshot/recipe, with all eight explicit vectors in
the pre-evaluation handoff. The large cell avoids the earlier all-bins cap,
and the isolated atom remains exactly at the ambiguous fractional boundary
under continuous VV. No state reset or velocity modification occurs in timers.

The frozen manifest now declares boundary=>fractional_boundary and
grid_cap=>grid_cap. For forced-fractional controls, binary and analyzer require
actual Scalar selection and the expected initial reason. Every new build in
each block must have that exact reason, verified from selected/fallback counts
as well as the last timed report; a wrong reason cannot hide behind a correct
last report. Initial diagnostics are serialized as actual requested/selected/
fallback triples for both replicas. The analyzer independently validates them.
Zero measured relevant builds gives `inconclusive_control_no_measured_builds`,
with fallback coverage and overhead acceptance false. Auto separately checks
scalar-dispatch overhead and reports `auto_scalar_dispatch_only`, never forced
fractional safety-path coverage and never a fabricated fallback reason.

Fix verification: 23 Python tests pass, including meaningful wrong/missing
reason RED/GREEN, one-hidden-wrong-build rejection, independently wrong initial
reason, no-measured-build inconclusiveness and Auto dispatch separation. Five
Rust tests pass on stable and Rust1.90; real fixtures evolve continuously for
250 steps at both skins with exact pairs/E/F/W/state parity. Actual diagnostics
match the intended paths at initialization and every new build. Atom0's
position, velocity and force remain exactly zero without reset. An actual
AllBins substitute is rejected by the native initial-control guard.

A new immutable release/source epoch and eight eight-atom CLI correctness
smokes are retained under
`target/b04-evidence/task7-harness-20261106/fix1-boundary-20261106/`:
boundary and grid_cap each run fractional step/skin0/thread1, fractional
rebuild/skin.3/thread4, fractional step/skin.3/thread4 and Auto step/skin0/thread4.
Forced step/skin0 and rebuild/skin.3 each observe200 measured builds per
replica, every candidate build with the expected reason. Natural skin.3 sees
zero measured builds, remains inconclusive and is not extended. Auto selects
Scalar with no fallback and makes no fractional coverage claim. All five
physical/state error maxima are zero. A CLI negative using unchanged old
boundary bytes retains actual AllBins diagnostics and fails before kernels.
The new fifteen-snapshot collection is generated only; the large ladder is
not numerically evaluated. Focused Clippy/rustfmt and release compilation pass.
All smoke durations are excluded from G4 performance acceptance. No CPU lease,
selector promotion or performance run is claimed by this correction.


## G4 — frozen campaign decision, 2026-11-06

**Not passed / incomplete; Auto remains Scalar.** The first accepted campaign
ran once at clean `f56f28854865e8068402eebf7a0548f60da8ecfd` with the approved
Task7/I1 controls. No selector or other production source changed. The required
moderate N>=2048 skin=.3 step regimes each observed zero measured rebuilds per
replica, below the fixed requirement of two. They remain inconclusive; their
trajectories were not extended. No narrower post-hoc promotion is made.

The driver finished with exit 1 after 6344.874002933502 seconds. Eighty cases
completed; `unreduced-8192-rebuild-skin0-t1` hit the fixed 600-second parent
timeout and retained an empty result file plus its child failure receipt.
The fail-closed driver then stopped. Thirty-five later cases were not started;
all 36 invalid/incomplete rows remain in the 116-case analysis. The 7200-second
total cap was not reached. The boundary/grid-cap and basis-change campaign
controls were among the unstarted cases, so their performance/reason coverage
is unavailable. Earlier Task7 correctness smokes remain a separate epoch.

| Analysis status | Cases |
| --- | --- |
| rebuild_gain_only | 26 |
| inconclusive_order_or_noise | 17 |
| no_gain | 12 |
| inconclusive_rebuilds | 12 |
| eligible_gain | 13 |
| invalid | 36 |

All 80 completed cases (4,000 raw warmup/measured records) passed exact pair
comparison and the specified E/F/W/position/velocity component checks; all five
maximum normalized errors are zero. This claim excludes the timed-out and
unstarted cases. Initial overlap/cutoff-margin checks, source/binary/input
identity, AB/BA order, cadence and counter guards passed for every accepted
case. The controller independently audited the 4,000 records before freezing
the data and releasing the team CPU lease.

### Required moderate full-step regimes

Ratios are paired forced-request/scalar block medians from the instrumented
binary, not ratios of independently pooled timings. The fixed 10,000-resample
bootstrap uses seed 20261106 and order statistics 249/9749. Raw full precision
governs eligibility; values below are rounded for display.

| Full-step case | Paired median ratio | 95% bootstrap CI | Measured rebuilds S/C | Outcome |
| --- | --- | --- | --- | --- |
| moderate-2048-step-skin0-t1 | 0.24017377 | [0.23609382, 0.24586326] | 200/200 | eligible_gain |
| moderate-2048-step-skin0-t4 | 0.24147685 | [0.23692591, 0.24456883] | 200/200 | eligible_gain |
| moderate-2048-step-skin0.3-t1 | 0.99344971 | [0.99265936, 0.99429888] | 0/0 | inconclusive_rebuilds |
| moderate-2048-step-skin0.3-t4 | 0.99978947 | [0.99953754, 1.0001539] | 0/0 | inconclusive_rebuilds |
| moderate-8192-step-skin0-t1 | 0.072528684 | [0.072351043, 0.072616933] | 200/200 | eligible_gain |
| moderate-8192-step-skin0-t4 | 0.065565116 | [0.065368546, 0.065785975] | 200/200 | eligible_gain |
| moderate-8192-step-skin0.3-t1 | 1.0011997 | [0.9978052, 1.0020715] | 0/0 | inconclusive_rebuilds |
| moderate-8192-step-skin0.3-t4 | 0.99869308 | [0.99826903, 0.99999206] | 0/0 | inconclusive_rebuilds |

The zero-skin rows satisfy the numerical gain criterion for this executable
and host. That does not close G4: all four required moderate natural-skin rows
are inconclusive. Seventeen other rows fail the fixed order/drift/CI-width
noise screens and remain inconclusive; no thermal cause is assigned. Forced
rebuild gains are labeled separately and cannot establish full-step acceptance.

The small/fallback requirement is also unresolved/unsatisfied. Only six of the
sixteen completed N=256 step rows satisfy the scalar-selection/<=5% overhead
rule. The other ten include noise, absent rebuilds and actual fractional
selection in elongated cases; favorable timing does not satisfy required
scalar selection. The declared fallback controls were not reached. No losing,
noisy, missing or unexpected-selection row was filtered out. Complete rebuild
and step ratios, intervals, absolute durations and every status remain in the
full analysis and the detailed Task8 report.

### Counter scopes and measured platform

Across the 3,200 accepted measured records, candidate checks range 0..335,503,360,
temporary neighbor requested bytes 0..132,160, process-live requested heap
562,525..22,066,359 bytes and timed allocation attempts 0..20 per block.
New-build selected counts across both roles are Scalar 18,240 and Fractional 7,380;
actual candidate fallback counts are AllBins 3,618 and FractionalBoundary 1,812.
These incidental ladder fallbacks do not replace the unrun declared controls.

`allocation_count` counts process-wide Rust alloc/alloc_zeroed/realloc attempts
inside timers, including failures. `retained_bytes` includes BOTH replicas,
pool, scratch, benchmark maps and prior records; it is neither per-strategy
ownership nor usable size or OS RSS. `requested_bytes` covers temporary neighbor
workspace only, and zero on cached/scalar records is not zero process memory.
Both arms use identical relaxed-atomic counting instrumentation. Its overhead
is unquantified; uninstrumented production throughput is not established.

Host: Mac15,6, Apple M3 Pro, 11 physical/logical cores, 19,327,352,832 bytes RAM;
macOS 26.6.2 arm64 (Darwin 25.6.0). The native executable is Mach-O arm64.
Compiler: `rustc 1.97.1 (8bab26f4f 2026-07-14)`. Cargo release profile is opt3,
thin LTO, one codegen unit, debug=false; both arms use `neighbor-diagnostics`.
The dedicated target directory preserved prior binaries. Build receipts record
the actual B04 revision/source/compiler variables and selected environment;
RUSTFLAGS, encoded flags, compiler wrappers, target and release-profile overrides
were absent. No applicable Cargo configuration file was found during the
provenance inspection. The controller reserved team CPU/write access; no user
process was stopped and the lease marker is not a machine-level lock.

### Frozen identities and local archive

| Artifact | SHA256 |
| --- | --- |
| Source identity (tracked Rust/Cargo and harness bytes) | `59f1544ad6351c8330bbbf66faaeefd3cb9228cf63eb18ea1121c246667f1e21` |
| Named release executable | `24ab6722a6fe4578907bedcdb5cffdf1405784bcda47673f09be0332791e4f91` |
| Campaign manifest (15 snapshot hashes /116 cases) | `8908fc21bc294de14a97e5a470ab3bcc6a341d25e6d4c5b1b3d56e9a202dbc9c` |
| Integrated analysis | `8a6117f056806269a782926fc9ceaac8fa2af3d585a75adc4b3975af3954fc34` |
| Controller freeze manifest | `d73e1e4d68f994d2a162789ee57e72aa0864a9b99f251a670d90a3da7da7b3bc` |
| Additive raw artifact manifest | `b7080c8fdf973dfd0891c694af005c5f137b4930401f1211714f4249382db0bd` |
| Source tar at clean f56f288 | `b7a2120d2554b623cd40372cedd04e35c6a9be825adaaa2f9f2ea89d91ab8d34` |
| Raw/source/executable archive | `e04e877c522bbac5fd4fb1796e7a2a47ec3723ed6af0ba628198175413825159` |

Stable local archive:
`/Users/ricardoleal20/Documents/Code/Misc/modena/target/b04-evidence/task8-g4-f56f288-20261106/raw-campaign-f56f288.tar.gz`.
Raw files remain in `target/b04-triclinic-campaign/`; receipts, independent
controller freeze, source tar, frozen executable and additive manifests remain
under `target/b04-evidence/task8-g4-f56f288-20261106/`. These ignored local
artifacts are unavailable in a fresh clone. Older prepared/smoke/failed epochs
were preserved, not overwritten or relabeled.

The integrated analysis was independently repeated using `--analyze-only` on
a byte-for-byte mirror. All numerical/status fields match; only35 missing-file
error directory prefixes differ. The archive helper's initial full-JSON equality
assertion failed on those expected path differences; its failed helper and both
analyses are retained. A path-only normalization then verified exact equality.
This correction did not change the harness, raw data, cases or timings.

Actual command vectors, raw stdout/stderr, UTC starts, durations and exits are
in numbered JSON receipts. The retained private wrappers are also archived:

```sh
rtk proxy python3 /private/tmp/b04-task8.py build
rtk proxy python3 /private/tmp/b04-task8.py run
rtk proxy python3 /private/tmp/b04-task8-freeze.py
rtk proxy python3 /private/tmp/b04-task8-checks.py
```

Build invokes `cargo build --locked --release -p modena-bench --features
neighbor-diagnostics --bin triclinic_neighbors` with actual B04 compiler
variables and a new CARGO_TARGET_DIR. Generation/run use the dedicated Python
driver with the frozen binary and fresh campaign directory. `--run` executes
once; standalone analysis uses `offline-analysis/` because writes are exclusive.
No timing restart occurred after the execution session became unavailable:
recovery read the terminal receipt and existing analysis. The controller froze
those files and released the lease before documentation and final checks.

Task7 preparation and I1 correction are independently approved through `f56f288`;
that approval adds no G4 acceptance. G1/G2/G3 retain their individual approved
software scopes and historical execution epochs. G2 remains initialized MACE
software parity only, with its disclosed dependency warnings; no FAIR-Chem
task/checkpoint, trained chemical accuracy or B10 acceptance is added.
The documentation commit changes revision while preserving measured source
bytes; the retained binary's runtime revision guard remains bound to `f56f288`.
Final whole-branch review and PR delivery remain controller-owned.

## Controller decisions and costs preserved from execution

These eleven decisions were made at the recorded execution stages, before
the affected evaluation. They preserve the controller's literal ruling scope,
rationale and costs from the private execution ledger for later public review.
They add no performance acceptance; G4 remains not passed.

1. **Independent B04 start and later B03 synchronization.** B04 numerical work
   began in a separate checkout while B03 closed derived reporting/review,
   because the immutable reviewed core was identical and B03's measurement
   lease had been explicitly released. The cost was mandatory synchronization
   of the unpublished B04 branch with final reviewed B03 before publication,
   with source changes verified, no force-push and no unauthorized merge.
   Only one B04 implementer could work at a time; independent B03 metadata
   closure could not edit B04 files. The later source/patch identity and commit
   mapping receipts above record the completed synchronization.

2. **Bounded f64 geometry with explicit rejection.** The graph certifies its
   executable f64 predicate and all-image coverage for accepted requests,
   rejecting uncertainty instead of silently omitting images. The cost is
   that some huge, skew or ill-conditioned finite inputs reject, and interval
   arithmetic adds implementation cost. Search bounds never widen the final
   physical cutoff predicate or return a truncated graph.

3. **Initialized official MACE software scope.** An actual initialized official
   MACE instance tests graph/API/math software parity because B04 is not a
   trained-chemistry validation. The cost is that this establishes no chemical
   accuracy or FAIR-Chem checkpoint parity. B10 still requires its own actual
   domain, model and reference validation.

4. **Unconditional tuple guards under Python optimization.** The plan's Python
   assert snippet was replaced with an unconditional mismatch exception and
   an optimization-enabled negative control: exact integer tuple identity
   remains mandatory under `PYTHONOPTIMIZE`/`-O`. The cost is a small explicit
   verifier/test divergence from the snippet and maintenance of the error path.
   Historical ASE evidence was retained, and the changed helper received a new
   source/evaluation epoch; no performance retiming was authorized by this fix.

5. **G2 NumPy compatibility revision before evaluation.** Only NumPy 2.2.6 was
   revised to 1.26.4 because actual matscipy 1.1.1 requires NumPy<2 and official
   NumPy 1.26.4 supports Python 3.12. The 44-package hashed candidate and direct
   binary dry-sync checks preceded evaluation. The cost is a different NumPy
   ABI/arithmetic baseline requiring actual import, E/F and finite-difference
   revalidation; no numerical comparison with the superseded, unexecuted NumPy
   2 candidate is claimed. Failed original pins/logs remain preserved. All other
   direct pins, model, seed, fixtures and gates stayed fixed. The published,
   hashed python-hostlist 2.3.0 sdist build was authorized with build provenance;
   no implicit framework/Python revision was allowed. Spec/plan/brief were
   updated before execution; compatibility resolution itself was not G2.

6. **Scoped `builtins.slice` allowlist for pinned e3nn constants.** Only
   `builtins.slice` is temporarily allowed through
   `torch.serialization.safe_globals` around e3nn 0.4.4 imports of its
   hash-recorded packaged `constants.pt`; Torch 2.6's restricted weights-only
   loader otherwise rejects that library constant. The cost is trust in the
   pinned packaged constant/source provenance and future import compatibility
   maintenance. Restoration of the prior allowlist was verified. Unrestricted
   pickle loading, global `TORCH_FORCE_NO_WEIGHTS_ONLY_LOAD` and additional
   allowed types without review were excluded. Pins/model/gates stayed fixed;
   imports were not model acceptance. Package, constant-file SHA, source and
   runtime receipts are retained. The audited primary API source is PyTorch
   v2.6.0 `torch/serialization.py` (`safe_globals`); no retrieved model/checkpoint,
   authorization or license gate was involved in this library-constant import.

7. **Deterministic change-shift mutation repair.** The first nonzero shift
   component now increases in absolute value. The original first edge
   `[-1,0,0]` plus x=1 produced a forbidden zero self edge and nonfinite MACE
   forces, making that control invalid rather than demonstrating model-parity
   failure. The cost is coverage of this fixed away-from-zero mutation, not
   arbitrary shift mutations. The first failed gate and preflight model digest
   were retained; the replacement rule was fixed before rerun. No seed/model/
   fixture/tolerance tuning occurred, and all three mutation classes plus
   finite differences, repeat stability and official/native parity remained
   mandatory.

8. **Private bounded builder and orthorhombic snapshot facade.** The bounded
   builder stays crate-private and the public snapshot facade stays
   orthorhombic, preserving existing visibility/domain while proving actual
   triclinic invariants. Real triclinic bounded, zero-skin, cap, failure and
   recovery checks live in private neighbor units; integration controls use
   the actual public orthorhombic snapshot path. The cost is coverage split
   across unit/integration locations and no new public bounded-triclinic
   snapshot capability.

9. **Explicit elongated matrix before generation.** The elongated family
   retains moderate off-diagonals `(.2L,.1L,.15L)` with diagonal factors
   `(2,.5,1)`, resolving the draft's unspecified shear so all four families
   have explicit frozen triclinic matrices. The cost is that anisotropic
   evidence covers this one declared shape/shear, not arbitrary elongated
   cells. This was fixed before snapshot generation/performance; determinant
   L^3, density, N, velocities, dt, blocks, caps and criteria stayed unchanged.

10. **Benchmark-local counting allocator.** Existing private capacities were
    not replaced by fabricated estimates or exposed through a new production
    API. Allocation count is actual timed process-wide Rust allocation and
    reallocation attempts; retained bytes are live Rust-requested heap after
    the block, including both replicas, pool and benchmark data; requested
    bytes are diagnostic temporary neighbor workspace only. The cost is
    instrumentation overhead and no per-model/per-strategy memory ownership
    or OS-RSS attribution. These scopes were frozen/disclosed before timing,
    with matched immutable timing/control mode/source hashes. No uninstrumented
    production claim or zero-total-memory inference from zero temporary
    workspace is allowed. Failed allocation attempts are included in the
    final instrumented implementation.

11. **Separate boundary recipe before accepted timing.** Only the separate
    boundary control was corrected: N=8, moderate L=32 cell; isolated origin
    atom 0 with zero velocity; seven mobile Cartesian cube corners, indices
    1..7, spacing 1.4 and anchor H*(.5,.5,.5); existing velocity formula applied
    only to that cluster, then its component means removed. The density-.8
    ladder stayed unchanged. The reason was to preserve FractionalBoundary
    under continuous dynamics without physical resets and avoid the earlier
    AllBins preflight. The cost is a separate low-density control, not dense
    production boundary-performance evidence. An unsupported recipe had to
    stop with diagnostics; all older prepared/smoke bytes were retained.
    Forced-fractional grid-cap/boundary controls require the intended actual
    initial reason and every measured new build's reason, distinct from Auto
    dispatch. A new immutable epoch was generated after the fix; no old timing
    was relabeled. Task7 correctness checks covered this rule, while Task8's
    later timeout left the campaign controls unrun and their G4 coverage absent.


### Final repository checks — actual Task8 execution

After the controller froze the raw campaign and released its lease, the full
default-enabled core suite plus `neighbor-diagnostics` passed **364 tests, zero
failures**, across 34 reported suites. It included the actual graph/triclinic,
bounded/skin/pair-policy and compile-fail borrow regressions. Four tests were
ignored in that main command: the three explicit official gates and one
pre-existing configuration-builder doc example. The three official tests were
then explicitly executed: **3 passed, zero ignored**, including 27 ASE fixtures,
the optimized wrong-tuple control and actual initialized MACE graph/autograd/
finite-difference/mutation guards. There were no optional `SKIP` markers in the
core output. GPU/MPI features were not enabled or validated by these checks.

Python harness: **23 passed**. Rust harness: **5 passed**. Core/harness Clippy
with `-D warnings` and full workspace formatting passed. The API guide's actual
Rust example passed **1/1** after binding to Cargo's exact stable library
artifact. Its first check selected a cached Rust 1.90 rlib by modification time
and failed E0514 against Rust 1.97.1; that failed receipt is retained. A narrow
`cargo build --message-format=json` identified the correct artifact, followed
only by the documentation recheck. No cargo clean, broad test repetition,
production source change or timing rerun was performed.

Native compiler and Clippy checks emitted no warning diagnostics. The retained
full-core proptest receipt contains five `FileFailurePersistence` diagnostics
(`proptest: FileFailurePersistence::SourceParallel set, but failed to find
lib.rs or main.rs.`); these do not change the check results. The historical upstream
MACE/TorchScript annotation warnings remain disclosed and were not resolved by
changing pins or suppressing warnings. The official driver captures helper
stderr internally on success; absence from its success log does not establish
that the dependency warning was removed. Existing pinned environments and
dependency/source receipts were reused; no weights/packages were downloaded.
Stable final checks are fresh; earlier Rust 1.90 evidence retains its own
historical execution epoch and is not relabeled as a Task8 rerun.

This documentation-only revision clarifies the warning receipt and diagnostic
API rustdoc. It does not change the measured frozen `f56f288` source/binary
epoch or add numerical validation; G1–G3 remain historical/fresh-at-Task8
results, and G4 remains not passed with Auto Scalar.

Unchanged thresholds: graph vectors atol=rtol=2e-12; MACE E/F atol=rtol=1e-10 and finite differences 1e-6; classical component normalization threshold 1e-10; initial cutoff margin >1e-8. No tolerance was widened.

Fresh initialized-MACE state SHA256 remains `1d1a81cf8cfaa7885469757c26b82d49944cd188034ba7f5c481ad020974066e`
(1968 parameters). Maximum native/official E error is
5.5511151231257827e-17; maximum force
component error is 1.717376241217039e-16;
the nine-component central-difference error is 1.6931400725894719e-11.
Repeat differences are zero. All three fixed mutation classes and optimized
tuple/fractional/schema negative controls passed. The helper/fixture/lock and
dependency provenance hashes remain explicit in receipt 11 and retained preflight
manifests. This is scoped software parity, with no trained-chemical/B10 claim.

Actual child command vectors below are rendered from the JSON receipts; every
agent shell entrypoint used RTK, while child/runtime argv have no RTK dependency.
Receipt directories are relative to
`target/b04-evidence/task8-g4-f56f288-20261106/`. Each JSON retains raw output,
exact cwd, UTC start, wall duration and exit. Receipt 17 is the preserved
documentation artifact-selection failure; receipt 19 is its narrow passing fix.

| Receipt | Actual inner argv | Exit |
| --- | --- | --- |
| `01-build.json` | `cargo build --locked --release -p modena-bench --features neighbor-diagnostics --bin triclinic_neighbors` | 0 |
| `02-generate.json` | `python3 bench/triclinic_neighbors.py --binary /Users/ricardoleal20/.codex/worktrees/modena-b04/modena/target/b04-evidence/task8-g4-f56f288-20261106/frozen/triclinic_neighbors --output-dir /Users/ricardoleal20/.codex/worktrees/modena-b04/modena/target/b04-triclinic-campaign --generate-only` | 0 |
| `03-run.json` | `python3 bench/triclinic_neighbors.py --binary /Users/ricardoleal20/.codex/worktrees/modena-b04/modena/target/b04-evidence/task8-g4-f56f288-20261106/frozen/triclinic_neighbors --output-dir /Users/ricardoleal20/.codex/worktrees/modena-b04/modena/target/b04-triclinic-campaign --run` | 1 |
| `04-analyze-only.json` | `python3 bench/triclinic_neighbors.py --output-dir /Users/ricardoleal20/.codex/worktrees/modena-b04/modena/target/b04-evidence/task8-g4-f56f288-20261106/offline-analysis --analyze-only` | 1 |
| `05-source-archive.json` | `git archive --format=tar --output=/Users/ricardoleal20/.codex/worktrees/modena-b04/modena/target/b04-evidence/task8-g4-f56f288-20261106/frozen/source-f56f288.tar f56f28854865e8068402eebf7a0548f60da8ecfd` | 0 |
| `06-host.json` | `sysctl hw.model machdep.cpu.brand_string hw.ncpu hw.physicalcpu hw.logicalcpu hw.memsize` | 0 |
| `07-binary-format.json` | `file /Users/ricardoleal20/.codex/worktrees/modena-b04/modena/target/b04-evidence/task8-g4-f56f288-20261106/frozen/triclinic_neighbors` | 0 |
| `10-full-core.json` | `env MODENA_ASE_PYTHON=/Users/ricardoleal20/Documents/Code/Misc/modena/target/b01-ase-venv/bin/python cargo test --locked -p modena-core --features neighbor-diagnostics -- --nocapture` | 0 |
| `11-official-gates.json` | `env MODENA_ASE_PYTHON=/Users/ricardoleal20/Documents/Code/Misc/modena/target/b01-ase-venv/bin/python MODENA_MACE_PYTHON=/Users/ricardoleal20/.codex/worktrees/modena-b04/modena/target/b04-mace-venv/bin/python MODENA_MACE_EVIDENCE=/Users/ricardoleal20/.codex/worktrees/modena-b04/modena/target/b04-evidence/task8-g4-f56f288-20261106/official-mace cargo test --locked -p modena-core --features neighbor-diagnostics --test periodic_graph_official -- --ignored --nocapture` | 0 |
| `12-python-harness.json` | `python3 -m unittest discover -s bench/tests -p test_triclinic_neighbors.py` | 0 |
| `13-rust-harness.json` | `cargo test --locked -p modena-bench --features neighbor-diagnostics --bin triclinic_neighbors` | 0 |
| `14-core-clippy.json` | `cargo clippy --locked -p modena-core --features neighbor-diagnostics --lib --tests -- -D warnings` | 0 |
| `15-harness-clippy.json` | `cargo clippy --locked -p modena-bench --features neighbor-diagnostics --bin triclinic_neighbors -- -D warnings` | 0 |
| `16-fmt.json` | `cargo fmt --all -- --check` | 0 |
| `17-api-doctest.json` | `rustdoc --edition=2021 --test docs/B04_PERIODIC_GRAPHS.md --extern modena_core=/Users/ricardoleal20/.codex/worktrees/modena-b04/modena/target/debug/deps/libmodena_core-94d2d062254895a0.rlib -L dependency=/Users/ricardoleal20/.codex/worktrees/modena-b04/modena/target/debug/deps` | 101 |
| `18-doc-library.json` | `cargo build --locked -p modena-core --features neighbor-diagnostics --lib --message-format=json` | 0 |
| `19-api-doctest-final.json` | `rustdoc --edition=2021 --test docs/B04_PERIODIC_GRAPHS.md --extern modena_core=/Users/ricardoleal20/.codex/worktrees/modena-b04/modena/target/debug/libmodena_core.rlib -L dependency=/Users/ricardoleal20/.codex/worktrees/modena-b04/modena/target/debug/deps` | 0 |

The build additionally set B04_BUILD_REV to clean f56f288 (full revision above),
B04_BUILD_SOURCE_SHA256 to 59f1544a… and B04_RUSTC_VERSION to the recorded full
rustc string; CARGO_TARGET_DIR was the dedicated evidence `build/` directory.
`build-environment.json` retains exact values and each source file's SHA256.
The run wrapper set MODENA_B04_CPU_LEASE=exclusive. The archive's first helper
assertion and its corrected path-only comparison are retained separately from
the expected analyzer exit 1 for an incomplete campaign.

Final independent rehashing verifies every controller-freeze entry, every
additive raw-manifest entry and every corresponding file inside the archive.
No original campaign or historical evidence file changed. Source bytes remain
identical to clean f56f288; only the three owned Markdown documents changed.

The final documentation diff check passed (receipt `20-diff-check.json`, actual inner argv `git diff --check`, exit0, empty stdout/stderr).

Final additive check/receipt manifest: `target/b04-evidence/task8-g4-f56f288-20261106/final-check-manifest.json`, SHA256 `b49ca8263604b738f879842da559240cb6977f7f2e4d01b9bc4703622e3bacd2`. It binds 233 evidence files, including all check receipts, official preflight/source receipts and preserved failed attempts. The prior raw/controller manifests remain unchanged.

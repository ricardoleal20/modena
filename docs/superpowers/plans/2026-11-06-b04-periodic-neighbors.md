# B04 Periodic Graphs and Triclinic Neighbors Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use `subagent-driven-development` to implement this plan task-by-task. The controller assigns one implementer per task and obtains independent review. B04 tasks remain sequential; B04 implementers must not run concurrently. Independent B03 reporting/review occurs in its separate checkout after its measurement lease release. The controller performs Task 0 as preflight, without a separate scaffold implementer/review. Steps use checkbox syntax for tracking.

**Revision 2:** mandatory B-independent global shift certificate closes the initial draft's circular range assumption. Read the revised Spec global proof before Task 1. No final-cutoff change or truncation is permitted.

**Goal:** Implement complete bounded directed periodic graphs, prove their software parity using independent geometry and official MACE, and enable a triclinic classical stencil only after correctness and complete-step cost gates.

**Architecture:** A stateless graph borrows an immutable System snapshot and publishes all qualifying directed integer-image tuples atomically. Private interval/certified-cell helpers bound enumeration and fractional binning. The existing classical half-list remains a separate contract with scalar fallback, exact pair ordering and a measured selector.

**Tech Stack:** Rust 2021/MSRV 1.90, existing thiserror/serde_json/Rayon, official ASE reference, isolated Python 3.12.8 + pinned official MACE 0.3.13/PyTorch CPU float64.

**Spec:** `docs/superpowers/specs/2026-11-06-b04-periodic-neighbors-design.md`, to be copied from `/private/tmp/modena-b04-spec-ready.md` only after final B03; temporary draft companion is `/private/tmp/modena-b04-plan-ready.md`.

## Global Constraints

- All new branches must use the user-required `ricardo/` prefix; this overrides the desktop default `ricardoleal20/`.
- Execution uses subagent-driven-development with one implementer per task and independent review; tasks are sequential and implementers must not run concurrently. Controller preflight owns reconciliation/worktree/document setup without a separate scaffold agent/review.
- Rust edition 2021, workspace rust-version 1.90; no fast-math or reassociation of specified geometry expressions.
- Physical positions/cutoff use the same units; reference fixtures are Angstrom and model energies/forces use their declared numerical units.
- Graph membership uses strict squared Cartesian cutoff with the exact operation order below; no physical epsilon shell, top-K, closest-image-only rule, classical exclusions or skin.
- Every emitted image shift is an exact `[i64; 3]`; output is lexicographically ordered and duplicate-free.
- Explicit finite resource budgets, checked integer arithmetic, fallible allocation and atomic graph publication are mandatory.
- Keep the existing classical half-list, scalar triclinic oracle, zero-skin rebuilding, B09/snapshot bounded pair builder and orthorhombic behavior valid.
- Actual official-model graph consumption and E/F parity are mandatory; EMT transport, LJ and invented toy networks cannot replace G2.
- Pin model settings, snapshots, tolerances and dependency manifest before numerical evaluation; record actual dependency hashes and platform.
- Triclinic stencil default selection requires G3 and G4; retain scalar selection for unsupported, inconclusive or unprofitable cases.
- Numerical tests may run after the documented B03 CPU lease release; acquire an exclusive measurement reservation before any performance campaign. Do not overlap team builds/tests with timed runs.

All paths below are relative to the isolated B04 worktree, never the B03 closure checkout. All commands below are future execution instructions, not checks performed during document preparation. Generated venvs, binaries, snapshots and timings live under ignored target/b04-* directories. Commit messages follow gitmoji + imperative Action; append a real Linear trailer only when an actual ticket ID is supplied, never invent one.

## File map

| File | Responsibility |
| --- | --- |
| `crates/modena-core/src/periodic_geometry.rs` | Private interval arithmetic, inverse residual certificate, image bounds, scalar-predicate envelopes. |
| `crates/modena-core/src/periodic_graph.rs` | Public graph/options/budgets/errors/stats and atomic two-pass construction. |
| `crates/modena-core/src/lib.rs` | Public graph exports and private geometry module declaration. |
| `crates/modena-core/tests/periodic_graph.rs` | Graph contract, resource and mutation integration tests. |
| `crates/modena-core/tests/support/periodic_oracle.rs` | Fixed-cube direct scalar oracle; no inverse/min_image reuse. |
| `crates/modena-core/tests/fixtures/periodic_graph/cases.json` | Frozen small geometry/model snapshots and explicit expected tuples. |
| `crates/modena-core/tests/fixtures/periodic_graph/README.md` | Independent numeric cube certificates and precision contract. |
| `crates/modena-core/tests/periodic_graph_official.rs` | Actual Rust edges to pinned official ASE/MACE subprocesses. |
| `crates/modena-core/tests/support/ase_graph_reference.py` | ASE tuple/vector parity and package/source identity. |
| `crates/modena-core/tests/support/mace_parity.py` | Official initialized MACE, actual graph consumption, E/F/autograd/mutation parity. |
| `bench/requirements-b04-mace.in`, `bench/requirements-b04-mace-py312.lock` | Exact initial direct pins and hashed resolved dependency environment. |
| `crates/modena-core/src/neighbor/triclinic.rs` | Scalar/stencil strategy, certified fractional bins/stencil/caps/fallback. |
| `crates/modena-core/src/neighbor.rs` | Narrow triclinic dispatch, existing skin/bounded/exclusion integration. |
| `crates/modena-core/tests/triclinic_neighbors.rs` | Scalar/stencil pair, E/F/W, zero-skin and snapshot regressions. |
| `crates/modena-core/Cargo.toml`, `crates/modena-bench/Cargo.toml` | Opt-in neighbor-diagnostics feature and dedicated benchmark bin feature requirement. |
| `crates/modena-bench/src/bin/triclinic_neighbors.rs` | One-build scalar/stencil rebuild and complete-step benchmark, exact input/output manifest. |
| `bench/triclinic_neighbors.py`, `bench/tests/test_triclinic_neighbors.py` | Campaign control, paired bootstrap analysis and result-schema validation. |
| `docs/superpowers/specs/2026-11-06-b04-periodic-neighbors-design.md`, `docs/superpowers/plans/2026-11-06-b04-periodic-neighbors.md` | Reviewed contract and tracked execution state. |
| `docs/B04_PERIODIC_GRAPHS.md`, `docs/B04_VERIFICATION.md` | User/API limits and separately evidenced G1-G4 outcomes. |

## Task 0: Controller preflight — bind reviewed core, isolate B04 and freeze the protocol

**Owner:** Controller preflight. This setup is not dispatched to an implementer and has no separate scaffold-agent review; independent implementation reviews begin with Task 1.

**Files:** Create the spec/plan documentation paths above and `docs/B04_VERIFICATION.md`. No production source change in this preflight.

**Consumes:** reviewed B03 core at 0c54db3, its identical core-tree certificate, and explicit release of the CPU lease. **Produces:** recorded BASE_REV, B04 branch/worktree, exact source reconciliation and frozen geometry/model/timing contract.

- [x] Confirm B03 measurement CPU lease is released and its reviewed core tree matches the B04 base before numerical implementation. Independent numerical work may proceed while B03 reports/reviews close; synchronize the unpublished branch with final reviewed B03 before publication. Read the applicable using-git-worktrees skill at that time. Create a fresh B04 worktree from reviewed B03 0c54db3, then synchronize with final reviewed B03 before publication on a branch prefixed `ricardo/` (the user-required override); never reset or reuse the live B03 worktree.
- [x] Record source reality and compare only the relevant symbols against this draft:

```sh
rtk git log -1 --format='%H%n%s'
rtk git status --short
rtk rg -n 'fn rebuild|rebuild_bounded|rebuild_impl|fn wrap_half|fn min_image|inv_basis_rows|pub struct System|positions_raw' crates/modena-core/src/neighbor.rs crates/modena-core/src/system.rs
rtk read crates/modena-core/Cargo.toml
rtk read crates/modena-bench/Cargo.toml
```

The source review must establish actual min_image/Vec3 operation order, scalar pair order, bounded failure behavior, zero-skin path and B03 neighbor changes. If any differs materially, update the spec/proof first. Do not assume the audited aec5bfb revision is final B03.

- [x] Copy the two draft files into the committed spec/plan locations, update their status to reviewed implementation contract after parent review, set Spec to the actual repository spec path, and record BASE_REV/source facts in verification. Freeze every fixture, tolerance, model constructor/seed and campaign criterion before collecting numerical evidence.
- [x] Record four independent statuses, initially pending:

```markdown
| Gate | Status | Required evidence |
| --- | --- | --- |
| G1 | pending | Fixed cube, resource failures, official ASE exact tuples |
| G2 | pending | Actual official initialized MACE graph and E/F parity |
| G3 | pending | Scalar/stencil pairs, E/F/W, skin/bounded regressions |
| G4 | pending | Paired complete-step benefit and fallback overhead |
```

- [x] Commit only reviewed documents:

```sh
rtk git add docs/superpowers/specs/2026-11-06-b04-periodic-neighbors-design.md docs/superpowers/plans/2026-11-06-b04-periodic-neighbors.md docs/B04_VERIFICATION.md
rtk git commit -m ':memo: Add: B04 periodic neighbor design contract'
```

**Gate:** no implementation, installation or timing starts before source/protocol reconciliation and CPU release. Spec changes resulting from a real incompatibility must be explicit, with no post-result tolerance broadening.

## Task 1: Implement certified inverse and explicit floating geometry

**Files:** Create `crates/modena-core/src/periodic_geometry.rs`; declare private module in `lib.rs`. Unit tests live in this file so no public geometry API is created.

**Consumes:** SimBox basis/inverse getters and exact spec formulas. **Produces:** `Interval`, `CertifiedCell`, certified_cell, global_shift_bound, image_bounds, edge_displacement, strict_radius and the shared scalar_min_image_search_radius used in Task 6. Public graph option/error types arrive in Task 2; define the independent private numeric core first with `GeometryError` and `GeometryOptions { cutoff, pbc, max_abs_shift }`, then map it into the public GraphBuildError in Task 2. Exact private types:

```rust
#[derive(Clone, Copy, Debug)]
pub(crate) struct Interval { pub(crate) lo: f64, pub(crate) hi: f64 }
#[derive(Clone, Copy, Debug)]
pub(crate) enum GeometryError { Invalid, Uncertified, ShiftRange }
#[derive(Clone, Copy)]
pub(crate) struct GeometryOptions {
    pub(crate) cutoff: f64,
    pub(crate) pbc: [bool; 3],
    pub(crate) max_abs_shift: i64,
}
```

Geometry functions return standard `Result<T, GeometryError>` internally, as defined in Spec. Task 2 maps the private errors at the public graph boundary. `CertifiedCell` fields and functions are crate-private; no source dependency on `periodic_graph` is introduced.

- [ ] Write unit cases `interval_contains_exact_dyadics`, `identity_inverse_certificate_encloses_zero`, `skew_inverse_certificate_encloses_triangular_inverse`, `uncertifiable_inverse_fails`, `normal_cutoff_required`, `integer_endpoints_checked_before_cast`, `global_bound_is_independent_of_configured_B`, `global_contraction_rejects_uncertifiable_cancellation`, `huge_coordinates_reject_before_enumeration`, `all_nonperiodic_global_shift_bound_is_zero`, `i64_conversion_error_is_in_global_proof`, and `cartesian_operation_order_boundary`. Hand-computable checks include identity, a=(2,0,0),b=(1,2,0),c=(0,0,4), whose U rows are (.5,-.25,0),(0,.5,0),(0,0,.25).

```rust
#[test]
fn interval_contains_exact_dyadics() {
    let a = Interval::point(0.5).unwrap();
    let b = Interval::point(0.25).unwrap();
    let sum = a.add(b).unwrap();
    assert!(sum.lo <= 0.75 && sum.hi >= 0.75);
    let product = a.mul(b).unwrap();
    assert!(product.lo <= 0.125 && product.hi >= 0.125);
    let subnormal = Interval::point(f64::from_bits(1)).unwrap();
    let half = subnormal.mul(a).unwrap();
    assert!(half.lo <= 0.0 && half.hi >= f64::from_bits(1));
}
```

- [ ] Run red:

```sh
rtk cargo test --locked -p modena-core --lib periodic_geometry
```

Expected failure: missing private module/types/functions. Do not treat a dependency-resolution failure as the intended red test.

- [ ] Implement point/add/sub/mul/div/sqrt/hull/abs_upper with outward operations and finite/zero-domain checks; square handles intervals spanning zero. Implement residual rho<.5/eta/inverse-row norm certificate and the Spec B-independent global range proof first: U_P/H_P, gamma_10 conversion-aware error, upper alpha<.5, K_global, and G<=B<=2^52 before casting. Evaluate all rational expressions with lower denominator enclosures, including 1-alpha and 1-gamma_m; test that an upper-rounded denominator cannot be substituted. Then implement the exact scalar expression and local Q_enum. image_bounds must invoke the global certificate internally, so no caller can apply a B-dependent envelope without proving global exclusion. All comparisons deciding to omit an image must use conservative interval endpoints. Implement checked integer conversions before any `as i64`. Include proof comments with the Neumann-series inequality, global ten-rounding cast-aware inequality, contraction rearrangement and local seven/five operation counts. Explain why hypothetical accepted i64 shifts outside the exactly representable domain are covered by the global bound, and why physically required integer images outside i64 are excluded for accepted requests. Use no vectors in the numeric core.

```rust
pub(crate) fn checked_shift_endpoint(x: f64, max_abs: i64)
    -> Result<i64, GeometryError>
{
    if !x.is_finite() || max_abs < 0 || max_abs > (1_i64 << 52)
        || x < -(max_abs as f64) || x > max_abs as f64
        || x.fract() != 0.0
    {
        return Err(GeometryError::ShiftRange);
    }
    Ok(x as i64)
}
```

Keep this helper separate from ceil/floor so conversion cannot silently saturate. For nonperiodic axes return `(0,0)` after input/cell certification regardless of the periodic bound.

- [ ] Run the targeted command green. Manually reconcile operation counts against actual compiler-safe Rust expression order; test perturbations of the dyadic cell and exact next-up/down cutoff. R2=R*R must be finite positive normal; tests explicitly reject overflow to infinity, zero underflow and nonzero subnormal squares on empty and nonempty inputs. For global-bound regression use identical geometry with two large B values and assert the returned G is identical, then B=G-1 fails ShiftRange. Prove by an inequality assertion that every fixed-cube passing tuple has shift norm<=G. Exercise i64 cast endpoints analytically with values i64::MIN, i64::MAX and neighbors of 2^52/2^53; no exhaustive huge-shift loop or large allocation is needed. Test huge common translation that loses certification returns a numerical/range error rather than an empty/partial graph. Stop if residual/error enclosures cannot be established.
- [ ] Commit numeric module and source declaration:

```sh
rtk git add crates/modena-core/src/periodic_geometry.rs crates/modena-core/src/lib.rs docs/superpowers/specs/2026-11-06-b04-periodic-neighbors-design.md
rtk git commit -m ':sparkles: Add: certified periodic geometry bounds'
```

**Gate:** neither arbitrary EPS margins nor a B-only local error bound may substitute for the global exclusion proof. Accepted requests must prove all real physical images and every i64 tuple accepted by the f64 predicate fit G<=B<=2^52. Numerical/range rejection must be explicit and deterministic.

## Task 2: Implement atomic directed all-image graph and resource policy

**Files:** Create `crates/modena-core/src/periodic_graph.rs` and `tests/periodic_graph.rs`; modify `lib.rs`. Public types/signatures are exactly the spec API. Map GeometryError::Invalid to InvalidInput, Uncertified to NumericalGeometry and ShiftRange unchanged.

**Consumes:** certified numeric core. **Produces:** public borrowed-snapshot graph, stats and explicit failure classes; no classical neighbor changes.

- [ ] Add test helpers and red contract tests. These helpers are complete small Rust scaffolding:

```rust
use modena_core::{GraphBudget, PeriodicGraph, PeriodicGraphOptions, SimBox, System, Vec3};
fn options(cutoff: f64, include_self: bool) -> PeriodicGraphOptions {
    PeriodicGraphOptions {
        cutoff, pbc: [true; 3], include_periodic_self_images: include_self,
        budget: GraphBudget {
            max_atom_pairs: 1_000_000, max_image_candidates: 10_000_000,
            max_edges: 1_000_000, max_bytes: 134_217_728,
            max_abs_shift: 1_000_000,
        },
    }
}
fn system(cell: SimBox, positions: &[Vec3]) -> System {
    let mut s = System::new(cell);
    s.resize(positions.len());
    for (i, &r) in positions.iter().enumerate() { s.set_pos(i, r); }
    s
}
#[test]
fn multi_image_graph_has_exactly_four_directed_tuples() {
    let s = system(SimBox::cubic(2.0), &[Vec3::zero(), Vec3::new(0.6,0.0,0.0)]);
    let graph = PeriodicGraph::build(&s, options(1.5, true)).unwrap();
    let actual: Vec<_> = graph.edges().iter()
        .map(|e| (e.source,e.destination,e.edge_cell_shift)).collect();
    assert_eq!(actual, vec![(0,1,[-1,0,0]),(0,1,[0,0,0]),
        (1,0,[0,0,0]),(1,0,[1,0,0])]);
    assert_eq!(graph.stats().atom_pair_visits, 8);
}
```

Also name tests `periodic_self_policy_and_zero_self`, `all_pbc_masks`, `strict_cutoff_next_float`, `malformed_geometry_fails_empty_and_nonempty`, `budget_failure_preserves_source_and_prior_graph`, `shift_budget_is_not_top_k`, `global_range_failure_is_atomic`, `configured_B_never_defines_completeness`, `graph_borrow_prevents_mutation` (compile-fail doc test), and `fresh_borrow_after_identity_and_box_mutations`.

- [ ] Run red:

```sh
rtk cargo test --locked -p modena-core --test periodic_graph
```

- [ ] Implement validation and two passes with `2*N*N` pair charging. For each pair first compute the B-independent global shift certificate and reject any unsupported whole request, then calculate certified local integer spans, check `2*running_candidates` before enumeration, then membership/self policy. Count accepted edges before allocation. Enforce max_edges, requested storage and `try_reserve_exact`. Repeat against borrowed immutable input, emit, sort_unstable, verify stats and publish. No fallible mutation to caller state is performed. No pair-scale lookup or minimum-image half list is consulted.

```rust
let pairs_one_pass = (n as u128).checked_mul(n as u128)
    .ok_or(GraphBuildError::BudgetExceeded { resource: "atom_pairs", limit: budget.max_atom_pairs })?;
let pair_visits = pairs_one_pass.checked_mul(2)
    .ok_or(GraphBuildError::BudgetExceeded { resource: "atom_pairs", limit: budget.max_atom_pairs })?;
if pair_visits > u128::from(budget.max_atom_pairs) {
    return Err(GraphBuildError::BudgetExceeded { resource: "atom_pairs", limit: budget.max_atom_pairs });
}
```

Use n conversion checked on targets whose usize could exceed u64. Per-axis span arithmetic is i128/u128; an empty span never wraps to huge unsigned length. The byte cap is logical requested storage as defined in spec, not a guarantee about allocator RSS.

- [ ] Add cfg(test)-only `allocate_edges_with` seam accepting a closure returning `Result<Vec<PeriodicEdge>,()>`; unit test injects failure without reserving huge real memory. Test each error resource exactly, input positions/forces/atomic_numbers unchanged, and a separate prior graph's edge bytes unchanged. Valid empty graph has zero pair/candidate/edge counters and fixed requested header/workspace bytes.
- [ ] Run targeted green and compile-fail doc tests; commit:

```sh
rtk cargo test --locked -p modena-core --test periodic_graph
rtk cargo test --locked -p modena-core --doc
rtk git add crates/modena-core/src/periodic_graph.rs crates/modena-core/src/lib.rs crates/modena-core/tests/periodic_graph.rs
rtk git commit -m ':sparkles: Add: bounded directed periodic image graphs'
```

**Gate:** exact expected tuples, no duplicates, reversal and precise budget counts. Any partial publication, silent truncation or loss of self images rejects the task.

## Task 3: Add fixed-cube independent oracle and official ASE graph parity

**Files:** Create fixture JSON/README, `tests/support/periodic_oracle.rs`, `tests/periodic_graph_official.rs`, `tests/support/ase_graph_reference.py`; extend `tests/periodic_graph.rs`.

**Consumes:** actual graph API. **Produces:** exact G1 evidence tied to fixed fixture bytes and pinned ASE identity.

Fixture wire shape is exact:

```json
{"id":"multi4","cell":[[2,0,0],[0,2,0],[0,0,2]],
 "positions":[[0,0,0],[0.6,0,0]],"atomic_numbers":[6,8],
 "pbc":[true,true,true],"cutoff":1.5,"periodic_self":true,
 "cube":16,"official":true}
```

Other fixtures follow the spec table exactly. Boundary fixtures set official=false. Wire edges are objects `{source,destination,edge_cell_shift}`. A subprocess request is `{fixture:<fixture object>,edges:[<edge objects>]}`. For MACE finite differences the Rust driver adds `finite_difference`, a list of `{atom,axis,h:1e-5,plus:{fixture,edges},minus:{fixture,edges}}` for every atom/Cartesian axis of displaced; each plus/minus edge list is built by the actual Rust API before the helper runs. The helper must use those lists, not substitute official graph construction for the alternate route. Response is `{fixture_id,edge_count,versions,source_sha256,max_vector_error}` for ASE and additionally `{model_digest,parameter_count,energy_error,max_force_error,mutation_deltas,finite_difference_error}` for MACE. Missing properties or nonzero exit fail the Rust test; capture bounded stderr in failures.

- [ ] Write oracle code that independently loops all integer x/y/z from -16 through 16, skips disabled axes and zero-self/self policy, computes direct Cartesian displacement in declared scalar order, and compares sorted tuples. It must not import production geometry or call min_image. Add README numeric certificates:

```text
unreduced: |Sz|<4.8/40<1 => Sz=0;
|Sy|<(4.8+.1)/10=.49 => Sy=0;
|Sx|<(4.8+19.7)/10=2.45 => |Sx|<=2.
oblique: Sz=0, |Sy|<(4.8+3)/10=.78 => Sy=0;
|Sx|<(4.8+6.37)/10=1.117 => |Sx|<=1.
mace3: |Sz|<(3.1+1.5)/2.5=1.84 => |Sz|<=1;
|Sy|<(3.1+1.3+.3)/2.3<2.044 => |Sy|<=2;
|Sx|<(3.1+.9+.7*2+.2)/2.4<2.334 => |Sx|<=2.
gauge adds <=3,2,0, so every required component remains <16.
```

Use a +1e-8 search-only cushion for these numeric certificates to dominate finite evaluation error without changing physical membership. Add direct certificates for cubic/self/boundary and deformed/displaced values. Do not use the dynamic reciprocal production bound to derive the cube.

- [ ] Run oracle tests red with at least one deliberately incomplete initial graph path or newly asserted certificate/boundary property:

```sh
rtk cargo test --locked -p modena-core --test periodic_graph
```

If production graph already passes every new test, document genuine test addition with no manufactured product defect; do not intentionally commit an incorrect graph just to obtain red.

- [ ] Implement ASE helper and ignored named integration test `official_ase_images`. Exact primary path:

```python
import hashlib, inspect, json, sys
import ase, numpy as np
from ase import Atoms
from ase import neighborlist
request = json.load(sys.stdin)
f = request['fixture']
a = Atoms(numbers=f['atomic_numbers'], positions=f['positions'],
          cell=f['cell'], pbc=f['pbc'])
i, j, shifts, vectors = neighborlist.neighbor_list('ijSD', a, f['cutoff'], self_interaction=False)
records = sorted(((int(x), int(y), tuple(map(int,shift))), vector)
    for x,y,shift,vector in zip(i,j,shifts,vectors)
    if f['periodic_self'] or x != y)
reference = [edge for edge, vector in records]
actual = sorted((e['source'], e['destination'], tuple(e['edge_cell_shift']))
                for e in request['edges'])
assert actual == reference
source = inspect.getsourcefile(neighborlist)
with open(source, 'rb') as stream:
    digest = hashlib.sha256(stream.read()).hexdigest()
assert ase.__version__ == '3.29.0'
assert np.__version__ == '2.5.4'
def vector(edge):
    x,y,shift = edge
    value = np.asarray(f['positions'][y])-np.asarray(f['positions'][x])
    for axis in range(3):
        value = value + np.asarray(f['cell'][axis])*shift[axis]
    return value
ref_vectors = np.asarray([v for edge,v in records]).reshape(-1,3)
act_vectors = np.asarray([vector(e) for e in actual]).reshape(-1,3)
np.testing.assert_allclose(act_vectors,ref_vectors,rtol=2e-12,atol=2e-12)
vector_error = float(np.max(np.abs(act_vectors-ref_vectors),initial=0.0))
json.dump({'fixture_id':f['id'], 'edge_count':len(reference),
           'versions':{'ase':ase.__version__,'numpy':np.__version__},
           'source_sha256':digest,'max_vector_error':vector_error}, sys.stdout)
```

The snippet compares ASE's independent vector path in matching tuple order against the declared reconstruction. Rust tests generate actual graph edges and send exact fixture bytes; helper exit errors on missing env/packages and tuple mismatches. First verify >=1e-8 cutoff margin using the fixed cube before invoking any official helper.

- [ ] Run explicit official gate after the CPU lease release:

```sh
rtk proxy env MODENA_ASE_PYTHON=/Users/ricardoleal20/Documents/Code/Misc/modena/target/b01-ase-venv/bin/python cargo test --locked -p modena-core --test periodic_graph_official official_ase_images -- --ignored --nocapture
```

Minimal Rust suite may leave this ignored, but G1 remains pending unless this command succeeds and its manifest is retained.

- [ ] Commit fixtures/helpers/tests and G1 outcome:

```sh
rtk git add crates/modena-core/tests/fixtures/periodic_graph crates/modena-core/tests/support/periodic_oracle.rs crates/modena-core/tests/support/ase_graph_reference.py crates/modena-core/tests/periodic_graph.rs crates/modena-core/tests/periodic_graph_official.rs docs/B04_VERIFICATION.md
rtk git commit -m ':white_check_mark: Add: independent periodic image graph oracles'
```

**Gate:** cube/ASE integer tuples must match exactly, not only distances/counts. No unresolved boundary margin or missing interpreter may be marked passed.

## Task 4: Resolve isolated Python 3.12 and run actual official MACE software parity

**Files:** Create `bench/requirements-b04-mace.in`, `bench/requirements-b04-mace-py312.lock`, `tests/support/mace_parity.py`; extend official Rust tests and verification docs.

**Consumes:** frozen fixtures/Rust graph/official subprocess schema. **Produces:** a pinned CPU environment and actual official-model graph/E/F evidence; no trained chemical claim.

- [ ] Write the exact direct requirements before resolution:

```text
torch==2.6.0
mace-torch==0.3.13
e3nn==0.4.4
numpy==1.26.4
ase==3.26.0
matscipy==1.1.1
```

Use installed Python 3.12.8 first. Future commands, only after B03 CPU release:

```sh
rtk proxy uv venv --python 3.12.8 --no-python-downloads target/b04-mace-venv
rtk proxy uv pip compile --python target/b04-mace-venv/bin/python --generate-hashes bench/requirements-b04-mace.in -o bench/requirements-b04-mace-py312.lock
rtk proxy uv pip sync --python target/b04-mace-venv/bin/python --require-hashes bench/requirements-b04-mace-py312.lock
```

Record Python exact version/platform, wheel hashes and MACE release source identity. Failed Python 3.12 compatibility leaves G2 blocked for environment input. Prepare a separate exact Python 3.11 bootstrap manifest and revise the environment contract before downloading/installing it; no implicit runtime/framework substitution. Do not mutate B01's venv or install into global Python.

- [ ] Add ignored test `official_mace_initialized_graph_energy_forces` that errors if MODENA_MACE_PYTHON is absent. Run red before the helper exists:

```sh
rtk proxy env MODENA_MACE_PYTHON=target/b04-mace-venv/bin/python cargo test --locked -p modena-core --test periodic_graph_official official_mace_initialized_graph_energy_forces -- --ignored --nocapture
```

The Rust driver should canonicalize an absolute helper/interpreter path from the worktree root so cargo's crate working directory cannot break it. Expected red is missing parity implementation, not environment or import failure; establish imports/constructor compatibility separately before treating the test cycle as meaningful.

- [ ] Implement the exact model constructor from Spec (seed=20261106, r_max=3.1, two elements [6,8], two interactions, hidden 8x0e+8x1o, correlation=2, CPU float64, one thread). Use official Configuration/AtomicData/Batch and same initialized model object on fresh cloned input tensors. Alternate batch changes exactly graph fields, never atom/node data:

```python
from mace.data import AtomicData, Configuration
from mace.tools import AtomicNumberTable
from mace.tools.torch_geometric.batch import Batch
config = Configuration(atomic_numbers=np.asarray(f['atomic_numbers']),
    positions=np.asarray(f['positions']), properties={}, property_weights={},
    cell=np.asarray(f['cell']), pbc=(True,True,True))
reference = Batch.from_data_list([AtomicData.from_config(
    config, z_table=AtomicNumberTable([6,8]), cutoff=3.1)])
alternate = reference.clone()
edge_array = [[e['source'],e['destination']] for e in request['edges']]
alternate.edge_index = torch.tensor(edge_array, dtype=torch.long).T.contiguous()
alternate.unit_shifts = torch.tensor([e['edge_cell_shift'] for e in request['edges']],
                                    dtype=torch.float64)
alternate.shifts = alternate.unit_shifts @ torch.tensor(f['cell'], dtype=torch.float64)
ref_out = model(reference.clone().to_dict(), compute_force=True, training=False)
alt_out = model(alternate.clone().to_dict(), compute_force=True, training=False)
torch.testing.assert_close(alt_out['energy'],ref_out['energy'],rtol=1e-10,atol=1e-10)
torch.testing.assert_close(alt_out['forces'],ref_out['forces'],rtol=1e-10,atol=1e-10)
```

Before forward, assert exact tuple equality and identical non-graph fields. Require exact integral unit_shift values before converting their tensor dtype. Model recipe has nonempty edges; preserve explicit empty 2x0/0x3 shapes if helper becomes general. Do not use pretrained calculators or torch.no_grad.

- [ ] Implement state_digest using sorted state_dict keys, encoded name/dtype/shape and contiguous CPU tensor bytes. Emit actual energy/max-component force errors; compare official original edge order and canonical order. For each fixed mutation (drop one nonzero image, drop all self images, change one shift), evaluate through actual MACE and require suite sensitivity >1e-10. Reject an ignored-graph or zero-output path. Central differences at h=1e-5 rebuild actual graphs for displaced positions and require force atol=rtol=1e-6, with unchanged neighbor sets/clear margin. Repeat calls to verify stable output.
- [ ] Run the named gate green, preserve exact manifest/model digest, and commit only small source/lock/evidence:

```sh
rtk proxy env MODENA_MACE_PYTHON=target/b04-mace-venv/bin/python cargo test --locked -p modena-core --test periodic_graph_official official_mace_initialized_graph_energy_forces -- --ignored --nocapture
rtk git add bench/requirements-b04-mace.in bench/requirements-b04-mace-py312.lock crates/modena-core/tests/support/mace_parity.py crates/modena-core/tests/periodic_graph_official.rs docs/B04_VERIFICATION.md
rtk git commit -m ':white_check_mark: Add: official initialized MACE graph parity'
```

**Gate:** actual model code, official graph, actual Rust graph, autograd E/F and mutation controls must all run. Success means software parity for this MACE instance only. Missing framework/runtime remains pending G2; do not replace it with EMT/LJ/toy tests.

## Task 5: Extract the scalar triclinic path and preserve bounded/skin behavior

**Files:** Modify `src/neighbor.rs`, create `src/neighbor/triclinic.rs`, modify core Cargo.toml and feature-gated exports in `src/lib.rs`, create `tests/triclinic_neighbors.rs`. No stencil default change.

**Consumes:** actual post-B03 NeighborList and pair-scaling behavior. **Produces:** retained scalar oracle plus diagnostic strategy selection and exact regression baseline.

Internal interfaces:

```rust
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(super) enum TriclinicBuildStrategy { Auto, Scalar, FractionalStencil }
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(super) enum FallbackReason {
    UncertifiedGeometry, FractionalBoundary, ImageDomain, GridCap,
    StencilCap, WorkCap, MemoryCap, AllBins, CandidateCost,
}
pub(super) struct BuildReport {
    pub(super) selected: TriclinicBuildStrategy,
    pub(super) fallback: Option<FallbackReason>,
    pub(super) candidate_checks: u64,
    pub(super) requested_workspace_bytes: usize,
}
// pair cap and zero/zero exclusions are passed through; emission remains fallible.
pub(super) fn build_triclinic<const BOUNDED: bool>(system: &System, list_cutoff2: f64,
    exclusions: &PairScaleSet, max_pairs: usize,
    strategy: TriclinicBuildStrategy, out: &mut Vec<(usize,usize)>)
    -> crate::Result<BuildReport>;
```

Keep actual module visibility compatible with PairScaleSet; if nested `neighbor::triclinic` sees private parent fields, avoid exposing unrelated fields. Graph errors map to fallback, not a new global physics error for ordinary classical requests.

- [ ] Add feature `neighbor-diagnostics=[]` to core. Under it expose enum `NeighborBuildStrategy { Auto, Scalar, FractionalStencil }`, chainable `with_triclinic_strategy_for_diagnostics`, and getter returning diagnostics. Normal API remains unchanged. Default internal Auto chooses Scalar until Task 8.
- [ ] Add red tests for exact pair ordering and bounded/zero-skin invariants. Use these names: `scalar_order_is_lexicographic`, `zero_skin_rebuilds_every_build`, `bounded_zero_skin_matches_scalar`, `bounded_pair_overflow_still_errors`, `same_bounding_box_changed_basis_rebuilds`, `periodic_crossing_skin_is_complete`, `pair_scale_zero_zero_only_exclusion`. Run:

```sh
rtk cargo test --locked -p modena-core --features neighbor-diagnostics --test triclinic_neighbors
```

- [ ] Move scalar i<j enumeration without changing membership/order/exclusions/cap semantics into the private module. Keep finalization/reference positions/full box/rebuild_count in one success path in NeighborList. A failed bounded build must retain its documented failure state and never publish a partial success. No O(N^2) staging vector is added.
- [ ] Verify targeted regressions including existing force_hardening/snapshot/pair_scales and default no-feature build:

```sh
rtk cargo test --locked -p modena-core --features neighbor-diagnostics --test triclinic_neighbors --test force_hardening --test snapshot --test pair_scales
rtk cargo test --locked -p modena-core --test force_hardening --test snapshot
```

- [ ] Commit the extraction only after exact scalar behavior passes:

```sh
rtk git add crates/modena-core/src/neighbor.rs crates/modena-core/src/neighbor/triclinic.rs crates/modena-core/src/lib.rs crates/modena-core/Cargo.toml crates/modena-core/tests/triclinic_neighbors.rs
rtk git commit -m ':refactor: Extract: scalar triclinic neighbor reference'
```

**Gate:** scalar E/F/W and pair ordering remain unchanged; bounded and zero-skin failures cannot be weakened.

## Task 6: Implement proof-checked fractional bins/stencil with explicit caps

**Files:** Extend `src/periodic_geometry.rs`, `src/neighbor/triclinic.rs`, `tests/triclinic_neighbors.rs`; use fixed cube oracle from Task 3. No change to System::min_image predicate.

**Consumes:** certified cell/error bounds, scalar path, list radius and caps. **Produces:** safe FractionalStencil candidate strategy, fallbacks and G3 evidence.

- [ ] Add tests `stencil_matches_scalar_on_all_fixed_cells`, `stencil_matches_fixed_cube_membership`, `wrapped_offsets_are_unique_for_dims_1_2_3`, `uncertain_fractional_face_falls_back`, `caps_fall_back_without_missing_pairs`, `stencil_preserves_classical_force_energy_virial`, `skin_half_threshold_and_box_change`, and `bounded_stencil_overflow_never_succeeds`. Force diagnostics to prove at least one non-boundary moderate fixture actually selects stencil; tests that all fall back do not prove stencil correctness. Construct interior atom fractions using dyadic offsets away from the chosen bin faces.
- [ ] Run red feature tests:

```sh
rtk cargo test --locked -p modena-core --features neighbor-diagnostics --test triclinic_neighbors
```

- [ ] Implement scalar_min_image_search_radius exactly from Spec, with post-B03 operation-count reconciliation. Add numeric unit tests for t_error, |t|<2^52 wrapping, J<=4096, and containment of scalar-selected displacements versus real-image counterparts in all fixed cube fixtures. Any uncertainty or overflow selects scalar with reason, never an incomplete stencil.
- [ ] Implement interval-certified single-bin wrapping, m=floor(1/w), modulo deduplicated offsets and conservative offset filter. Grid<=1_000_000, offsets<=4096, N*offsets<=50_000_000, workspace<=134_217_728, existing i32 bounds and candidate_estimate<.8*all_pairs are hard eligibility checks. Compute estimate from occupancies and deduplicated unordered bin pairs with checked arithmetic; do not use observed timing to alter it mid-campaign.

```rust
fn fallback_for_cost(candidate_estimate: u128, n: usize) -> bool {
    let n128 = n as u128;
    let all_pairs = match n128.checked_mul(n128.saturating_sub(1)) {
        Some(value) => value / 2,
        None => return true,
    };
    match candidate_estimate.checked_mul(5).zip(all_pairs.checked_mul(4)) {
        Some((candidate_scaled, scalar_scaled)) => candidate_scaled >= scalar_scaled,
        None => true,
    }
}
```

The saturating subtraction handles N=0 only; all products are checked and any overflow is fallback `candidate_cost`. Use a cap-preflighted sorted occupied-bin list or existing linked arrays; no atom enters multiple bins in this first version. Pair emission is i<j, exact min_image prune, exclusions, bounded push; sort pairs lexicographically before publication. Deduplicate bin targets, not accepted pairs after unbounded overscan.

- [ ] Run feature and standard regressions green. Require all optimized fixture pair sets equal both retained scalar and independent cube. E/F/W atol=rtol=1e-10, existing tighter tests unchanged. Inspect diagnostics for every cap reason and requested workspace. Auto must still choose scalar.
- [ ] Commit correctness implementation:

```sh
rtk cargo test --locked -p modena-core --features neighbor-diagnostics --test triclinic_neighbors --test force_hardening --test snapshot --test pair_scales
rtk git add crates/modena-core/src/periodic_geometry.rs crates/modena-core/src/neighbor/triclinic.rs crates/modena-core/tests/triclinic_neighbors.rs docs/B04_VERIFICATION.md
rtk git commit -m ':sparkles: Add: certified fractional triclinic neighbor stencil'
```

**Gate:** real/finite-precision coverage proofs, actual stencil-selected tests, pair ordering, E/F/W, bounded/skin behavior and cap fallbacks all pass. G3 does not imply G4.

## Task 7: Add the dedicated full-rebuild/full-step cost campaign

**Files:** Create benchmark bin/harness/test paths in file map; modify bench Cargo.toml to add opt-in `neighbor-diagnostics=["modena-core/neighbor-diagnostics"]` and a `[[bin]]` named `triclinic_neighbors` requiring that feature. Do not edit B03 campaign files.

**Consumes:** strategy API, exact same serialized snapshots/velocities and CPU lease release. **Produces:** repeatable raw paired cost records and offline comparison report. Public CLI is fixed:

```text
triclinic_neighbors --snapshot PATH --strategy scalar|fractional|auto|paired --candidate fractional|auto --mode rebuild|step
  --skin 0|0.3 --threads 1|4 --warmup-blocks 5 --sample-blocks 20 --block-steps 10
  --dt 0.0005 --output PATH
python bench/triclinic_neighbors.py --binary PATH --output-dir PATH --generate-only
python bench/triclinic_neighbors.py --binary PATH --output-dir PATH --run
python bench/triclinic_neighbors.py --output-dir PATH --analyze-only
```

JSON records include schema_version=1, base_rev, binary_sha256, snapshot_sha256, model/precision, host/compiler/threads, strategy/requested/selected/fallback, mode, skin, sample_index, order, block_steps, elapsed_ns, rebuild_count, candidate_checks, n_pairs, requested/retained_bytes, allocation_count, max_energy_error,max_force_error,virial_error. Timing schema rejects missing/nonfinite/negative values, snapshot/revision mismatch, strategy mixing and insufficient paired samples.

- [ ] Add meaningful parser/statistical tests before benchmark wiring: reject mismatched snapshot/revision, missing scalar counterpart, force parity failures, unavailable sample, thermal/noise inconclusiveness, and verify known paired ratios/seeded bootstrap on a fixed synthetic timing table. These are timing arithmetic tests, not scientific model substitutes.
- [ ] Implement deterministic snapshot generation. For atom i with integer grid indices x,y,z on side M=ceil(cuberoot(N)), use fractional coordinate `(axis_index+.31+jitter(i,k))/M`, where `jitter=.06*((i*1103515245+(k+1)*12345+20261106)%65536/65535-.5)`. Compute hash arithmetic in u64 with declared modulo, serialize exact resulting f64 positions. Initial velocity component uses `.05*(2*((i*1664525+(k+1)*1013904223+20261106)%65536)/65535-1)`, then subtract component mean in atom order. Store resulting vectors; no regeneration by strategy.

Use spec N/shape/skin/thread ladder and LJ epsilon=sigma=mass=1, force cutoff=2.5, no tail correction, VelocityVerlet dt=.0005. Call the same existing PairForceModel/Integrator construction APIs identified in Task 0; retain actual default cutoff shift treatment identically in both strategies and record it. If post-B03 APIs require additional physical parameters, fix/record them before the first numerical run, not between strategies.

- [ ] Implement paired execution in one release build with scalar/candidate controls differing only strategy. The harness invokes `--strategy paired --candidate fractional` before any selector exists and `--candidate auto` to validate a promoted selector. The binary owns both replicas, alternates AB/BA by sample_index, initializes both from the identical serialized state, and evolves each continuously between paired blocks. Warmup 5 blocks, 20 paired blocks, 10 steps each; rebuild mode times the entire forced build, step mode times complete repeated steps. Include cell changes as separate forced rebuild cases. Use monotonic clocks and retain raw nanoseconds; no generation/parsing inside timed kernels. Enforce per-case wall-time cap 600 seconds and total campaign cap 7200 seconds; cap aborts are incomplete evidence, not passes. Natural skin=.3 amortized claims require at least two measured rebuilds per replica; if absent, report that regime inconclusive without extending the protocol after results. Count all conversion/grid/sort/reduction/fallback costs. No timing of an unvalidated strategy.
- [ ] Run small harness tests and compile the named binary once:

```sh
rtk proxy python3 -m unittest discover -s bench/tests -p test_triclinic_neighbors.py
rtk cargo build --locked --release -p modena-bench --features neighbor-diagnostics --bin triclinic_neighbors
```

Do not start performance trials if any of G1-G3 is pending or the CPU lease is unavailable. Generate-only can prepare input bytes outside the timed run once authorized.

- [ ] Commit the harness and predeclared manifest contract:

```sh
rtk git add crates/modena-bench/Cargo.toml crates/modena-bench/src/bin/triclinic_neighbors.rs bench/triclinic_neighbors.py bench/tests/test_triclinic_neighbors.py docs/B04_VERIFICATION.md
rtk git commit -m ':white_check_mark: Add: paired triclinic neighbor cost campaign'
```

**Gate:** campaign measures full rebuilding and complete steps from the same inputs; no graph-only result can satisfy default selection. Thermal/noise/incomplete data is unresolved evidence.

## Task 8: Measure, decide selection and publish limited acceptance evidence

**Files:** Small result manifests under `docs/B04_VERIFICATION.md`, user/API docs `docs/B04_PERIODIC_GRAPHS.md`; selection constants only in `src/neighbor/triclinic.rs` after passing G4. Generated raw samples stay in ignored `target/b04-triclinic-campaign` with hashes and stable archive location referenced by the small report.

**Consumes:** passed G1-G3, preregistered manifests, released CPU lease, one release binary. **Produces:** explicit G4 decision and honest B04 acceptance scope.

- [ ] Execute the frozen campaign once, preserving all samples and failures:

```sh
rtk proxy python3 bench/triclinic_neighbors.py --binary target/release/triclinic_neighbors --output-dir target/b04-triclinic-campaign --generate-only
rtk proxy python3 bench/triclinic_neighbors.py --binary target/release/triclinic_neighbors --output-dir target/b04-triclinic-campaign --run
rtk proxy python3 bench/triclinic_neighbors.py --output-dir target/b04-triclinic-campaign --analyze-only
```

Analyze median paired optimized/scalar ratio and seeded 10_000-resample 95% bootstrap interval. Moderate-tilt N>=2048 cases at every enabled thread/skin setting require median<=.95 and interval upper<1; small/fallback cases require <=5% median overhead. Retain full-rebuild results, memory/counters and cases that lose.

- [ ] If criteria fail, leave Auto=Scalar and record G4 failed/pending with actual reasons. Do not manufacture acceptance by narrowing the campaign after results. If criteria pass, add a selector restricted to measured N>=2048, shape/cost eligibility and measured hardware; record exact host/platform predicate and require scalar on unknown hosts. Re-run selector decisions against saved case manifests and small correctness regressions. Changing the selector or compiled production branch requires a fresh performance validation of that exact binary; do not reuse pre-selector timings as if identical.
- [ ] Write API docs showing borrowed graph construction, integer convention, all-image/self policy, budgets/errors, graph-only partial PBC and rebuilding. Write verification with every actual command/revision/hash/tolerance/outcome; label initialized MACE as software parity, preserve FAIR-Chem/B10 separation. No global statement that all FAIR-Chem/MACE/NequIP graph contracts match.
- [ ] Run required repository checks once after final changes, including full core suite, meaningful graph/triclinic regressions, explicit official gates and harness tests. Use the commands actually required by the post-B03 repository; if an unrelated baseline failure occurs, record it rather than mislabel the change complete. Do not repeatedly broaden tests after all required checks pass without a new change/concern.
- [ ] Commit the measured decision/evidence:

```sh
rtk git add crates/modena-core/src/neighbor/triclinic.rs docs/B04_PERIODIC_GRAPHS.md docs/B04_VERIFICATION.md docs/superpowers/plans/2026-11-06-b04-periodic-neighbors.md
rtk git commit -m ':memo: Document: B04 graph parity and triclinic cost decision'
```

**Final gate:** report G1/G2/G3/G4 individually. B04 optimization is not accepted when its required performance gate is unresolved. A correct graph with unprofitable stencil may remain useful, with scalar default and incomplete broader optimization acceptance clearly stated. No implementation or scientific acceptance is claimed by these draft documents.

## Self-review and remaining execution gates

- Geometry formula, B-independent global exclusion/contraction proof, finite-precision predicate/local envelope, independent cube certificates, self/multi-image/partial-PBC semantics, per-pass caps and atomic failures each map to Tasks 1-3.
- Actual pinned official MACE graph consumption/autograd/mutation sensitivity maps to Task 4; model/runtime compatibility remains a real execution gate, not an omitted design choice.
- Scalar ordering, B09/snapshot bounded behavior, zero skin and proof-checked stencil/cost fallbacks map to Tasks 5-6.
- Same-snapshot complete-step timings, precision freezes and reproducible-benefit selector map to Tasks 7-8.
- Initial source reconciliation is against reviewed B03 core at 0c54db3; synchronize the unpublished B04 base with final reviewed B03 before publication. The separate B03 checkout is preserved and its released timing lease is not consumed.

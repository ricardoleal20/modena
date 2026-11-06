# B04 periodic graphs and triclinic neighbors — draft design contract

Status: implementation contract adopted after controller design self-review; no implementation or acceptance is claimed. The B-independent global shift certificate is mandatory. B04 is isolated on `ricardo/b04-periodic-neighbors`, initially based on reviewed B03 source `0c54db355b830c49e687e97e2e9aa8502bbaa5d1`. Core tree `39b46c7496c6820219056fe38a821378a32633e8` is identical to reviewed prototype `66af0e4`; no core/CLI source changed in the B03 campaign. B03's measurement CPU lease has been explicitly released. Independent numerical work may start while B03 closes its reports/reviews; synchronize the unpublished B04 base with final reviewed B03 before publication, without force-push or an unauthorized merge. Companion plan: `docs/superpowers/plans/2026-11-06-b04-periodic-neighbors.md`.

## Goal and acceptance boundaries

Deliver a directed, integer-image periodic graph API; independent graph and official-model software parity; a separately proved fractional triclinic classical half-list implementation; and a reproducible complete-step measurement that controls default selection. These are gates G1 geometry, G2 official initialized MACE software parity, G3 classical correctness, and G4 full-step cost. All four outcomes must be recorded explicitly. No skipped external gate counts as passed.

G2 uses an actual official MACE architecture with initialized fixed weights. Its result covers that release, model instance, declared snapshots and arithmetic only. It does not establish chemical efficacy, a trained deployment, FAIR-Chem graph equivalence, FAIR-Chem checkpoint/task parity, or B10 acceptance. The B10 actual FAIR-Chem task/weights/material/reference gate remains separate.

Partial PBC belongs only to the graph's geometry interface. Runner, B01/B09 provider validation, MD wrapping, electrostatics and advertised MD PBC remain unchanged. Classical minimum-image half pairs are a different interaction contract from an all-image periodic graph. Do not interchange them.

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

## Public graph API

Create `crates/modena-core/src/periodic_graph.rs`, private `periodic_geometry.rs`, and exports in `lib.rs`. Use existing thiserror; no production dependency addition is required.

```rust
#[derive(Clone, Copy, Debug, Eq, PartialEq, Ord, PartialOrd)]
pub struct PeriodicEdge {
    pub source: usize,
    pub destination: usize,
    pub edge_cell_shift: [i64; 3],
}
#[derive(Clone, Copy, Debug)]
pub struct GraphBudget {
    pub max_atom_pairs: u64,
    pub max_image_candidates: u64,
    pub max_edges: usize,
    pub max_bytes: usize,
    pub max_abs_shift: i64,
}
#[derive(Clone, Copy, Debug)]
pub struct PeriodicGraphOptions {
    pub cutoff: f64,
    pub pbc: [bool; 3],
    pub include_periodic_self_images: bool,
    pub budget: GraphBudget,
}
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct GraphBuildStats {
    pub atom_pair_visits: u64,
    pub image_candidate_visits: u64,
    pub edges: usize,
    pub requested_bytes: usize,
}
#[derive(Debug, thiserror::Error)]
pub enum GraphBuildError {
    #[error("invalid periodic graph input: {0}")]
    InvalidInput(String),
    #[error("uncertifiable periodic geometry: {0}")]
    NumericalGeometry(String),
    #[error("periodic graph image shift outside supported range")]
    ShiftRange,
    #[error("periodic graph budget {resource} exceeded: {limit}")]
    BudgetExceeded { resource: &'static str, limit: u64 },
    #[error("periodic graph allocation failed")]
    Allocation,
}
pub struct PeriodicGraph<'a> {
    system: &'a System,
    options: PeriodicGraphOptions,
    edges: Vec<PeriodicEdge>,
    stats: GraphBuildStats,
}
impl<'a> PeriodicGraph<'a> {
    pub fn build(system: &'a System, options: PeriodicGraphOptions)
        -> Result<Self, GraphBuildError>;
    pub fn edges(&self) -> &[PeriodicEdge];
    pub fn system(&self) -> &'a System;
    pub fn options(&self) -> PeriodicGraphOptions;
    pub fn stats(&self) -> GraphBuildStats;
}
```

`System` in this interface is the existing type; `Result` is the two-argument standard result, not the crate's single-parameter alias. Borrowing the snapshot prevents position/identity/order/cell mutation while the graph is used. No cache is added. End the borrow and rebuild after any change in snapshot, PBC, cutoff or self-image policy. System atom identity is retained by index; do not infer atomic numbers from force-field species. Geometry requires `positions.len()==3*n`; other physical arrays are not read. Chemical validation remains the consuming model's responsibility.

With cell columns a,b,c, edge `(i,j,S)` has `d=r_j-r_i+a*Sx+b*Sy+c*Sz`. Nonperiodic lattice axis k permits only S[k]=0. Always omit `(i,i,0)`; include other self images iff requested. Distinct overlapping atoms remain geometric edges. Sort by source, destination, shift x/y/z. Every emitted edge has its reversed tuple `(j,i,-S)`. Translation `r_i'=r_i+H*T_i` transforms shifts by `S'=S+T_i-T_j`.

The scalar expression is fixed per component k:

```rust
let x0 = rj[k] - ri[k];
let x1 = x0 + basis[0].0[k] * shift[0] as f64;
let x2 = x1 + basis[1].0[k] * shift[1] as f64;
let d_k = x2 + basis[2].0[k] * shift[2] as f64;
// Compute the three components with the same order.
let n2 = (d[0] * d[0] + d[1] * d[1]) + d[2] * d[2];
let inside = n2 < cutoff * cutoff;
```

Reject nonfinite intermediate evaluations rather than emit a nonfinite vector. This is an explicit f64 predicate; no assertion of exact-real classification at every possible boundary is made. Analytically representable exact equality is excluded. Official parity fixtures must be at least 1e-8 Angstrom from cutoff; adjacent-f64/equality tests are separate.

## Complete enumeration with a certified numerical envelope

Revision 2: the global range certificate below closes the earlier circular dependence on configured max_abs_shift. A local rounding envelope valid only within a user-selected range is insufficient to establish global completeness. This revision narrows the supported numerical domain explicitly and rejects uncertifiable requests before enumeration; it does not truncate images or change the physical cutoff.

Let H be the exact-real matrix whose entries are the supplied f64 cell entries; U=H^-1 and delta=rj-ri in exact arithmetic. For any image with real norm <=Q, `|U_row[k]*delta+S[k]|<=Q*norm(U_row[k])`. Therefore the complete integer search has bounds `ceil(-s[k]-w[k])..floor(-s[k]+w[k])`, with nonperiodic axes fixed at zero. This proof supports unreduced/skew cells, multiple images and self images. It does not justify rounding fractional displacement or a constant ±1 image cube.

Implement a private `Interval { lo:f64, hi:f64 }` and `CertifiedCell { basis:[Vec3;3], inverse:[Vec3;3], inverse_error_inf:f64, inverse_row_norm_upper:[f64;3] }`. Exact point inputs use singleton intervals. Each operation encloses the real operation with `next_down` and `next_up` after correctly rounded +,-,*,/,sqrt; multiplication/division use min/max of all four endpoints; division rejects denominators containing zero. Overflow/nonfinite intervals fail certification; underflow is enclosed by neighboring representable values. Square/norm helpers explicitly handle intervals containing zero. Do not use subnormal-flushing arithmetic.

`certified_cell(SimBox)->Result<CertifiedCell,GeometryError>` performs:

1. Validate finite basis, positive finite volume and finite approximate inverse A from `inv_basis_rows()`. Negative orientation is allowed through positive absolute volume.
2. Compute interval E=I-AH, and outward upper rho for its maximum absolute row sum. Require rho<0.5; larger residual returns NumericalGeometry, even if the cell could mathematically be invertible.
3. Compute outward `eta=rho/(1-rho)*norm_inf(A)`. Since U=(I-E)^-1 A, `norm_inf(U-A)<=eta`. Bound row norm of U by outward `norm_2(A_row)+eta`.
4. Enclose exact delta by interval subtraction. Enclose A_row*delta and inflate by outward `eta*norm_inf(delta)` to obtain `[s_lo,s_hi]`.

Completeness must also cover the specified f64 membership predicate, whose rounded d can pass slightly outside the exact-real R sphere. Define EPS=f64::EPSILON, T=f64::from_bits(1), gamma_m=m*EPS/(1-m*EPS), all upper calculations outward. Compute `R2=R*R` once and require finite positive **normal** R2 (`R2.is_normal()`): zero/subnormal underflow and infinity/NaN overflow are InvalidInput, including on empty inputs. This same unchanged R2 is used by final membership.

### Global shift/range certificate, independent of configured B

This certificate is required **before** using any B-dependent local envelope. It covers all integers in the periodic lattice for exact-real membership, and all i64 shifts for the stated f64 expression, including shifts whose conversion to f64 is inexact. Let P be the enabled lattice axes. Shifts on other axes are zero. Let K=norm_inf(S) in real arithmetic and define outward upper bounds

```
U_P = max_{k in P} (sum_l abs(A[k,l]) + eta)
H_P = max_l sum_{k in P} abs(H[l,k])
p = norm_inf(ri) + norm_inf(rj)
D = norm_inf(exact_delta_interval)
Q_fp = sqrt((R2 + 5*T)/(1-gamma_5))
Q_0 = max(R, Q_fp)
alpha = gamma_10 * U_P * H_P
C = U_P * (D + Q_0 + gamma_10*p + 10*T)
K_global = C/(1-alpha)
G = ceil(K_global)
```

For no enabled periodic axis set U_P=H_P=K_global=G=0. Otherwise require finite upper alpha<0.5 and finite bound arithmetic; failure is NumericalGeometry. All displayed quotients use full interval evaluation: in particular 1-alpha and 1-gamma_m use lower denominator enclosures before the result is rounded upward. Upper-rounding a denominator before division would not certify the upper bound and is forbidden. The eta matrix-infinity error bounds each inverse-row 1-norm error, so U_P is a certified upper bound even for an unreduced cell. None of U_P, H_P, p, D, Q_0, alpha, C, K_global depends on configured max_abs_shift.

**Proof for the f64 predicate over the entire i64 domain.** Let x=delta+H*S be the exact-real displacement for an arbitrary integer i64 shift S with disabled components zero. The cast of each i64 component to f64 is finite and satisfies `abs(fl(S)-S)<=EPS*abs(S)`, including |S|>2^52. Count the cast as an additional rounding operation in each multiplication branch. One subtraction, three casts, three products and three additions are conservatively bounded by gamma_10:

```
norm_inf(d_computed-x)
    <= gamma_10*(p + H_P*K) + 10*T.
```

The additive underflow allowance is conservative under gradual underflow; conversions of nonzero integers cannot underflow. Any accepted finite squared-norm result implies all the specified component/product/sum intermediates were finite: in this expression an infinity or NaN cannot later become a finite passing norm by addition or nonnegative squared summation. Thus the finite-operation bound applies to every accepted i64 tuple, including hypothetical tuples outside B. Positive squared-sum rounding gives `norm_inf(d_computed)<=norm_2(d_computed)<=Q_fp<=Q_0`. Applying the enabled rows of S=U*(x-delta) yields

```
K <= U_P*(D + Q_0 + gamma_10*p + 10*T) + alpha*K
  = C + alpha*K,
```

and alpha<1 implies K<=K_global. This proves exclusion globally; it does not assume |S|<=B. The stricter alpha<0.5 is the supported-domain policy.

**Proof for physical images beyond the i64 domain.** For any integer lattice shift, with no machine representation assumption, exact-real norm(x)<R gives `K<=U_P*(D+R)<=K_global`. Consequently the same certificate bounds every physically required image, including one whose hypothetical index would exceed i64. There is no silent omission of such images: unsupported requests are rejected below.

Before converting G to any integer, require G finite, nonnegative and `G<=B<=2^52`, where B is configured max_abs_shift. The right endpoint 2^52 is exactly representable and lies strictly inside i64. If G exceeds B or 2^52, return ShiftRange for the entire request; do not enumerate only the supported subset. If bound arithmetic itself cannot be certified, return NumericalGeometry. Since every physically required image and every i64 tuple accepted by the f64 predicate lies within G, accepted requests have no relevant image outside B or the exact-conversion domain. Only after this global proof may the tighter existing B-dependent reciprocal search run. Conservative G can reject a request whose actual neighbors would fit; that explicit supported-domain restriction is preferable to truncation.

Huge unwrapped coordinates are handled by the same certificate. They can make D or p large, overflow a bound, make G>B, or lose a useful range certificate; these produce explicit errors. A huge common translation can therefore be rejected despite unchanged real geometry. Do not silently wrap/recenter input or coerce shifts, because that changes the requested integer-image gauge. A caller may submit a separately recentered snapshot and accept its changed gauge. The fixed gauge fixtures remain supported because their certified bounds are far below fixture B=1_000_000.

### Local reciprocal bounds after global certification

Let B be configured max_abs_shift, now known to contain every required and f64-accepted image. Per Cartesian component use outward

```
e[k] = gamma_7*(abs(ri[k])+abs(rj[k])+B*sum_l abs(H[k,l])) + 7*T
Q_fp = sqrt((R2+5*T)/(1-gamma_5))
Q_enum = Q_fp + norm_2(e)
w[k] = Q_enum * inverse_row_norm_upper[k]
lo[k] = ceil(lower(-[s_lo,s_hi]-w[k]))
hi[k] = floor(upper(-[s_lo,s_hi]+w[k]))
```

The first bound encloses the seven scalar arithmetic operations producing d_k; positive sum-of-squares rounding supplies Q_fp. Q_enum is a search-only radius. Final membership still uses R2. All formula operations are interval/outward, including gamma denominators; if uncertain, choose larger operation counts, never smaller. Exact cells/boundary tests and an operation-count review must verify this implementation before G1. This is a concrete sufficient envelope; conservative over-enumeration is acceptable within budgets.

Cutoff must be finite positive with finite normal R2 as checked before any pair work. Require `0<=B<=2^52`; the independent global certificate must prove G<=B before local bounds. Local endpoints outside B or outside an exactly representable integer range return ShiftRange before casting. Either conservative global or local bounds may fail even if the actual accepted images would fit. Require every integer endpoint/range length to be representable; preflight span products in u128. Do not saturating-cast f64 endpoints.

Private geometry interfaces consumed by graph/stencil tests:

```rust
pub(crate) enum GeometryError { Invalid, Uncertified, ShiftRange }
pub(crate) struct GeometryOptions {
    pub(crate) cutoff: f64,
    pub(crate) pbc: [bool;3],
    pub(crate) max_abs_shift: i64,
}
fn certified_cell(cell: SimBox) -> Result<CertifiedCell, GeometryError>;
fn global_shift_bound(cell: &CertifiedCell, ri: Vec3, rj: Vec3,
    options: GeometryOptions) -> Result<i64, GeometryError>;
fn image_bounds(cell: &CertifiedCell, ri: Vec3, rj: Vec3,
    options: GeometryOptions) -> Result<[(i64,i64);3], GeometryError>;
fn edge_displacement(cell: &CertifiedCell, ri: Vec3, rj: Vec3,
    shift: [i64;3]) -> Result<Vec3, GeometryError>;
fn strict_radius(d: Vec3, cutoff2: f64) -> Result<bool, GeometryError>;
```

global_shift_bound computes the B-independent global envelope, verifies G<=B<=2^52 before conversion, and returns its exactly representable integer G. image_bounds must call it before computing a local B-dependent Q_enum; callers cannot bypass the certificate. Map GeometryError::Invalid to GraphBuildError::InvalidInput, Uncertified to NumericalGeometry, and ShiftRange unchanged at the public graph boundary. This keeps private geometry independent of public graph types. Use crate-private visibility so unit tests can exercise certificates. A future exact-real cutoff policy would require a separate exact/adaptive predicate and a new contract; do not silently substitute it here.

## Resource policy and atomic publication

No unbounded default constructor. Fixture driver budget is exactly `(max_atom_pairs=1_000_000, max_image_candidates=10_000_000, max_edges=1_000_000, max_bytes=134_217_728, max_abs_shift=1_000_000)`.

Two enumeration passes are intentional. `max_atom_pairs` counts both visits, so preflight `2*N*N` before any pair loop; `max_image_candidates` counts both passes, so each first-pass span product is charged twice before it is enumerated. Count candidates even if zero-self, self policy or radius rejects them. Nonperiodic singleton axes count one. Empty spans count zero. `max_edges` counts final accepted edges. All counters use checked arithmetic. Error resources are exactly `atom_pairs`, `image_candidates`, `edges`, `bytes`. Error precedence: malformed options/input including unusable R2, cell/numerical certification, base-pair budget, per-pair global certificate (NumericalGeometry or ShiftRange) then local shift ranges, candidate budget, edge budget, byte budget, allocation.

Pass 1 uses O(1) geometry workspace and counts edges, checking the running budgets; pass 2 allocates exactly the accepted edge count with `try_reserve_exact`, emits, `sort_unstable`s and returns. Account requested edge-capacity bytes plus `size_of::<PeriodicGraph>()` and `size_of::<CertifiedCell>()` in max_bytes; no hidden candidate/interval vectors. Cap is logical requested storage, not OS allocator RSS. Check actual Vec capacity against budget before publication. Numeric and allocation failures return Err; no graph escapes, no System/forces/virial mutation occurs, and an existing independent graph remains intact. Add an allocation-failure seam under unit-test cfg only. Validate options/cell even for N=0; a valid empty graph returns zero pair/candidate/edge counters and fixed requested header/workspace bytes.

## Independent fixtures and G1

Create `tests/fixtures/periodic_graph/cases.json` and README under modena-core. Fixed cube is `[-16,16]^3` for every listed fixture and every PBC subset. The oracle does direct Cartesian enumeration of all ordered pairs and this fixed cube; it must not call production bounds, inverse helpers or min_image. Displacements/membership follow the declared scalar order. Filter nonperiodic axes and self policy. Sorted tuple equality is exact.

Fixture certificates use explicit triangular equations, independently of production reciprocal bounds. For H with a=(A,0,0),b=(B,C,0),c=(D,E,F) and bounded coordinates, accepted images imply `|Sz|<=(R+|delta_z|)/|F|`, then `|Sy|<=(R+|delta_y|+|E*Sz|)/|C|`, then `|Sx|<=(R+|delta_x|+|B*Sy|+|D*Sz|)/|A|`. Substitute fixture constants and round the displayed upper bounds outward; all are <16. Gauge translation adds at most (3,2,0) to these bounds. These fixed numeric certificates go in the fixture README, rather than a dynamic inverse-derived cube.

| ID | Exact fixture data / assertions |
| --- | --- |
| empty | cell diag(2,2,2), positions=[], R=1.5; empty valid graph, invalid options still fail. |
| self6 | same cell, one atom at origin, R=2.1: exactly six axial self shifts when enabled; none when disabled. With R=1.9 none. |
| multi4 | same cell, positions origin,(.6,0,0), R=1.5: four tuples (0,1,-ex),(0,1,0),(1,0,0),(1,0,+ex). |
| unreduced | a=(10,0,0),b=(49,10,0),c=(0,0,40), positions origin,(19.7,.1,0), R=4.8: contains (0,1,(-2,0,0)); retain prior force_hardening regression. |
| oblique | a=(10,0,0),b=(4.9,10,0),c=(0,0,40), positions origin,(6.37,3,0), R=4.8: contains shift (-1,0,0). |
| boundary | diag(8,8,8), origin and (2,0,0), R=2, plus next_down(2),next_up(2): equality/outside excluded, inside included. These do not enter official parity. |
| masks | multi4 and mace3 repeated with all eight explicit PBC masks; disabled axes' shifts exactly zero. |
| mace3 | a=(2.4,0,0),b=(.7,2.3,0),c=(.2,.3,2.5), numbers=[6,8,6], positions=(.1,.2,.3),(1,.7,.9),(.6,1.5,1.8), R=3.1, full PBC, periodic self enabled. |
| gauge | mace3 with atom 1 moved by 3a-2b; compare gauge-transformed tuples, vectors and forces. |
| deform | mace3 with b=(.8,2.35,0), same positions. |
| displaced | mace3 with atom 2 x increased by .037. |

Each official fixture preflight must establish no image is within 1e-8 of cutoff using the fixed cube before any official/model comparison. A failed margin is a fixture-contract failure, not permission to shift the cutoff after seeing model errors. In separate snapshot borrows test changes of positions, N, index ordering, atomic_numbers, cutoff, PBC and self policy. Equal-bounding-box alternate basis must rebuild the classical list. Invalid cases cover NaN/Inf coordinates/basis/cutoff, nonpositive and unusable cutoff square, singular cell, uncertifiable inverse, malformed positions length, endpoint/product overflow, a well-formed cell failing alpha<.5, huge coordinates with G>B, no-periodic-axis G=0, all i64 conversion boundary proof cases, each zero/tiny budget, and allocation failure. Atom/edge indices must never wrap.

ASE oracle uses installed interpreter `/Users/ricardoleal20/Documents/Code/Misc/modena/target/b01-ase-venv/bin/python`, observed ASE=3.29.0, NumPy=2.5.4, Python 3.14. It calls `neighbor_list('ijS', atoms, R, self_interaction=False)`, removes i==j iff periodic self is disabled, and canonicalizes exact tuples. Wire cell rows are [a,b,c], the transpose convention of column matrix H. Record package versions and neighborlist.py SHA256. Do not compare adjacency without integer shifts.

## Actual official model software parity, G2

Environment revision before any model evaluation (2026-11-06): the initial
candidate NumPy2.2.6 conflicted with matscipy1.1.1's declared runtime
`numpy>=1.16.0,<2.0.0`; original failed requirements/resolver logs are preserved.
Revise only NumPy to1.26.4; retain Python3.12.8 and all other direct pins, model
constructor/seed/fixtures/tolerances. Exact candidate resolution produced44
packages; hash-required dry sync and binary availability checks succeeded for
all six direct dependencies on this macOS arm64 host. python-hostlist2.3.0
requires its published source distribution; retain its hash and actual build/
installed-package provenance. Resolution is not import/model acceptance.
Primary metadata: [matscipy1.1.1 contract](https://raw.githubusercontent.com/libAtoms/matscipy/v1.1.1/pyproject.toml),
[MACE0.3.13 package metadata](https://pypi.org/pypi/mace-torch/0.3.13/json),
[NumPy1.26.4 supported Python](https://numpy.org/doc/stable/release/1.26.4-notes.html).
No numerical result exists under the superseded NumPy2 candidate; future G2
claims bind only to the actual revised hashed environment.


Preferred future environment is isolated `target/b04-mace-venv` using installed **Python 3.12.8**. Root inventory also found 3.14.5; neither known ASE venv contains MACE/Torch. Candidate pins: torch=2.6.0, mace-torch=0.3.13, e3nn=0.4.4, numpy=1.26.4, ase=3.26.0, matscipy=1.1.1. Resolve all transitive packages with hashes under 3.12.8 and record source release SHA. Compatibility is a gate, not a present assertion. Use CPU float64, one thread. If 3.12.8 resolution/import/actual forward fails due to supported-runtime incompatibility, stop G2; prepare a separate pinned Python 3.11 bootstrap manifest and revise the environment contract before installing it. Do not reuse/mutate B01 or silently change framework versions.

Pin the official model constructor below before evaluation; initialized parameters are not trained chemical weights. Record sorted state_dict tensor names/dtypes/shapes/raw CPU bytes SHA256, parameter count, package lock and platform.

```python
torch.set_default_dtype(torch.float64)
torch.set_num_threads(1)
torch.manual_seed(20261106)
np.random.seed(20261106)
model = modules.MACE(
    r_max=3.1, num_bessel=4, num_polynomial_cutoff=5, max_ell=1,
    interaction_cls=modules.RealAgnosticResidualInteractionBlock,
    interaction_cls_first=modules.RealAgnosticInteractionBlock,
    num_interactions=2, num_elements=2,
    hidden_irreps=o3.Irreps('8x0e + 8x1o'),
    MLP_irreps=o3.Irreps('8x0e'), atomic_energies=np.zeros(2),
    avg_num_neighbors=6.0, atomic_numbers=[6,8], correlation=2,
    gate=torch.nn.functional.silu, radial_MLP=[8,8],
).double().cpu().eval()
```

Use official `Configuration`, `AtomicData.from_config(...,AtomicNumberTable([6,8]),cutoff=3.1)` and official vendored `Batch.from_data_list` as reference. Build the alternate by cloning that batch and replacing only edge_index (torch.long 2xE), unit_shifts (float64 exact integers), and shifts=unit_shifts@cell from actual Rust serialized edges. Preserve node attributes, positions, cell, ptr, batch and head. Run fresh cloned positions through the same actual model with gradients enabled, `compute_force=True,training=False`. Never use no_grad. Compare original official ordering and canonical ordering; no accidental second official graph may stand in for Rust output.

mace3/gauge/deform/displaced must first have exact tuple parity, then finite/nontrivial E and every force component within atol=rtol=1e-10. Verify repeated call stability. Central differences h=1e-5 on displaced fixture with rebuilt graphs require force atol=rtol=1e-6 and no boundary crossing. Mutation controls drop a nonzero image edge, remove all periodic self edges, and change one nonzero shift. Across the fixed multi-atom suite each applicable mutation must change E or F by >1e-10; otherwise fail fixture adequacy. Do not train or select a new seed to mask a failure.

`MODENA_MACE_PYTHON` selects the interpreter. Dedicated official-gate tests are explicitly ignored in the minimal Rust suite and run by named `--ignored` commands; their driver errors if environment variables/packages are absent. The acceptance manifest lists pending/failed explicitly. No pretrained download, model hub access, GPU, CUDA or MPS is required for this software gate. Pinned MACE changes aperiodic cell rows, so the model gate uses full PBC only; geometry masks remain G1.

## Triclinic fractional stencil and G3

Create private `src/neighbor/triclinic.rs`. Retain the exact scalar i<j path as oracle; preserve PairScaleSet zero/zero exclusion, list radius=cutoff+skin, final min_image predicate, zero-skin rebuilding and bounded max_pairs behavior. Emit lexicographically sorted half pairs so serial triclinic accumulation order equals the existing scalar loop. No graph-image summation enters the classical force kernel.

For true inverse rows u and a search radius Q that encloses every scalar-accepted pair, let w[k]=Q*norm(u[k]). Fractional wrapped coordinates are binned into m[k]=max(1,floor(1/w[k])) subdivisions, capped below. A necessary bin-offset condition is `dist(q[k]/m[k],integers)<=w[k]+1/m[k]`. Enumerate integer q from `-ceil(m[k]*w[k])-1` through its positive endpoint, modulo m[k], deduplicate wrapped targets, retain every offset satisfying the conservative condition, then exact-min-image prune actual atom pairs. This proves complete coverage without a Cartesian bounding-box stencil. Bounds and rejection decisions are interval-safe; uncertain bins stay included.

The scalar min_image uses approximate inverse wrapping before reconstructing an image. Its floating discrepancy must be included, not hidden by the real-arithmetic proof. Implement `scalar_min_image_search_radius(cert, positions, R)->Result<f64,FallbackReason>` as a conservative bound for this exact existing implementation:

- Let P[k]=max_i abs(r_i[k]), D[k]=2*P[k], Dmax=max D. Bound coordinate subtraction error by eps*D[k]+T.
- For row k, bound the difference between computed t=dot(A_row,fl(delta)) and exact U_row*delta by outward `eta*Dmax + sum_l abs(A[k,l])*subtraction_error[l] + gamma_5*sum_l abs(A[k,l])*(D[l]+subtraction_error[l]) + 5*T`. Call this t_error[k]. Require computed |t| plus its enclosure <2^52 so wrapping rounds to an exactly representable integer; otherwise scalar fallback.
- wrap_half(t)=t-round(t) is exact in this finite domain. Let C0[k]=0.5 and Rinit be an outward upper bound on Cartesian norm of H*C0 evaluated with the scalar expression. Bound the absolute inner enumeration index by `J[k]=ceil(0.5+norm(A_row[k])*Rinit+32*EPS*(1+norm(A_row[k])*Rinit))+1`; enclose every operation and require J[k]<=4096.
- Candidate coefficient subtraction s[k]-n[k] contributes at most `eps*(0.5+J[k])+T`. Bound Cartesian evaluation error by summing abs(H[k,l]) times t_error plus coefficient error and adding the gamma_5 bound for three products/two sums at coefficients <=0.5+J. This gives e_scalar[k]. Set Q=Q_fp+norm(e_scalar), with Q_fp defined above at list R2. If any required domain/operation proof fails, return scalar fallback.
- Reconcile these counts with the actual post-B03 min_image/Vec3 arithmetic before implementing; if source differs, derive a new sufficient upper bound and update the proof before enabling the stencil. Do not alter min_image's physical predicate to simplify the proof.

Certify a fractional interval for every atom via A plus eta. After wrapping and multiplying by m, require its entire interval to occupy one bin (using closed enclosures; a boundary-straddling interval triggers scalar fallback). This intentionally conservative first version avoids duplicate atom insertion and uncertified rounding at bin faces. Interval-safe range reduction handles a periodic seam; ambiguity falls back.

Cost caps: grid bins<=1_000_000, distinct stencil offsets<=4096, `N*stencil_count<=50_000_000`, requested triclinic grid/bin/staging workspace<=134_217_728 bytes, N fits existing i32 cell links. Compute exact candidate pair estimate from occupancies/deduplicated bin targets using checked arithmetic, capped before allocation; if estimate>=0.8*N*(N-1)/2, choose scalar. If every wrapped grid bin is searched, choose scalar. Reason codes are `uncertified_geometry`, `fractional_boundary`, `image_domain`, `grid_cap`, `stencil_cap`, `work_cap`, `memory_cap`, `all_bins`, `candidate_cost`. Bounded builder still stops at its max_pairs; a cap error cannot be converted into success by fallback. Triclinic fallback does not affect orthorhombic behavior.

Internally use `TriclinicBuildStrategy::{Auto,Scalar,FractionalStencil}`. Auto remains Scalar until G4. Add cargo feature `neighbor-diagnostics`; only under it expose `NeighborBuildStrategy` and chainable `with_triclinic_strategy_for_diagnostics(self,strategy)->Self`, plus diagnostics getters. Forced FractionalStencil still applies every safety fallback and cost cap. The benchmark and its scalar control use the same feature/build. No production environment variable is added. Unit/integration strategy tests use that feature explicitly.

## Tolerances, cost campaign and G4

Freeze before evaluation: exact tuple/pair-set equality; reconstructed geometry atol=rtol=2e-12; classical fixture E/F/W atol=rtol=1e-10 in reduced units; MACE atol=rtol=1e-10 and finite difference 1e-6; official cutoff margin 1e-8 Angstrom. Preserve all tighter existing regression tolerances. No silent widening after a failure.

Create dedicated `crates/modena-bench/src/bin/triclinic_neighbors.rs`, `bench/triclinic_neighbors.py`, benchmark tests and documentation. Do not edit live B03 campaign files. After the CPU lease is released, compare scalar and stencil in one release binary, same snapshots/logical order/velocities/compiler/profile, with thread counts 1 and fixed 4. CPU host, OS, compiler, full revision, binary hash, input bytes/hash, precision and environment go into the report.

Predeclared N=(256,2048,8192), number density 0.8, Rforce=2.5 and skin=(0,0.3); fractional stratified positions from grid side ceil(cuberoot(N)) with deterministic rational per-axis jitter in [-0.03,0.03], coordinate offset 0.31, and seed=20261106; serialize snapshots before runs. Use cubic-volume cell families: moderate b_x/L=.2,c_x/L=.1,c_y/L=.15; substantial (.7,.2,.4); unreduced (4.9,0,0); elongated diag factors (2,.5,1). L=(N/.8)^(1/3). Include explicit grid-cap/boundary fallbacks and basis-change rebuild tests separately. Record exact generated positions and require no initial overlap or unresolved cutoff fixture before timings.

Matrix clarification before G4 input generation: basis columns are
`a=(L,0,0), b=(.2L,L,0), c=(.1L,.15L,L)` for moderate;
`a=(L,0,0), b=(.7L,L,0), c=(.2L,.4L,L)` for substantial;
`a=(L,0,0), b=(4.9L,L,0), c=(0,0,L)` for unreduced; and
`a=(2L,0,0), b=(.2L,.5L,0), c=(.1L,.15L,L)` for elongated.
Thus elongated retains the moderate off-diagonal components with the declared
diagonal factors; every determinant is L^3. This resolves the previously
unspecified elongated shear before any timing, without changing density,
velocities, tolerances, sample counts or eligibility criteria.


Measure forced full rebuild and complete VelocityVerlet+neighbor+PairForceModel step/block, not graph construction alone. Include transformation/bin/stencil creation, allocation, sorting, fallback, force reductions and rebuild cadence. For each case create scalar/candidate replicas from identical state and evolve them through 5 warmup paired blocks and 20 measured paired blocks, maintaining each replica continuously between blocks; alternate AB/BA, block size 10 complete steps, dt=.0005 in reduced LJ units, identical deterministic initial velocities with zero center of mass. Serialize initial velocities and all controls. Validate pair/E/F/W before timed blocks and final E/F/W parity; count rebuilds, candidate checks, emitted pairs, fallback causes, requested/retained memory and allocations. Compute paired median ratios and 10_000-resample bootstrap 95% interval with seed=20261106; retain absolute durations and raw paired records.

Promotion requires G1-G3 passed, moderate-tilt N>=2048 intended cases at each enabled thread count and skin regime with median full-step optimized/scalar<=.95 and interval upper<1.0; small/fallback cases must select scalar with <=5% median overhead. Inconclusive noise, unavailable hardware or failed criteria leaves Auto=Scalar and G4 pending/failed. Do not generalize a host-specific selector to unmeasured hardware. For the natural skin=.3 amortized claim, require at least two measured rebuilds per replica; if the fixed campaign does not observe them, label that regime inconclusive rather than infer a rebuild benefit. Per-case wall-time cap is 600 seconds, whole campaign cap is 7200 seconds; a cap abort is incomplete evidence and never a successful performance gate. A candidate default selector is restricted to validated shape/cost conditions on the measured platform; publish exact predicate and re-run it against raw cases before enabling. A fast rebuild with no complete-step gain does not pass G4.

## Required artifacts and official sources

Committed implementation-time docs: `docs/superpowers/specs/2026-11-06-b04-periodic-neighbors-design.md`, `docs/superpowers/plans/2026-11-06-b04-periodic-neighbors.md`, `docs/B04_PERIODIC_GRAPHS.md`, `docs/B04_VERIFICATION.md`. Verification lists G1-G4 separately with commands, revision, snapshot/lock/model hashes, tolerances, outcomes and limits. Large generated timings/venvs remain under ignored target paths with small linked manifests.

Official sources already audited, consulted 2026-11-06: [ASE API](https://docs.ase-lib.org/ase/neighborlist.html), [ASE source](https://docs.ase-lib.org/_modules/ase/neighborlist.html), [MACE constructor v0.3.13](https://raw.githubusercontent.com/ACEsuit/mace/v0.3.13/mace/modules/models.py), [MACE graph](https://raw.githubusercontent.com/ACEsuit/mace/v0.3.13/mace/data/neighborhood.py), [MACE AtomicData](https://raw.githubusercontent.com/ACEsuit/mace/v0.3.13/mace/data/atomic_data.py), [MACE autograd/vector helpers](https://raw.githubusercontent.com/ACEsuit/mace/v0.3.13/mace/modules/utils.py), [MACE dependencies](https://raw.githubusercontent.com/ACEsuit/mace/v0.3.13/setup.cfg), [Torch archived versions](https://pytorch.org/get-started/previous-versions/), [LAMMPS neighbor rationale](https://docs.lammps.org/Developer_par_neigh.html). Mutable [FAIR-Chem calculator](https://raw.githubusercontent.com/facebookresearch/fairchem/main/src/fairchem/core/calculate/ase_calculator.py) includes a capped external graph path; it is not an all-images oracle. No claim that FAIR-Chem, NequIP and MACE graph conventions are interchangeable.

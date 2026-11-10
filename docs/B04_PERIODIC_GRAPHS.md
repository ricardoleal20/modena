# Periodic image graphs and triclinic neighbors (B04)

`PeriodicGraph` builds a directed graph of every supported periodic image
inside a strict Cartesian cutoff. It borrows one immutable `System` snapshot.
The classical `NeighborList` remains an unordered minimum-image half list;
choose the interface whose semantics match the consumer.

## Build a borrowed graph

```rust
use modena_core::{GraphBudget, PeriodicGraph, PeriodicGraphOptions, SimBox, System, Vec3};

fn main() -> Result<(), modena_core::GraphBuildError> {
    let mut system = System::new(SimBox::cubic(8.0));
    system.resize(2);
    system.set_pos(0, Vec3::new(0.1, 0.2, 0.3));
    system.set_pos(1, Vec3::new(1.0, 0.7, 0.9));
    let options = PeriodicGraphOptions {
        cutoff: 3.1,
        pbc: [true, true, false],
        include_periodic_self_images: true,
        budget: GraphBudget {
            max_atom_pairs: 8,
            max_image_candidates: 100_000,
            max_edges: 10_000,
            max_bytes: 1_000_000,
            max_abs_shift: 100,
        },
    };
    {
        let graph = PeriodicGraph::build(&system, options)?;
        for edge in graph.edges() {
            // Integer-image identity is part of the edge, not just its vector.
            println!("{} -> {} {:?}", edge.source, edge.destination, edge.edge_cell_shift);
        }
        println!("{:?}", graph.stats());
    } // end the snapshot borrow before changing the system
    system.set_pos(1, Vec3::new(1.1, 0.7, 0.9));
    let _rebuilt = PeriodicGraph::build(&system, options)?;
    Ok(())
}
```

Budgets above are application choices for this small example, not recommended
limits for arbitrary systems. Positions and cutoff use the same length unit.
The geometry API does not validate chemical arrays or assign energy units.

## Edge and cutoff convention

For cell basis columns `H = [a b c]`, an edge `(i,j,S)` has vector
`r_j - r_i + H*S`, with exact integer `S: [i64; 3]`. The emitted array is
lexicographically sorted by `(source,destination,edge_cell_shift)` and has no
duplicate tuples. Distinct images of one atom pair are separate edges; reverse
directed edges are included when they satisfy the same predicate. There is no
top-K, maximum-neighbor truncation, skin, pair-scale or exclusion lookup.
Positions may be unwrapped when the numerical and shift-domain certificates
succeed; image integers must be retained consistently with that snapshot.

Membership uses the specified f64 Cartesian expression and strict
`(d.x*d.x + d.y*d.y) + d.z*d.z < cutoff*cutoff`. Equality is excluded. Search
bounds conservatively enclose that executable predicate; they do not add a
physical tolerance shell. Keep integer shifts when converting to a model's
graph representation: matching edge counts or Cartesian vectors alone does
not establish tuple identity.

The zero-image self edge `(i,i,[0,0,0])` is always absent. Nonzero periodic
self images are included only when `include_periodic_self_images` is true.
Distinct indices at coincident positions remain distinct graph nodes and can
have zero-length edges. Consumers must apply their own model validity rules.

`pbc[k] = false` requires `S[k] = 0`. All eight axis masks are supported by
this graph API. A full finite, certifiable cell is still required. This mask
is graph metadata; it does not add partial-PBC classical MD, electrostatics,
thermostat or barostat support. B04's actual MACE model gate uses full PBC.

## Resource limits and errors

Construction counts before allocating, then emits in a second pass. Work
counters include both passes: atom-pair visits are `2*N*N`, and candidate
visits count every enumerated image candidate before the final self/cutoff
filter. `max_edges` limits accepted edges, while `max_bytes` charges requested
edge storage, graph header and certified-cell workspace. `requested_bytes`
is a logical storage charge, not allocator usable size or operating-system RSS.

`max_abs_shift` must be in `0..=2^52`. It declares a supported image domain:
the builder proves whole-request coverage within it or rejects the request.
It never clips an otherwise larger image range to return an incomplete graph.
Checked endpoint/product arithmetic and fallible allocation enforce limits.

Errors distinguish invalid input (`InvalidInput`), uncertifiable floating-point
geometry (`NumericalGeometry`), unsupported shifts (`ShiftRange`), exhausted
work/storage budgets (`BudgetExceeded`) and allocation failure (`Allocation`).
Nonfinite coordinates/cells, malformed positions, singular or uncertifiable
cells and unusable cutoff squares cannot produce a graph. Some extremely
large, skewed or ill-conditioned finite inputs intentionally reject. Failed
construction publishes no partial graph and does not mutate the system.

## Rebuild after snapshot changes

The graph stores an immutable borrow of the exact system and options. End the
borrow and rebuild after changing positions, atom count/order/identity, cell,
cutoff, PBC mask or self policy. There is no skin cache or incremental update.
Rust prevents mutation while a graph that borrows the snapshot is in use.

Classical `NeighborList` retains its separate skin, pair-policy and bounded
failure behavior. A zero skin rebuilds on each required update. Triclinic cache
invalidation compares the complete cell basis, including changes with the same
axis-aligned bounding box. The bounded builder remains private; public snapshot
support remains orthorhombic.

## Triclinic strategy and acceptance scope

Auto selects Scalar. The opt-in Cargo feature `neighbor-diagnostics` exposes
`NeighborBuildStrategy` and `with_triclinic_strategy_for_diagnostics` for an
explicit Scalar or FractionalStencil request, plus build diagnostics. Forced
FractionalStencil still applies numerical safety, work, memory and candidate
cost fallbacks. A request is not a guarantee that the stencil executes. There
is no production environment-variable selector.

G1 independent geometry/official ASE, G2 initialized official MACE software
parity and G3 classical correctness have separate approved evidence. The G4
campaign decision and exact scope are in [B04_VERIFICATION.md](B04_VERIFICATION.md).
A successful geometry graph or faster forced rebuild does not establish a
complete-step benefit. The benchmark counts process-wide Rust allocations in
both arms; its atomic overhead is unquantified, so it does not establish
uninstrumented production throughput.

The pinned initialized MACE 0.3.13 instance establishes software parity on the
declared snapshots and arithmetic. It establishes no trained chemical accuracy,
FAIR-Chem checkpoint/task parity, universal NequIP/MACE/FAIR-Chem graph contract
or B10 acceptance. B10 must execute its own model/material/reference gate.

## Recorded contract choices

The verification record preserves failed attempts and the controller's
pre-evaluation decisions, including these costs:

| Choice | Limit retained in the evidence |
| --- | --- |
| Certified f64 coverage with explicit rejection | Some large/skew/ill-conditioned finite inputs reject; interval arithmetic adds implementation cost. |
| Unconditional official tuple guards under Python optimization | The verifier intentionally replaces the original assert-based snippet; count/vector agreement is insufficient. |
| G2 NumPy 1.26.4 compatibility revision | A different ABI/arithmetic baseline was revalidated; no numerical comparison with the unexecuted NumPy 2 candidate is claimed. |
| Scoped `builtins.slice` allowlist for pinned e3nn constants | Trust remains tied to the recorded packaged constant/source, and future import compatibility needs maintenance; unrestricted pickle loading is not enabled. |
| Change-shift mutation increases the first nonzero component's magnitude | The control covers this fixed away-from-zero mutation, not every possible shift mutation; the failed original control remains historical evidence. |
| Private bounded classical builder | Unit and integration coverage remain separate; no public bounded-triclinic snapshot capability is added. |
| Explicit elongated shear and isolated-boundary recipe | Each result covers its declared matrix/control, not arbitrary elongated cells or all boundary geometries. |
| Benchmark-local allocation counters | Overhead is unquantified; retained heap is process-wide across both replicas, with no per-model, per-strategy or RSS attribution. |

The unpublished B04 branch was synchronized onto delivered reviewed B03 with
source/patch identity receipts before later implementation. Historical run
epochs retain their original revisions rather than being relabeled as new
executions. Upstream MACE/TorchScript annotation warnings remain disclosed.

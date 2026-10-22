# Task 2 report: conservative forces, electrostatic pressure and owned caches

Status: complete within Task 2 scope; independent review belongs to the orchestrator.

Implementation commit: `9fa3b9e6fefa326a16a473ccade05709d7e76585` — `:bug: Fix: conservative force models and reusable SPME workspaces`. Post-commit `rtk proxy git status --short` is clean. The local report is stored under the ignored task-evidence directory.

Worktree: `/Users/ricardoleal20/.codex/worktrees/modena-scientific-hardening/modena`; branch `ricardoleal20/modena-scientific-hardening`; immutable task BASE `fdbe886886e0eb143f99276c46debcbea836e6d0`. No subagents, external publication, or modifications to the original checkout/audit evidence. Read the brief, full plan/spec, force handoff, Task 1 integration report, original PHYSICS-NOTES/PERFORMANCE-NOTES and relevant probe references. Applied the test-driven-development and verification-before-completion skills. All shell commands started with `rtk`.

## Delivered behavior

- Tabulated U and radial F come from one cubic Hermite polynomial with slopes `-F` at both endpoints, so the shifted-force pair still satisfies `dU/dr=-F`. Below the first node, the endpoint tangent is extended consistently; the previous flat U/nonzero F defect is gone. Arrays must be finite. `InterpKind::Linear` is retained as a documented legacy selector for the corrected joint interpolation. Node values and shifted cutoff U/F are preserved. Tabulated is now included in the existing cross-potential derivative properties; the midpoint unit expectation was updated to the Hermite endpoint-slope result.
- Direct Ewald uses reciprocal isotropic strain virial `sum_k U_k*(1-k²/(2*kappa²))`. The coordinate-force `sum r_i.F_i` expression was removed from that contribution. Both direct Ewald and SPME add `W_bg=3*U_bg` for the uniform neutralizing background; self energy is independent of strain at fixed kappa. Existing SPME/direct virial tests now require agreement at their finite mesh tolerance, replacing the old test asserting the audited defect.
- `PairScale` supports independent LJ/Coulomb weights. The list removes only pairs with both weights zero; pair providers scale their own energy, force and scalar virial. Binary exclusion constructors remain available. A normalized HashMap index replaces the global linear exclusion scan; ordinary empty-policy lookups return before hashing. Bounds are checked against the live atom count by force providers.
- Electrostatic model/list policies reconcile once per pass without allocating a merged index. Identical overlapping Coulomb weights are applied once; contradictory values are errors. A model-only or neighbor-only Coulomb policy works. Reciprocal correction removes the Gaussian fraction `(1-s)` inside the real cutoff. Outside it, the correction removes the full bare fraction (the omitted erfc tail is accounted for), keeping `U_scaled-U_plain=-(1-s)*q_i*q_j/r` relative to the finite-cutoff model at every tested distance. The analogous force and virial differences use the bare pair derivative. This convention is documented in the correction source. Energy-only direct Ewald passes an empty optional force slice safely and returns the same energy as compute with special pairs.
- Cutoff providers cover PairForceModel, EwaldForceModel, SpmeForceModel, LongRangeForceModel, CompositeForceModel, DomainForceModel and GpuPairForceModel. Composite returns max(pair, real electrostatic cutoff). GPU returns its actual f32 cutoff converted to f64. Removed the Task 1 temporary delegating LJ adapter from ensembles tests. A real two-replica simulation now checks nonzero pair energy and cached forces at `r=3.2` with actual `rc=3.7`, which would be lost by guessing 2.5.
- Triclinic lists use guaranteed-correct all-pair enumeration with Task 1's rigorous minimum image. The bounding-box Cartesian ±1 stencil remains only for orthorhombic cells. Neighbor cache keys compare the complete SimBox, including the basis; an empty System invalidates the cached list before repopulation. Audited skew and nonadjacent-bin fixtures are covered.
- Direct Ewald's k-vector cache compares the full cell basis and k-space cutoff and returns an immutable Arc snapshot. Ensure+use can no longer reacquire a different simulation's vectors after unlocking. Concurrent evaluations in different boxes match separately constructed model references.
- DomainForceModel and CompositeForceModel now honor the trait's accumulation contract for existing nonzero forces and virial. Their drivers must zero once before a complete pass. The manual domain NVE test was corrected to follow this contract explicitly; integrators already zero. Empty pair/SPME models no longer erase existing virial. Domain backend rejects special-pair policies and changes from its construction box; unsupported results are not silently returned.
- SPME has a per-model Mutex-owned workspace holding six FFT plans, both complex grids, reusable per-worker line/FFT scratch and influence factors. The reciprocal lock covers key validation through the complete pass. A new mesh rebuilds owned storage; a changed box or kappa invalidates influence factors even when the rounded mesh is unchanged. The single global influence slot is removed. FFT rows/slabs are grouped into bounded tasks with persistent worker scratch. Standalone reciprocal helpers use temporary local workspaces and have no global cache. Their energy-only route skips inverse/gather work.
- The refactor exposed finite-difference cancellation on the existing tiny-displacement crystalline test. Compensated reciprocal energy/virial summation fixed that numerical sensitivity; the original derivative and strain tolerances remain unchanged. Forces and deterministic serial pair order remain unchanged within existing tolerances. No tensor or constraint impulse stress is claimed.

## Exact public interfaces and compatibility

```rust
pub struct PairScale {
    pub i: usize,
    pub j: usize,
    pub lj: f64,
    pub coulomb: f64,
}
// Re-exported from modena_core root, defined in neighbor.
// Derives Debug, Clone, Copy, PartialEq.

impl NeighborList {
    pub fn with_pair_scales(cutoff: f64, skin: f64, scales: &[PairScale]) -> crate::error::Result<Self>;
    pub fn pair_scale(&self, i: usize, j: usize) -> (f64, f64);
    pub fn pair_scales(&self) -> &[PairScale];
    pub fn validate_pair_scales(&self, n: usize) -> crate::error::Result<()>;
    // Existing binary 0/0 compatibility APIs retain Self signatures:
    pub fn new(cutoff: f64) -> Self;
    pub fn with_skin(cutoff: f64, skin: f64) -> Self;
    pub fn with_skin_exclusions(cutoff: f64, skin: f64, exclusions: &[(u32,u32)]) -> Self;
}
impl EwaldForceModel {
    pub fn with_pair_scales(self, scales: Vec<PairScale>) -> crate::error::Result<Self>;
    pub fn with_excluded_pairs(self, pairs: Vec<(usize,usize)>) -> Self;
    // new(kappa,real_cutoff,kspace_cutoff)->Result<Self> unchanged.
    // energy(&System,&mut NeighborList)->Result<f64> unchanged.
}
impl SpmeForceModel {
    pub fn with_pair_scales(self, scales: Vec<PairScale>) -> crate::error::Result<Self>;
    pub fn with_excluded_pairs(self, pairs: Vec<(usize,usize)>) -> Self;
    // new(kappa,real_cutoff,kspace_cutoff)->Result<Self> unchanged.
}
// ForceModel::neighbor_cutoff(&self)->Option<f64> seam unchanged.
```

Constructor policy is explicit: pair order normalizes to `(min,max)`; distinct indices and finite nonnegative weights are required; amplification >1 is allowed (there is no undocumented upper bound). Identical normalized duplicates coalesce; contradictory duplicates return Config errors. Missing pairs return `(1,1)`. New neighbor constructors require finite positive cutoff, finite nonnegative skin and finite sum. Electrostatic constructors now reject infinite/NaN/nonpositive parameters. The new Result APIs are intended for production input paths. Legacy Self-returning constructors preserve their signature and panic on invalid arguments; bounds require a live atom count and are checked during force evaluation or explicitly through `validate_pair_scales`.

The neighbor list is authoritative for LJ scaling. The electrostatic builder consumes Coulomb weights; its LJ member is validated but does not reweight the separate pair provider. Passing the same complete policy to the neighbor constructor and electrostatic builder is supported and does not double-scale. The former mutable public `excluded_pairs` vectors on electrostatic models are replaced with private validated indexed policy; existing builder callsites remain compatible. No Task 1 production files were edited. Only its ensembles test adapter/integration fixture changed.

Task 3/4 should use the new Result constructors, obtain the max cutoff from the actual model, and propagate input errors. Task 3 must retain the constraint-pressure limitation and guard constrained/anisotropic NPT. Full runner/config checkpoint and graph/topology policy wiring remain the later tasks' ownership.

## Observed RED evidence

Before the corresponding production physics fixes:

```text
rtk proxy cargo test -p modena-core --test force_hardening -- --nocapture
exit 101: 0 passed, 9 failed
Tabulated r=0.7: force=24, -dU/dr=26.9999999957804
PairForceModel cutoff: None instead of Some(3.7)
Ewald energy-only special pair: panic indexing empty force slice
Background W=0, independent expected -0.017069481689064907
Domain existing virial erased
Skew nonadjacent-bin list: [] instead of [(0,1)]
Ewald same-bounding-box/different-basis cache mismatched fresh model
Neutral direct W=-0.3444518957140921, strain -0.3977306489988397
Nonneutral SPME W=-0.42625083625577076, strain -0.43119997664153864
```

Weighted-policy tests were first compiled against missing PairScale/with_pair_scales APIs (exit 101: unresolved import/methods). After introducing validated constructors and list policy but before force scaling/corrections, observed behavioral RED:

```text
rtk proxy cargo test -p modena-core --test pair_scales -- --nocapture
exit 101: 1 constructor-policy test passed, 4 failed
Independent LJ energy failed to scale
Conflicting model/list Coulomb policies wrongly accepted
Ewald r=1.3 c=0: U=-0.5505602811895959, expected -0.012098742728057443
SPME weighted bare-fraction correction failed
```

Owned-cache regression before workspace changes:

```text
rtk proxy cargo test -p modena-core --test force_cache -- --nocapture
exit 101: 0 passed, 1 failed
SPME warm compute allocations=66; required 0
```

Additional behavioral RED before defensive fixes:

```text
rtk proxy cargo test -p modena-core --test force_hardening -- --nocapture
exit 101: 9 passed, 2 failed
Electrostatic constructors accepted infinity
Empty then repopulated System reused [] instead of [(0,1)]
```

Covering tests initially caught a numerical regression during the FFT refactor:

```text
rtk proxy cargo test -p modena-core --lib --test force_cache --test force_hardening --test pair_scales --test property_consistency --test electrostatics --test decomp --test ensembles -- --nocapture
exit 101: library 79 passed, 1 failed (integration harnesses not reached)
Tiny displaced crystal atom0/axis0: mesh force 0.000022890638422921417,
finite-difference force 0.00002297895207448164.
```

Compensated summation repaired this without loosening the original test. Subsequent six SPME reference tests and final library suite passed.

## Independent expectations and numerical gates

- Cubic table: raw U=r³, F=-3r², with nodes at r=[1,2,3] and a nonzero cutoff force. Hermite must reproduce the cubic exactly on the table, and central finite differences with h=1e-6 must agree with force within 1e-7, including lower tangent extrapolation. Existing arbitrary table nodes also pass gradient/finiteness/cutoff properties.
- Isotropic strain uses box AND positions scaled by lambda, fixed kappa, reciprocal cutoff/mesh fixed away from truncation crossings. Scalar W must be `-dU/dlambda` at lambda=1. Neutral direct and charged full SPME tests use h=1e-6 and absolute 1e-7. Single charged direct system with no permitted nonzero k-vectors isolates `W_bg=3U_bg` to absolute 1e-14.
- Special pairs at r=1.3 and r=5 (>real rc4.2), weights c=[0,0.4,1], compare independently known bare `q_i*q_j/r`, its pair derivative and strain virial to total energy/force/virial differences, with absolute 1e-12. Both model+list and neighbor-only policies are checked; SPME additionally uses finite strain at h=1e-6, tolerance 1e-7.
- Skew fixture a=(10,0,0), b=(49,10,0), c=(0,0,40), displacement=(19.7,0.1,0), has independently known image=(-0.3,0.1,0), while Cartesian bins are nonadjacent. Audited b=(4.9,10,0), displacement=(6.37,3,0) gives image=(-3.63,3,0), r=4.709235606762524 <4.8. Both lists include the pair.
- Shared-model evaluations in boxes of edges 9 and 11 run concurrently 12 times each for direct and SPME, compared exactly against separate-model references. Cache same-mesh edge change 9->9.01 and kappa 0.8->0.9 compare exact full U/W/force tuples against fresh models.

## GREEN verification and allocation evidence

```text
rtk proxy cargo test -p modena-core --lib --test force_cache --test force_hardening --test pair_scales --test property_consistency --test electrostatics --test decomp --test ensembles
exit 0: 125 passed, 0 failed.
Counts: lib80; decomp5; electrostatics9; ensembles4;
force_cache1; force_hardening11; pair_scales6; property_consistency9.

rtk proxy cargo test -p modena-core --test water
exit 0: 4 passed, including excluded-pair free flight and full rigid-water/SPME NVE. Harness 66.28 s.
Combined final default coverage: 129 passed (125 above + water4).

rtk proxy cargo test -p modena-core --test force_cache -- --nocapture
exit 0: 1 passed
SPME warm compute allocations=0
SPME warm parallel compute allocations=0

rtk proxy cargo test -p modena-core --lib spme::tests -- --nocapture
exit 0: 6 passed (original derivative/strain tolerances unchanged).

rtk proxy cargo test -p modena-core --features gpu --test gpu_pair gpu_rejects_unsupported_lj_scales_before_dispatch -- --nocapture
exit 0: 1 passed, 6 unrelated filtered. Actual local adapter obtained; no skip message.

rtk proxy cargo clippy --workspace --all-targets -- -D warnings
exit 0.
rtk proxy cargo clippy --workspace --all-targets --features gpu -- -D warnings
exit 0.
rtk proxy cargo clippy --workspace --all-targets --features mpi -- -D warnings
exit 0. Homebrew mpicc at /opt/homebrew/bin/mpicc.
rtk proxy cargo fmt --all -- --check
exit 0.
rtk proxy git diff --check
exit 0.
```

A dedicated integration binary installs a counting global allocator and contains one test, avoiding concurrent test-harness allocation contamination. After one warmup compute, the second stable-box/mesh pass counts every allocation/reallocation. Serial case uses one Rayon thread, edge9/kcut6.4 ->20³ mesh; parallel case uses four Rayon threads, edge20/kcut6.4 ->45³ mesh (above the parallel threshold). Both count zero in the complete compute including the reusable skin neighbor list. Before caching, the serial case counted 66. This is allocation/reuse evidence, not a wall-time speedup claim. Different mesh/box/kappa invalidation legitimately allocates/rebuilds owned cache data. Parallel pair-force scratch was not speculatively rewritten: no new profile/matched throughput evidence justified extending its existing thresholded path in this task. Empty-scale fast path avoids hashing and builds no index buckets.

## Restrictions and remaining scope

Pressure is still scalar isotropic pressure; constraint impulse stress and a full stress tensor are absent. Task 3 guards constrained and anisotropic coupling and reports the diagnostic pressure limitation. SPME remains orthorhombic; direct Ewald and pair models use rigorous triclinic images, with intentionally O(N²) triclinic list fallback. Ill-conditioned lattice searches can be expensive. Ordinary finite erfc/kspace truncation and finite-mesh errors remain numerical controls; bare-fraction special-pair differences do not certify every chemical force-field convention or arbitrary low-accuracy mesh.

DomainForceModel rejects all special-pair/exclusion policies and box changes. GPU rejects fractional/independently excluded LJ weights; its CellList gather accepts binary 0/0 omissions, while AllPairs rejects any LJ weight differing from 1. Coulomb weights are not consumed by a GPU LJ provider. GPU is still f32 host-offload LJ; no new electrostatics/resident-trajectory claim. MPI transport API is unchanged (it does not accept pair policies); compile/lint verified, no multi-rank execution claimed in Task 2. Production acceleration integration and full branch benchmarks remain later work.

SPME models serialize shared-instance reciprocal evaluation via their workspace mutex. Separate model instances no longer thrash a global influence cache. An intentionally shared model with alternating mesh sizes may rebuild its single owned workspace; correctness is covered but no shared-instance parallel speedup is promised. Ewald snapshots allow safe concurrent reciprocal evaluation. Exact split continuation belongs to Task 3; tabulated interpolation and weighted policies change the Hamiltonian relative to the audited baseline, while corrected strain/background terms change its reported pressure.

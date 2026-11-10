# B07 Tensor Pressure Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use `subagent-driven-development` after controller handoff, or `executing-plans` if the controller executes inline. User-authorized execution is sequential: one owner for each task followed by independent specification and code review before advancing. No parallel repository writers. The controller adopts spec/plan and performs source preflight directly; there is no separate scaffolding-agent task.

**Goal:** deliver independently validated six-component force virials and constraint-inclusive endpoint/interval pressure with exact stage accounting, continuation, and explicit missing-capability semantics.

**Architecture:** retain the scalar force interface and add an additive full evaluation seam; share CPU physics kernels across scalar/tensor paths. Record actual per-edge projection impulses and solve endpoint constraint accelerations independently. Keep pressure as a versioned sibling of legacy scalar Energies; validate analytic dynamics and sampling against independent energy and external references.

**Tech Stack:** existing Rust std/serde/serde_json/libm/rayon/rustfft, Python standard library for campaign orchestration, pinned CPU double LAMMPS with RIGID+KSPACE; no new Rust dependency.

**Spec:** `docs/superpowers/specs/2026-11-06-b07-tensor-pressure-design.md` (all scientific equations, exact DTOs, API migrations, tolerance tables, and budgets are binding).

## Global Constraints

- Adopted 2026-11-10 on reviewed integrated main `d6cd33ca3418b1617ab71fa2cf7416a24d27f405`; preparation itself performed no implementation or measurements. Current user explicitly authorizes B07 implementation, PR and merge. Historical preparation source is retained below.
- Source epoch inspected: B04 `7ea653641c813bc2f01476805db2990e0de4677a`; provisional until final reviewed B04/B05 synchronization. Future branch `ricardo/b07-tensor-pressure`.
- Preserve scalar `Virial {w}`, `ForceModel::compute`, legacy `Energies`, conservative force caches, thermostat RNG sequence, B02 scratch and B03 packed paths.
- No new Rust dependency; production SPME remains orthorhombic and analytic order4. No inferred isotropic tensor or W/3 anisotropic controller path.
- Exact B05 twelve profile fields and existing compute_profiled remain; final B05 synchronization precedes all implementation.
- All shell commands begin `rtk`; all numerical campaigns require release of other streams' exclusive measurement lease.
- All supported-case analytic, six-strain, impulse, timestep, statistical, and external gates are mandatory. Pending/reference failures cannot count as acceptance.
- Freeze protocol/snapshots/tolerances/work limits before outcomes; preserve all failures. No adaptive tolerance relaxation, seed replacement, or unbounded campaign.
- The spec's explicit Result-returning Barostat::apply migration is intentional; no panic/no-op substitute. Production anisotropic/constrained controller rejection remains.
- All tasks include one fresh failing-test cycle, focused implementation, the indicated verification, and a Conventional Commit per current COMMIT_RULES.md. Do not re-run broad tests after unchanged code without a new concern.

---

## Controller adoption (not an implementation task)

- [ ] Wait for final reviewed B04 and B05 handoff before touching code. Read their final plan/protocol, `AGENTS.md`, RTK, relevant git status/diffs; record source hashes and ownership in spec §1. Do not move B04/B05 HEAD or edit their active checkout.
- [ ] Create/reuse the B07 isolated worktree on the integrated reviewed base, branch `ricardo/b07-tensor-pressure`. Copy these two drafts to their future repository paths. Record source epoch, compiler/profile/features, and exact B05 twelve field names. Adoption does not claim scientific validation.
- [ ] Review the explicit migrations (Barostat::apply Result, RunSpec optional field, RunEvent variant), Nyquist mathematical extension, cost caps, and test/reference domain. User already authorized sequential owner+review; begin Task1 after adoption without another optional permission gate.
- [ ] Use read-only inventory commands:

```sh
rtk git status --short
rtk git rev-parse HEAD
rtk rg -n 'compute_profiled|struct SpmeProfile|impl Barostat|\.anisotropic\(|Energies \{|RunSpec \{' crates python
rtk rg -n 'current_virial|RunEvent::|ContinuationState \{' crates/modena-core
```

## File ownership and order

| Area/files | Responsibility / first task |
|---|---|
| `src/stress.rs`, `src/lib.rs`, `src/force_model.rs` | tensor and additive evaluation, Task1 |
| `src/bonded.rs`, `src/composite.rs` | connected-image bonded tensors and composition, Task2 |
| `src/ewald.rs`, `src/external.rs` | electrostatics/cell derivative and provider six stress, Task3 |
| `src/spme.rs` | mesh tensor plus B05 compatibility, Task4 |
| `src/constraint_stress.rs`, `src/constraints.rs`, `src/vec3.rs` | periodic edge geometry and actual impulse report, Task5; endpoint solve Task6 |
| `src/integrator.rs` | shared step engine + VV Task7; BAOAB Task8 |
| `src/runner.rs`, `src/thermo.rs`, `src/io/input.rs`, `src/config/validation.rs`, `src/barostat.rs` | typed output/cache/controller compatibility, Task9 |
| `src/io/restart.rs` | versioned sibling stress checkpoint, Task10 |
| `tests/stress_{tensor,bonded,electrostatics,mesh,constraints,endpoint,vv,baoab,runner,restart,statistics}.rs` | production-facing tests in their respective tasks |
| `tests/support/stress.rs` | small reproducible fixture helpers; first Task1 and extended only by owning task |
| `examples/stress_validation.rs`, `examples/stress_validation/{mod,fixtures,mesh_oracle,statistics,protocol}.rs` | independent oracles, actual deterministic/statistical campaign; introduced with first need Tasks4/11 |
| `bench/b07_pressure_reference.py`, `bench/fixtures/b07/`, `docs/scientific-hardening/B07-PRESSURE-VALIDATION.md` | external executable adapter, frozen snapshots/manifest/results, Task12 |

Rust integration-test snippets import their named modena_core modules/types and the explicitly defined support helpers with `#[path="support/stress.rs"] mod support; use support::*;`; ordinary `use` declarations are not new APIs. Unit-test snippets define the same local `near(a:f64,b:f64,t:f64)` assertion when test-support imports are not available.

All `src`, `tests`, `examples` paths above are under `crates/modena-core/`. No edits to B05's independent Ewald oracle; consume its final mathematical protocol as context, never call it as the mesh energy oracle. No generic repository restructuring.

### Task 1: Symmetric tensor and additive pair evaluation

**Files:** create `crates/modena-core/src/stress.rs`, `crates/modena-core/tests/stress_tensor.rs`, `crates/modena-core/tests/support/stress.rs`; modify `src/lib.rs`, `src/force_model.rs`.
**Consumes:** existing System/Vec3/NeighborList/PairPotential and exact spec §6.
**Produces:** SymTensor6, MissingTensorReason/MissingTensor/TensorAvailability, ModelEvaluation, `ForceModel::evaluate(&self,&mut System,&mut NeighborList)->Result<ModelEvaluation>`, pair shared pass. Later record types are added only in their owning task.

- [ ] Write the exact numeric failing test; helper `complete` must panic on Missing rather than treating it as zero. These helpers go in `tests/support/stress.rs` and are imported by path in later tests.

```rust
use modena_core::{system::{System,SimBox},vec3::Vec3};
use modena_core::potential::PairPotential;
use modena_core::stress::{SymTensor6,TensorAvailability};
pub struct Harmonic;
impl PairPotential for Harmonic {
    fn cutoff(&self)->f64 { 8.0 }
    fn cutoff2(&self)->f64 { 64.0 }
    fn energy_force(&self,r2:f64,_:u32,_:u32)->(f64,f64) { (0.5*r2,-1.0) }
}
pub fn atoms(p:&[[f64;3]],m:&[f64])->System {
    let mut s=System::new(SimBox::cubic(20.0)); s.resize(p.len());
    s.positions=p.iter().flatten().copied().collect(); s.masses=m.to_vec(); s
}
pub fn complete(x:&TensorAvailability)->SymTensor6 {
    match x { TensorAvailability::Complete(t)=>*t, _=>panic!("tensor unavailable") }
}
pub fn near(a:f64,b:f64,t:f64) { assert!((a-b).abs()<=t*(1.0+b.abs()),"{a} != {b}"); }
#[test]
fn pair_tensor_is_additive_and_not_isotropic() {
    use modena_core::{force_model::{ForceModel,PairForceModel},neighbor::NeighborList};
    let mut s=atoms(&[[0.,0.,0.],[1.,2.,3.]],&[1.,3.]);
    s.forces.fill(7.); let model=PairForceModel::new(Harmonic);
    let mut nl=NeighborList::new(8.);
    let a=ForceModel::evaluate(&model,&mut s,&mut nl).unwrap();
    assert_eq!(a.energy,7.); assert_eq!(a.scalar_virial,-14.);
    assert_eq!(complete(&a.tensor).0,[-1.,-4.,-9.,-2.,-3.,-6.]);
    assert_eq!(s.forces,vec![8.,9.,10.,6.,5.,4.]);
    let b=ForceModel::evaluate(&model,&mut s,&mut nl).unwrap();
    assert_eq!(b,a); assert_eq!(s.forces,vec![9.,11.,13.,5.,3.,1.]);
}
```

- [ ] Run `rtk cargo test -p modena-core --test stress_tensor`; expect compile failure for the new type/method. Add a legacy ForceModel implementation with only compute returning E2/W3 and verify evaluate is Missing while preloaded forces remain additive; known empty pair returns Complete(ZERO).
- [ ] Implement SymTensor6 directly: diagonal a_i*b_i; off-diagonal .5*(a_i*b_j+a_j*b_i). Add default trait method exactly as spec. Refactor pair serial/parallel/deterministic passes to one tensor-aware accumulator; add each dyad alongside the existing scalar summation, preserving half-list and scratch leases. Do not compute tensor by reading total system.forces.

```rust
let contribution=SymTensor6::dyad(d,f_i);
partition.tensor=partition.tensor.plus(contribution);
// ModelEvaluation owns one pass; caller W/forces remain independently additive.
```

- [ ] Run focused tests plus `rtk cargo test -p modena-core --test force_hardening` and final-base B02/B03 pair scratch/packed tests. Add six-strain pair FD helper using independent cell+coordinate affine deformation and the spec eps ladder; assert all components, trace, 1/4 threads, cold/reused scratch, finite values.
- [ ] Commit explicit files: `rtk git add crates/modena-core/src/stress.rs crates/modena-core/src/lib.rs crates/modena-core/src/force_model.rs crates/modena-core/tests/stress_tensor.rs crates/modena-core/tests/support/stress.rs`; `rtk git commit -m 'feat(core): additive tensor force evaluation'`. Reviewer confirms shared arithmetic/no 3N scratch increase.

### Task 2: Bonded tensors and completeness-preserving composition

**Files:** modify `src/bonded.rs`, `src/composite.rs`, `src/stress.rs`; create `tests/stress_bonded.rs`.
**Consumes:** Task1 SymTensor6/ModelEvaluation and existing Bond/Angle/Dihedral types.
**Produces:** built-in bonded/long-range/composite evaluate overrides and `ModelEvaluation::combine(self,other)->Self`.

- [ ] Write production-facing angle test with exact expected nonzero shear and zero scalar virial; an isotropic reconstruction must fail.

```rust
#[test]
fn right_angle_has_zero_trace_and_nonzero_shear() {
    use modena_core::{bonded::{Angle,BondedTopology,BondedForceModel},force_model::ForceModel,neighbor::NeighborList};
    let mut s=atoms(&[[1.,0.,0.],[0.,0.,0.],[0.,1.,0.]],&[1.,1.,1.]);
    let mut top=BondedTopology::new();
    top.angles.push(Angle{i:0,j:1,k:2,theta0:std::f64::consts::PI/3.,ktheta:3.});
    let m=BondedForceModel::new(top);
    let e=ForceModel::evaluate(&m,&mut s,&mut NeighborList::new(8.)).unwrap();
    let expected=[0.,0.,0.,std::f64::consts::PI/2.,0.,0.];
    for (a,b) in complete(&e.tensor).0.into_iter().zip(expected) { near(a,b,2e-12); }
    near(e.scalar_virial,0.,2e-12);
}
```

- [ ] Run `rtk cargo test -p modena-core --test stress_bonded`; expect Missing or shear mismatch before override. Add combine test: complete[-1,-4,-9,-2,-3,-6] plus Missing(legacy) keeps energy/scalar sums but Missing tensor; empty known zero preserves all six values.
- [ ] Refactor term returns to `(energy,SymTensor6)` or a private named contribution record; preserve the original local connected vectors. Sum `(q_a-q_anchor)⊙f_a` only from that term's forces. Do not use accumulated global forces or independently wrap each dihedral atom. Add full harmonic bond, rotated angle, generic proper/improper dihedral, integer-image/translation tests and six-strain FD away from branch changes.

```rust
let tensor=local_positions.iter().zip(term_forces.iter())
    .fold(SymTensor6::ZERO,|acc,(q,f)|acc.plus(SymTensor6::dyad(*q-anchor,*f)));
```

- [ ] Override current concrete CompositeForceModel and LongRangeForceModel without changing their struct shape; combine each active evaluation exactly once. Run `rtk cargo test -p modena-core --test stress_bonded` and existing bonded/composite tests located in the final base; expected all tensor/trace/additive gates pass.
- [ ] Commit owned files with `:sparkles: Add: bonded and composite tensor virials`; independent review checks image anchors, term-local forces, empty/missing semantics.

### Task 3: Direct Ewald, correction/background, and provider stress

**Files:** modify `src/ewald.rs`, `src/external.rs`; create `tests/stress_electrostatics.rs`; extend existing external provider tests without starting real workers.
**Consumes:** Task1 evaluation seam; exact existing Ewald real/correction policy and ASE protocol.
**Produces:** six direct reciprocal, real, self/background and scaled-exclusion tensors; trait provider evaluate retaining all six wire stress values. Preserve public provider inherent evaluate.

- [ ] Add internal `fn reciprocal_mode_tensor(u:f64,k:Vec3,kappa:f64)->SymTensor6` used by production and this numeric test; the independent integration test remains the six-strain energy derivative (the helper test alone cannot certify physics).

```rust
#[test]
fn reciprocal_mode_has_exact_trace_and_oblique_shear() {
    let s=reciprocal_mode_tensor(2.0,Vec3::new(1.,2.,3.),0.5);
    for (a,b) in s.0.into_iter().zip([-16./7.,-106./7.,-256./7.,-60./7.,-90./7.,-180./7.]) {
        near(a,b,2e-14);
    }
    near(s.trace(),-54.,2e-14);
}
```

- [ ] Run `rtk cargo test -p modena-core reciprocal_mode_has_exact_trace_and_oblique_shear`; expect missing helper. Add charged background test with U_bg=-π*.5²/(2*720*.8²): each diagonal U_bg, each shear zero, trace3U_bg. External conversion mock σ=[[1,.4,.5],[.4,2,.6],[.5,.6,3]], V_A³=8, energy_scale=2 must yield [-16,-32,-48,-6.4,-8,-9.6], W=-96.
- [ ] Implement real/correction dyads from the actual force returned in each smooth branch; make `exclusion_correction` return tensor alongside its existing energy/trace through one shared kernel. Implement reciprocal U_k*[I-Akk], retain existing mode multiplicities/cutoff membership, self zero, background U_bg I. Keep scalar numeric contracts.

```rust
let a=2.0/k.norm2()+1.0/(2.0*kappa*kappa);
let tensor=SymTensor6([u,u,u,0.,0.,0.])
    .plus(SymTensor6::dyad(k,k).scaled(-u*a));
```

- [ ] Convert provider stress using all six entries before mutation; preserve symmetry1e-12/finiteness/protocol identity and unit checks, overflow/session rollback. Tests include missing wire stress error, nonfinite/shear-asymmetric rejection, repeated additive calls and no second worker request for tensor.
- [ ] Run `rtk cargo test -p modena-core --test stress_electrostatics`, `rtk cargo test -p modena-core --test pair_scales`, existing external session/force-cache tests. Six-strain fixtures include neutral/net-charge, scales/exclusions inside/outside rc, cold/warm cache and triclinic direct Ewald with fixed integer modes. Expected gates spec §9.
- [ ] Commit `:sparkles: Add: reciprocal and external stress tensors`; review checks no double scaling, provider sign, fixed-κ derivative, reciprocal multiplicity.

### Task 4: Independent arbitrary-cell mesh energy and production SPME tensor

**Files:** modify `src/spme.rs`; create `examples/stress_validation/mod.rs`, `examples/stress_validation/mesh_oracle.rs`, `tests/stress_mesh.rs`; retain B05 profile tests.
**Consumes:** Tasks1/3, final B05 compile-time profile observer and twelve fields.
**Produces:** public `spme_reciprocal_tensor(kappa:f64,box_size:&SimBox,positions:&[f64],charges:&[f64],forces:&mut[f64],mesh:[usize;3])->Result<(f64,SymTensor6)>`; `SpmeForceModel::evaluate_profiled(&self,&mut System,&mut NeighborList)->Result<(ModelEvaluation,SpmeProfile)>`; independent test/example-only functions below.

```rust
pub fn mesh_energy(charges:&[f64], fractional:&[[f64;3]], cell:[[f64;3];3],
    mesh:[usize;3],kappa:f64)->Result<f64,String>;
pub fn grid_energy(grid:&[f64], cell:[[f64;3];3], mesh:[usize;3], kappa:f64)
    ->Result<f64,String>;
pub fn strain_cell(cell:[[f64;3];3],component:usize,eta:f64)->[[f64;3];3];
pub fn energy_strains<F:Fn([[f64;3];3])->Result<f64,String>>(cell:[[f64;3];3],eta:f64,f:F)
    ->Result<[f64;6],String>;
```

Cell arrays are columns. mesh_energy spreads charges via explicit order4 polynomials w0=(1-t)^3/6,w1=(3t³-6t²+4)/6,w2=(-3t³+3t²+3t+1)/6,w3=t³/6 to indices floor(Ns)-1..+2. grid_energy loops every grid point/mode using sin/cos DFT, independently inverts H, enumerates equal ± Nyquist lifts; no production math/FFT helper calls. strain_cell applies spec's half-shear matrix. All errors are finite/shape/singular-cell/work-limit strings; deterministic tests unwrap after validation.

- [ ] Write odd/even oracle tests before production tensor. The exact Nyquist synthetic grid test must differentiate actual energy and assert nonzero surviving yz component in the second fixture:

```rust
#[test]
fn nyquist_extension_cancels_only_affected_shear() {
    let mesh=[4,6,8]; let cell=[[8.,0.,0.],[0.,9.,0.],[0.,0.,10.]];
    let grid:Vec<f64>=(0..4).flat_map(|x|(0..6).flat_map(move|y|(0..8).map(move|z|
        (std::f64::consts::PI*x as f64+2.*std::f64::consts::PI*(y as f64/6.+z as f64/8.)).cos()))).collect();
    let d=energy_strains(cell,1e-5,|h|grid_energy(&grid,h,mesh,.8)).unwrap();
    near(d[3],0.,2e-6); near(d[4],0.,2e-6); assert!(d[5].abs()>1e-4);
}
```

- [ ] Run `rtk cargo test -p modena-core --test stress_mesh`; expect missing oracle/public reciprocal tensor. Add second mixed mode(Nx/2,Ny/2,1) with all three offdiag derivatives zero. Test conjugate symmetry under general shear by comparing the direct complex potential's imaginary norm<=1e-12*scale, not only real-energy summation.
- [ ] Implement independent oracle then prove H0 energy/coordinate force/trace against production on grids[5,5,5],[5,9,5],[4,6,8],[6,6,6] and spec charged positions. Add all six strains four eps with fixed integer grid and κ. Only after H0 passes implement orthorhombic cache six factors, zeroing ab only when a or b is Nyquist. Use same stored fac and |Q|² energy, no new force differentiator.

```rust
let shear=if nyquist[a] || nyquist[b] {0.0} else {-u*factor*k[a]*k[b]};
```

- [ ] Thread tensor through B05 ordinary/profile common implementation; compute_profiled returns old tuple and evaluate_profiled retains tensor. Preserve all12 fields; include tensor time in existing phases. Add ordinary/profile E/F/W/S parity, additive prior forces, cold/reuse/cell-change influence invalidation. Run `rtk cargo test -p modena-core --test stress_mesh`, `rtk cargo test -p modena-core --test spme_profile`, `rtk cargo test -p modena-core --test force_cache`, `rtk cargo test -p modena-core --lib spme::tests`.
- [ ] Run the fixed mesh-translation/Ewald convergence gate spec §9, not a benchmark sweep. If no independent energy derivative meets tolerance, stop this task as failed; no arbitrary convention switch after outcomes.
- [ ] Commit `:sparkles: Add: independently differentiated mesh stress`; review proof of conjugacy/orthorhombic equivalence plus actual odd/even derivative data. Record cache memory+six accumulation cost.

### Task 5: Periodic constraint geometry and exact projection impulses

**Files:** create `src/constraint_stress.rs`, `tests/stress_constraints.rs`; modify `src/constraints.rs`, `src/lib.rs`, `src/vec3.rs`; extend support fixtures.
**Consumes:** Task1 dyads, B04 periodic geometry; spec §3 and exact geometry/report methods §6.
**Produces:** ConstraintGeometry/ConstraintEdge/ConstraintImpulse/ProjectionReport/ProjectionStage, constraint_geometry, Constraints::shake_report/rattle_report; wrappers retain signatures and skip ledger storage when unused.

- [ ] Add rotor fixture helper `rotor(m1:f64,m2:f64)->(System,Constraints)` to support: use exact e,t,d,omega,R,V in spec §9, place masses about COM, velocities V plus mass-weighted relative u, one bond. Add failing per-atom impulse reconstruction test:

```rust
#[test]
fn rattle_report_reconstructs_actual_momentum_change() {
    let (mut s,c)=rotor(1.,3.); let g=constraint_geometry(&s,&c).unwrap();
    let r=g.edges[0].reference; s.add_vel(1,r*0.2);
    let before=s.velocities.clone();
    let rep=c.rattle_report(&mut s,&g,ProjectionStage::VvFinal,1e-13,100).unwrap();
    let mut dp=vec![Vec3::zero();s.n()];
    for j in &rep.impulses { let e=&g.edges[j.edge];dp[e.j]+=j.impulse_on_j;dp[e.i]-=j.impulse_on_j; }
    for i in 0..s.n() { for a in 0..3 {
        near(dp[i].0[a],s.masses[i]*(s.velocities[3*i+a]-before[3*i+a]),5e-10);
    }}
    near(rep.tensor_impulse.trace(),-0.2535,5e-10); // -mu*d²*.2
}
```

- [ ] Run `rtk cargo test -p modena-core --test stress_constraints`; expect missing API. Add position projection test: unwrapped predicted q=q0+τv, call shake_report; reconstructed J_i=m_i(q_after-q_pred)/τ, not /τ². Check every iteration count/list entry from actual updates; no dummy count assertion.
- [ ] Build graph images before prediction, validate ties/cycle closure/nonzero winding/duplicates/positive masses, carry frozen shifts and current continuous edge vectors. Relative stage motion criterion spec §6 excludes image-losing drift while allowing whole-molecule translation. Implement rollback on report failure; test overflow/nonconvergence leaves source arrays unchanged and returns no report.

```rust
// Loose SHAKE actual update: atom i += g*r_ref/m_i; atom j -= g*r_ref/m_j.
let impulse_on_j=r_ref*(-g/interval);
report.impulses.push(ConstraintImpulse{edge,lever:r_ref,impulse_on_j});
report.tensor_impulse=report.tensor_impulse.plus(SymTensor6::dyad(r_ref,impulse_on_j));
```

- [ ] Instrument every loose/water Newton/velocity increment with its actual signs; water lam uses +lam*r_ref/τ. Test water Newton versus3 explicit constraints, unequal H masses, face/corner straddles, translation, full unsymmetrized torque and momentum; analytic force comes Task6, not falsely claimed by path agreement.
- [ ] Run `rtk cargo test -p modena-core --test stress_constraints` and existing constraint/integrator tests. Require geometry and impulse tolerances independently; old settle initialization returns no physical impulse. Commit `:sparkles: Add: periodic constraint impulse ledgers`; review actual update/ledger correspondence and rollback.

### Task 6: Independent endpoint constraint reaction and mechanical pressure

**Files:** modify `src/constraint_stress.rs`, `src/stress.rs`; create `tests/stress_endpoint.rs`; extend support water fixtures.
**Consumes:** Task5 geometry; final conservative system.forces; final tangent velocities; ModelEvaluation.
**Produces:** EndpointReaction, PressureObservation, ValidationStatus, kinetic_tensor, endpoint_reaction and endpoint_pressure exact signatures spec §6. No conservative cache mutation.

- [ ] Add test checking analytic unequal-mass rotor tensor and nonzero COM pressure, not merely distance residual:

```rust
#[test]
fn endpoint_rotor_cancels_internal_scalar_pressure() {
    let (s,c)=rotor(1.,3.); let g=constraint_geometry(&s,&c).unwrap();
    let f=s.forces.clone(); let reaction=endpoint_reaction(&s,&g).unwrap();
    near(reaction.multipliers[0],.3675,2e-9); // mu*omega²
    near(reaction.tensor.trace(),-.621075,2e-9);
    let p=endpoint_pressure(&s,&ModelEvaluation::zero(),Some(&c)).unwrap();
    near(p.scalar.unwrap(),.56/(3.*8000.),2e-9);
    assert_eq!(s.forces,f);
}
```

- [ ] Run `rtk cargo test -p modena-core --test stress_endpoint`; expect missing solver/observation. Add singular duplicate-edge and tangent-invalid errors; no NaN pressure or guessed zero. Empty constraints return zero exactly; missing model tensor gives total Missing/scalar None even though kinetic and constraint are known.
- [ ] Assemble block A and b per spec, row/column scale, symmetric pivoted Cholesky and permutation, reject rank/residual thresholds exactly. Build Fc=-Gᵀλ and dyad sum without writing system.forces.

```rust
b[e]=(v_j-v_i).norm2()+r.dot(f_j*(1.0/m_j)-f_i*(1.0/m_i));
// After A*lambda=b, reaction on j is -lambda[e]*r.
constraint=constraint.plus(SymTensor6::dyad(r,r).scaled(-lambda[e]));
```

- [ ] Add water analytic principal/nonprincipal Euler acceleration fixture from spec; solve pair tensions from atom force balance in test independently (different formulation, not production A builder), test all6 values and each Fc. Forced rotor trap tests b-force term; PBC/translation tests use exact analytic expectation. Tight tests1e-13 geometry; normal output preserves existing1e-10 tolerances.
- [ ] Run `rtk cargo test -p modena-core --test stress_endpoint`; expected analytic2e-9 and linear1e-11 gates pass. Commit `:sparkles: Add: endpoint mechanical constraint pressure`; review rank detection, units, force signs, volume and missing behavior.

### Task 7: VV exact interval estimator and compatible step result

**Files:** modify `src/integrator.rs`, `src/stress.rs`; create `tests/stress_vv.rs`.
**Consumes:** Tasks1/5/6 and exact spec §4a VV stages.
**Produces:** ConstraintSolverOptions, StepEvaluation, StepStress, IntervalPressure, BathImpulse and default/override Integrator::step_with_stress plus step_with_stress_options; one shared VV engine, legacy step unchanged return type.

- [ ] Add free-rotor test exercising the real production integrator and comparing both correct and intentionally wrong time conventions:

```rust
#[test]
fn vv_interval_uses_drift_kinetic_not_endpoint_kinetic() {
    let (mut s,c)=rotor(1.,3.); let h=.08/.7;
    let mut vv=VelocityVerlet::with_constraints(c);
    let mut nl=NeighborList::new(1.);
    let out=vv.step_with_stress_options(&mut s,&ZeroModel,&mut nl,h,&ModelEvaluation::zero(),ConstraintSolverOptions::TIGHT).unwrap();
    let stress=out.stress.unwrap();
    near(stress.interval.scalar.unwrap(),.56/(3.*8000.),5e-10);
    let wrong=stress.endpoint.kinetic.trace()+stress.interval.constraint_average.trace();
    let c=(1.-(.7*h).powi(2)).sqrt();
    near(wrong-.56,-.621075*(1.-c)/(1.+c),5e-10);
    assert!((wrong-.56).abs()>1e-4);
}
```

`ZeroModel` is a test-only ForceModel with no force additions, compute returns0, evaluate returns ModelEvaluation::zero; no synthetic pressure data. Define/import it in support here.

- [ ] Run `rtk cargo test -p modena-core --test stress_vv`; expect missing method or interval. Add exact ConstraintSolverOptions and per-call options API from spec §6; keep default1e-10/100, campaign TIGHT1e-13/100, unchanged public integrator struct shapes and legacy trajectories. Pass options through the shared engine, validate before mutation, and serialize actual options with the stress record. Instrument pre-position projection, corrected drift velocity, final force evaluation and final tangent projection exactly once. K_A must be sampled before final B/RATTLE.

```rust
let kinetic_average=kinetic_tensor(&drift_state);
let constraint_average=position.tensor_impulse.plus(final_projection.tensor_impulse).scaled(1.0/dt);
let model_average=average_model_tensors(&initial.tensor,&final_model.tensor);
```

Define private `average_model_tensors(&TensorAvailability,&TensorAvailability)->TensorAvailability` by combining complete values with .5; preserve missing reasons otherwise. `drift_state` above means sampling current System at that exact stage, not a required full clone.

- [ ] Add isolated-cluster full ΔG identity and exact finite rotor r1/u1/Ipos/Ifinal formulas from spec; all6 components, COM torque, translated and wrapped fixtures, masses2/5. Final cached force remains F(q1). No constraints reduces K_A+(S0+S1)/2. Legacy custom integrator default returns stress None.
- [ ] Run `rtk cargo test -p modena-core --test stress_vv` and existing VV/conservation tests. Compare old/new step trajectories and force call counts with actual model invocations (one post-drift evaluation/step); legacy outputs unchanged. Commit `:sparkles: Add: time-consistent Verlet pressure`; review discrete identity and negative control.

### Task 8: BAOAB complete stage accounting and conditional OU law

**Files:** modify `src/integrator.rs`, `src/stress.rs`; create `tests/stress_baoab.rs`; extend support.
**Consumes:** Task7 StepStress and Task5 projection reports; exact spec §4a BAOAB ordering and OU exclusion.
**Produces:** BAOAB step_with_stress override with separate bath ledger and interval pressure; legacy RNG checkpoint unchanged.

- [ ] Add actual gamma0 BAOAB rotor test with two h/2 maps, not one h map. Freeze h=.08/.7, compose e/t by rotations angle asin(.7*h/2) twice; assert final r/u to5e-10 and interval trace .56. Require six non-bath stage tags B0,A1Pre,A1Post,A2Pre,A2Post,B1 plus separate O reaction report; empty impulses allowed for truly zero stages.

```rust
#[test]
fn baoab_interval_excludes_ou_bath_exchange() {
    let (mut s,c)=rotor(1.,3.); let before=moment(&s);
    let mut b=LangevinBAOAB::with_constraints(1.,1.,1001,c);
    let h=.02;
    let out=b.step_with_stress_options(&mut s,&ZeroModel,&mut NeighborList::new(1.),h,&ModelEvaluation::zero(),ConstraintSolverOptions::TIGHT).unwrap();
    let st=out.stress.unwrap(); let bath=st.bath.unwrap();
    let lhs=moment(&s).plus(before.scaled(-1.));
    let rhs=complete(&st.interval.total).scaled(8000.*h)
        .plus(bath.unprojected_moment).plus(bath.reaction_moment);
    for (a,e) in lhs.0.into_iter().zip(rhs.0) { near(a,e,5e-10); }
}
```

`moment(&System)->SymTensor6` test helper computes Σunwrapped q⊙m v. Test isolated cluster does not wrap in its small step; PBC identity uses image tallies separately.

- [ ] Run `rtk cargo test -p modena-core --test stress_baoab`; expect absent BAOAB report/identity. Preserve number/order of StandardNormal samples and checkpoint state; sample K_A1/K_A2 before post-A projection, average .5; conservative endpoint average .5; divide sum of mechanical impulses by outer h once.

```rust
let kinetic_average=k_a1.plus(k_a2).scaled(.5);
let reaction_sum=mechanical_reports.iter().fold(SymTensor6::ZERO,|a,r|a.plus(r.tensor_impulse));
let constraint_average=reaction_sum.scaled(1./dt); // excludes BaoabO
```

- [ ] Add fixed-q conditional OU harness using actual production O-update/projector path (factor a private function used by BAOAB rather than copy code). At gamma*h=.2,kT1 draw65,536 cases; compare zero-mean reaction and full mass-matrix projected covariance from independently assembled dense projector against the frozen .02/4SE gate. Include unequal masses and water nonprincipal orientation; exact h0/gamma0 limits without altering RNG draws.
- [ ] Run `rtk cargo test -p modena-core --test stress_baoab`; expected discrete identities and actual conditional statistics pass. Existing seeded continuation tests must remain identical for mechanical state. Commit `:sparkles: Add: BAOAB stage-resolved pressure`; review bath exclusion is explicit and no ensemble-order claim yet.

### Task 9: Runner pressure observations, full caches, and controller failure semantics

**Files:** modify `src/runner.rs`, `src/thermo.rs` docs, `src/io/input.rs`, `src/config/validation.rs`, `src/barostat.rs`; create `tests/stress_runner.rs`; update final-base RunSpec literals/RunEvent exhaustive matches only as compiler requires.
**Consumes:** StepEvaluation/StepStress, endpoint_pressure, all shipped model overrides.
**Produces:** Simulation pressure getters, new RunEvent::Pressure, optional RunSpec.pressure_jsonl, Result-returning Barostat::apply and experimental apply_with_pressure exact spec §8. Legacy Energies/CSV meaning unchanged.

- [ ] Add failing guard test that proves old fabricated anisotropic scaling cannot mutate a state:

```rust
#[test]
fn anisotropic_scalar_call_fails_before_mutation() {
    let mut s=atoms(&[[1.,2.,3.]],&[1.]); let before=s.positions.clone();
    let old_box=s.box_size;
    let mut b=BerendsenBarostat::new(1.,1.,.1).anisotropic();
    let error=b.apply(&mut s,&Virial{w:9.},.01).unwrap_err();
    assert!(error.to_string().contains("tensor pressure"));
    assert_eq!(s.positions,before);assert_eq!(s.box_size,old_box);
}
```

- [ ] Run `rtk cargo test -p modena-core --test stress_runner`; expect unwrap_err on () compile failure. Change trait return type and all4 implementations, propagate in default apply_with_model, update known in-module calls549/567 and final-base callers. Add typed endpoint-diagonal direct experimental path; missing tensor/constrained/wrong-volume invalid input errors leave piston/cell/state unchanged. Keep runner config rejects.
- [ ] Runner owns `current_model:ModelEvaluation` plus private cache generation; calls step_with_stress with matching initial evaluation, carries final_model without rebuilding from e.virial. After velocity thermostat recompute K/reaction only; after a successful box move refresh complete model evaluation and mark old interval unavailable. A no-op barostat must not add a force call. Capture actual invocations in tests.

```rust
let out=self.integrator.step_with_stress(&mut self.system,&*self.model,&mut self.neighbors,dt,&self.current_model)?;
self.current_model=match out.stress.as_ref() {Some(s)=>s.final_model.clone(),None=>legacy_missing(out.energies)};
```

Private `legacy_missing(Energies)->ModelEvaluation` preserves E/W and records Missing(LegacyForceModel); it never invents tensor. Publish pressure after all output-time velocity edits, retain pre-controller interval provenance separately.

- [ ] End-to-end rotor config with thermostat velocity scale2 must multiply kinetic+constraint endpoint contribution by4 while preserving conservative F; legacy pressure field keeps prior meaning. Check scalar-only custom provider JSON has null new scalar/Missing cause; no zero. Check units: full energy tensor uses energy conversion, pressure total energy/volume; JSON all6 and step/time stamped. Add malformed output path failure.
- [ ] Run `rtk cargo test -p modena-core --test stress_runner`, `rtk cargo test -p modena-core --test runner_hardening`, `rtk cargo test -p modena-core --test ensembles`. Document downstream Barostat, RunSpec, RunEvent migrations. Commit `:sparkles: Add: pressure observations and controller capability guards`; review no controller is scientifically promoted.

### Task 10: Exact continuation with versioned stress sibling

**Files:** modify `src/io/restart.rs`, `src/runner.rs`; create `tests/stress_restart.rs`; preserve existing v2 compatibility fixture and all14 continuation fields.
**Consumes:** runner full cache, all pressure DTOs, exact spec §8 restart schemas.
**Produces:** StressCheckpoint/RestartStress, save_with_stress/load_with_stress; runner v3 exact resume; old save/load APIs/v2 wire unchanged.

- [ ] Add actual runtime continuation test, not JSON-key assertions: run seeded constrained BAOAB through17 steps, checkpoint, continue23; load new runner from checkpoint and run23; compare positions/velocities/forces/RNG/controller states and endpoint/interval/ledger values. Run same test with classic v2 fixture: mechanical continuation agrees, immediate interval Missing(LegacyRestart), first missing-S0 interval Missing, next complete interval valid. No forced scratch evaluation for live external sessions.

```rust
#[test]
fn legacy_v2_reader_rejects_v3_explicitly() {
    let cp=make_stress_checkpoint_after_rotor_step();
    let path=unique_test_path("b07-v3");
    save_with_stress(&path,&cp).unwrap();
    let restored=load_with_stress(&path).unwrap();
    assert_eq!(restored.mechanical.positions,cp.mechanical.positions);
    assert_eq!(restored.stress.as_ref().unwrap().model,cp.stress.as_ref().unwrap().model);
    assert!(load(&path).unwrap_err().to_string().contains("load_with_stress"));
    std::fs::remove_file(path).unwrap();
}
```

Test helper `make_stress_checkpoint_after_rotor_step()->StressCheckpoint` performs Task7 rotor step and copies actual mechanics/report into exact spec DTO, not hand-authored stress. `unique_test_path(&str)->PathBuf` combines temp dir/processID/atomic sequence and cleans files; define these in this test.

- [ ] Run `rtk cargo test -p modena-core --test stress_restart`; expect missing APIs. Implement v3 nested serializer with existing atomic write/fsync logic shared privately; parse v1/v2 through untouched readers and construct stress None. Validate estimator1, dt>0, volume/step/physical signature and finite shapes/scalars before restore; reject partial stage or mismatch.
- [ ] Maintain old `save` writing v2, old `load` accepting only its versions. Add actual old-reader compiled fixture test from B02/B01 baseline; do not replace with serde-key count. Add malformed tensor, future estimator, changed model pin/controller/order, unknown JSON fields, nonfinite numeric and atomic-preservation tests.
- [ ] Run `rtk cargo test -p modena-core --test stress_restart`, `rtk cargo test -p modena-core --test restart`, existing exact continuation tests. Commit `:sparkles: Add: versioned pressure continuation state`; review full state/RNG fidelity and explicit unavailable first interval.

### Task 11: Actual analytic dynamics and BAOAB distribution campaign

**Files:** create `examples/stress_validation.rs`, `examples/stress_validation/{fixtures,statistics,protocol}.rs`, extend `mod.rs`; create `tests/stress_statistics.rs`; create validation doc with frozen protocol before running outcomes.
**Consumes:** Tasks4–10; every fixture/tolerance/seed/duration/budget in spec §9.
**Produces:** bounded executable phases prepare/deterministic/ou/ensemble/report; typed protocol/results with actual measurements, never fabricated status/counters. This task includes real physical runs after source review and timing lease release.

Example-only exact records/functions:

```rust
pub enum CaseStatus { Prepared,Passed,Failed,Inconclusive,Unsupported,Pending,BudgetExceeded }
pub struct GateResult { pub name:String,pub observed:Vec<f64>,pub expected:Vec<f64>,
    pub absolute_error:Vec<f64>,pub tolerance:Vec<f64>,pub status:CaseStatus }
pub struct Sample { pub seed:u64,pub step:u64,pub time:f64,pub dt:f64,
    pub endpoint:[f64;6],pub interval:[f64;6],pub com_energy:f64,
    pub rotation_energy:f64,pub orientation_second:[f64;6] }
pub struct Estimate { pub mean:f64,pub se:f64,pub ci95:[f64;2],pub n_eff:f64,
    pub block_samples:usize,pub retained_blocks:usize }
pub fn estimate_replicates(samples:&[Vec<f64>],sample_interval:f64)->Result<Estimate,String>;
pub fn run_deterministic(out:&std::path::Path)->Result<Vec<GateResult>,String>;
pub fn run_ou(out:&std::path::Path)->Result<Vec<GateResult>,String>;
pub fn run_ensemble(out:&std::path::Path)->Result<Vec<GateResult>,String>;
```

Serialize these records with serde derive and finite validation. Use exact ValidationFixture/WorkLimits/ValidationProtocol/ExecutionProvenance/ValidationReport in spec §12, including actual source/binary/snapshot hashes and measured work/RSS. CLI validates a newly created output directory and immutable protocol hashes at every phase. Hash with host `shasum -a256` outside timed loops, require availability. Exit0 complete pass,1 internal/I/O,2 unsupported/budget,3 physical failure,4 inconclusive/pending. Preserve raw rows on any exit.

- [ ] Write failing analytic convergence integration test that executes VV hω=.08,.04,.02,.01 to Tω8 and uses exact continuum e(t)=cos(ωt)e0+sin(ωt)t0; compare all final q/v and normalized pressure finite-window ΔG/T against exact rotor and Euler water, ratios3–5. Test statistical estimator against deterministic eight-replicate nonidentical correlated time series with known replica means; require nonzero between-replica SE and increasing uncertainty for autocorrelation, not an IID false pass.

```rust
#[test]
fn actual_rotor_pressure_campaign_passes_analytic_gates() {
    let path=unique_test_path("b07-deterministic");std::fs::create_dir(&path).unwrap();
    let gates=run_deterministic(&path).unwrap();
    assert!(gates.iter().all(|g|matches!(g.status,CaseStatus::Passed)));
    assert!(gates.iter().any(|g|g.name=="vv_rotor_mixed_time_negative_control"
        &&g.observed[0]< -1e-4));
    std::fs::remove_dir_all(path).unwrap();
}
```

The negative-control gate's observed value is actual wrong numerator minus expected COM trace at hω=.08, matching spec's negative analytic number; it is not a static label/counter test. Add failing mocks for budget/file/provenance errors only for orchestration, never substitute mocks for physics gates.

- [ ] Run `rtk cargo test -p modena-core --test stress_statistics`; expect absent harness. Implement exact bounded fixture matrices, raw logs and analytic oracles. For principal-axis water use exact rigid rotation; for nonprincipal water integrate independent body Euler equations and quaternion kinematics with RK4 at h_min/32 and h_min/64, require their final state differences<=1e-10, and use the finer state as continuum reference. Verify analytic Euler Fc at each sampled state; never fit a state against itself.
- [ ] Implement gas statistics: ideal mean1/54, COM/rotational χ² distributions, orientation second moments, two-half autocorrelation blocking, eight-seed Student-t/bootstraps, predeclared duration extension and h² extrapolation per spec. Use actual disjoint-gas production integrator; no post hoc COM removal. Record endpoint/interval separately and all four h rows. Missing precision at max run is Inconclusive.
- [ ] Freeze manifests before physics, then execute sequentially under released lease:

```sh
rtk cargo run --release -p modena-core --example stress_validation -- prepare --out /private/tmp/modena-b07-validation-v1
rtk cargo run --release -p modena-core --example stress_validation -- deterministic --out /private/tmp/modena-b07-validation-v1
rtk cargo run --release -p modena-core --example stress_validation -- ou --out /private/tmp/modena-b07-validation-v1
rtk cargo run --release -p modena-core --example stress_validation -- ensemble --out /private/tmp/modena-b07-validation-v1
```

- [ ] Inspect real gate values, convergence ratios, covariance, distribution moments/CDFs, N_eff/CI halfwidth and work/RSS. Fix implementation errors with a fresh affected test cycle; if the projected BAOAB method fails physical sampling/order, report the concrete method deficiency to controller before changing the integrator/spec. A weak-order claim is withheld when slope is unresolved; actual distribution gates still must pass.
- [ ] Commit code and small protocol/result summary `:white_check_mark: Add: analytic and ensemble pressure validation`; large generated output remains outside git. Independent review checks preregistration, actual runs and uncertainty, not report prose alone.

### Task 12: Matched external validation, acceptance matrix, and final verification

**Files:** create `bench/b07_pressure_reference.py`, `bench/test_b07_pressure_reference.py`, `bench/fixtures/b07/`; update validation doc with actual outcome links/hashes; modify status promotion only after applicable gates pass. No physics edits inside the measurement phase.
**Consumes:** all preceding gates, frozen spec §10 exact external roster and protocol; external binary installation/authorization resolved by controller.
**Produces:** actual independent external parity/convergence/ensemble data, manifest+raw records, case-scoped validation status or explicit blocked/failed/inconclusive rows; complete B07 only if every required supported case passes.

- [ ] Write adapter tests that parse an actual small retained LAMMPS sample with pressure order/units, plus a synthetic asymmetric-column negative fixture; test process failure, timeout, missing packages and unexpected source hash. These test the adapter only. Write physics adapter test to load reference q/v/F and compare analytic rotor/water Fc and tensors on those exact snapshots; expected constraints from independent analytic formulas, not our production solver.

```python
import unittest
from b07_pressure_reference import parse_pressure_line

class PressureParserTests(unittest.TestCase):
    def test_pressure_columns_preserve_shear(self):
        row = parse_pressure_line("12 0.12 1 2 3 .4 .5 .6")
        self.assertEqual(row["step"], 12)
        self.assertEqual(row["pressure"], [1., 2., 3., .4, .5, .6])
```

Define `parse_pressure_line(line: str)->dict` as exactly eight fields step,time,xx,yy,zz,xy,xz,yz; reject wrong count/nonfinite. Use the exact eight-token line above. This is a parsing gate, not external acceptance.

- [ ] Run `rtk proxy python3 -m unittest discover -s bench -p 'test_b07_pressure_reference.py'`; expect missing adapter/test module. Implement CLI prepare/run/report with strict --binary,--protocol,--out arguments, subprocess timeout, CPU/thread caps, retained stdout/stderr/returncode and actual row counts. Never generate a success record when binary absent.
- [ ] Prepare immutable external pin and snapshot arrays; verify package support and conventions. Freeze all20 trajectory inputs plus<=4 reference preflight jobs, exact shifted-force LJ table, charge/exclusion maps, initial q/v, κ/rc and convergence ladder before running. Ensure force/energy/constraint convention is exactly matched; table/long-range mismatch blocks dependent runs.
- [ ] Run actual external campaign sequentially after analytical source review and lease release:

```sh
rtk proxy python3 bench/b07_pressure_reference.py prepare --binary "$B07_LAMMPS_BINARY" --protocol /private/tmp/modena-b07-validation-v1/protocol.json --out /private/tmp/modena-b07-reference-v1
rtk proxy python3 bench/b07_pressure_reference.py run --binary "$B07_LAMMPS_BINARY" --protocol /private/tmp/modena-b07-validation-v1/protocol.json --out /private/tmp/modena-b07-reference-v1
rtk proxy python3 bench/b07_pressure_reference.py report --binary "$B07_LAMMPS_BINARY" --protocol /private/tmp/modena-b07-validation-v1/protocol.json --out /private/tmp/modena-b07-reference-v1
```

`B07_LAMMPS_BINARY` is resolved/recorded by prepare, not guessed. If preparation requires a physical scheme/convention choice not covered by spec, stop that dependent gate and report exact evidence; do not silently choose a favorable reference.

- [ ] Match same-snapshot U/F/conservative all6, rotor/water analytic endpoint forces, NVE continuum-window pressure/state convergence and interacting ensemble means with required uncertainty. Finite-step LAMMPS constraint pressure is not asserted identical to our endpoint estimator. Reference NVT is the pinned LAMMPS rigid/nvt route and must first pass its two free-gas preflight jobs; actual LAMMPS NVE cases remain mandatory. No scheme switch after outcomes.
- [ ] Complete artifact report: all provider/integrator/constraint/status rows, raw files/hashes, failed/pending cases, measured work/RSS/ledger/tensor overhead, source migration and exact validated domain. Set validation only for tuples with every applicable gate passed; no global Validated switch or self-asserted restart promotion.
- [ ] Run required final checks once after last source change:

```sh
rtk cargo fmt --all -- --check
rtk cargo test -p modena-core
rtk cargo test -p modena-core --test spme_profile
rtk cargo clippy -p modena-core --all-targets -- -D warnings
rtk proxy python3 -m unittest discover -s bench -p 'test_b07_pressure_reference.py'
```

The full core suite includes all focused tests; the explicit B05 profile command records integration parity. Run existing repository-prescribed workspace checks found at adoption if stricter. No numerical campaign rerun unless source/protocol changes invalidate it.
- [ ] Commit documentation/status with `:memo: Document: validated pressure domains and compatibility`; keep status Diagnostic if any required row is pending. Final response states what passed, what failed/pends, exact source and artifact paths, public API migrations, costs and remaining scientific limits. No acceptance from compilation/JSON shape alone.

# B07 tensor pressure design — adopted contract

Prepared 2026-11-06 without implementation or measurements. Adopted 2026-11-10 for user-authorized B07 implementation, PR and merge.

Adoption source: integrated reviewed main `d6cd33ca3418b1617ab71fa2cf7416a24d27f405` (B04 PR59 and B05 PR60 complete). Branch `ricardo/b07-tensor-pressure`, existing isolated workspace `/Users/ricardoleal20/.codex/worktrees/modena-b04/modena`. Source hashes/toolchain are retained in `target/b07-evidence/adoption-20261110/source-epoch.json`. All B05 twelve Duration fields and compute_profiled signature exist unchanged in this base. No measurement lease is active. Earlier preparation epoch below remains historical. Current Conventional Commits rules override old plan gitmoji examples. All scientific gates and budgets below remain unchanged; adoption grants no scientific acceptance.
Future repository path: `docs/superpowers/specs/2026-11-06-b07-tensor-pressure-design.md`.
Future plan: `docs/superpowers/plans/2026-11-06-b07-tensor-pressure.md`.
Future branch: `ricardo/b07-tensor-pressure` (explicit user override).

## 1. Authority, epoch, and scope

Implements B07 in `docs/research/2026-11-05-core-improvements/phase3_analysis/ACCEPTANCE-BACKLOG.md:46`.
The source inspection used `/Users/ricardoleal20/.codex/worktrees/modena-b04/modena`, HEAD `7ea653641c813bc2f01476805db2990e0de4677a`.
At that historical preparation epoch B04/B05 were still owned by other work, and this draft was provisional.
The 2026-11-10 adoption above supersedes that dependency: final B04/B05 are integrated, their owners have released the source,
and the exact current hashes are recorded. B05's actual `compute_profiled` interface is present and checked on the adopted base.

The approved architectural direction is complete force tensors, actual constraint impulses, physically distinct endpoint/interval pressure,
and actual analytic, equilibrium-distribution, timestep and independently matched external gates. All cases in this specification are mandatory;
a missing reference, exceeded budget, failed gate, or incomplete provider leaves the affected case pending/failed, never accepted.
AlFCC belongs to B10 and is not the material fixture here. No broad triclinic SPME production work is a prerequisite.
No new Rust dependency or silent change to pair potentials, electrostatics parameters, FFT differentiation, thermostat random sequence,
legacy scalar `Virial`, `Energies`, or conservative force cache. Explicit API migrations are listed below.

This document supersedes the preparation report where it resolves alternatives: even-Nyquist convention, exact APIs, barostat failure,
serialization, fixture/precision/work budgets, and acceptance status. Scientific derivations in §§2–5 remain normative except for the
explicit decisions in §§6–11. The controller adopts the documents directly; sequential owner + independent review begins at Task 1.

## 2. Tensor convention and six strain derivatives

Use symmetric tensor order `[xx,yy,zz,xy,xz,yz]`, energy units, compressive-positive configurational virial S.
Define `A⊙B = (A⊗B+B⊗A)/2`, `K = Σ_i m_i v_i⊗v_i` (twice the usual kinetic-energy tensor).
For cell H (vectors are columns), fractional positions s fixed, deform `H(η)=(I+ηE)H`, `q(η)=(I+ηE)q`.
Then `dU/dη|0 = -S:E`; tensile-positive configurational stress is `σ = -S/V`.
For normal components use `E=e_a⊗e_a`; for shear use `E=(e_a⊗e_b+e_b⊗e_a)/2`.
These shear matrices give `dU/dη=-S_ab`, not `-2S_ab`. A matrix with both off-diagonals equal to 1 needs the factor 2.
`W=tr(S)`, `P=(K+S+S_c)/V`, `p=tr(P)/3`. Never infer any component of S from W alone.
Numerical FD: `S:E = -[U(+η)-U(-η)]/(2η)`; deform the cell AND every coordinate, keeping parameters fixed.
Do not project constrained geometry during this unconstrained Hamiltonian derivative test.
Hold κ, charges, mesh, integer reciprocal modes, topology/scales, and cutoff convention fixed.
Choose geometries away from real cutoffs, nearest-image ties, spline knots, and reciprocal cutoff membership changes.
Do not differentiate an adaptive tuning rule. Record it separately if it exists.

Contributions:
- Pair: choose `d=q_i-q_j` in the same periodic image as force `f_i`; add `d⊙f_i` once per half-list pair.
- Bond: same. Angle/dihedral: reconstruct the same connected local images used by that term;
  anchor atom j at zero and add `Σ_a (q_a-q_j)⊙f_a`. No wrapped global `Σq⊗f`.
  Angular energies can have zero trace while nonzero shear; include such fixtures.
- Ewald real and exclusion/scaling correction: outer product of the actual returned pair force and its image vector.
  `exclusion_correction` currently changes analytic branch outside `real_cutoff`; FD each smooth branch and avoid the seam.
  Validate fractional Coulomb scales as well as full exclusions, without double scaling.
- Direct reciprocal: for each existing weighted mode `U_k ∝ V^-1 exp(-k²/(4κ²)) |ρ_k|²/k²`,
  `k=2π H^-T n`, fixed s implies fixed `k·q`. Thus
  `S_k = U_k [I -(2/k²+1/(2κ²)) k⊗k]`.
  Derivation: `δlog V=tr ε`, `δk²=-2kk:ε`; insert these in `δlog U_k`.
  Trace is exactly `U_k[1-k²/(2κ²)]`, matching the current scalar formula.
  Preserve the existing ±mode multiplicity; do not introduce another half factor.
- Self energy has zero cell derivative at fixed κ. Uniform neutralizing background `U_bg=-C/V`
  contributes `S_bg=U_bg I`, trace `3U_bg`; charge-neutral-only tests would miss this.
- External ASE stress: preserve all six inputs, `S=-V_angstrom³ σ_eV/angstrom³ * EV_TO_KCAL_MOL/epsilon`.
  Keep the current force/energy unit mapping, symmetry and finiteness checks; do not subtract provider kinetic stress.
  Protocol must promise configurational tensile-positive stress. Unknown stress is unavailable, never zero.
- Composite: add tensors only when every nonempty contribution supports them; partial support is explicitly incomplete.

### SPME shear is a real design gate

Current `SpmeForceModel::compute` and reciprocal force helper reject triclinic cells; `InfluenceKey`
uses κ and three box edges. Merely adding six cached factors cannot establish six strain derivatives.
At fixed fractional coordinates and integer mesh, deposited B-spline charges and their DFT are fixed.
For the existing reciprocal mesh energy, differentiate only V and physical k in the influence factor:
`S_mesh = Σ_k U_k[I-(2/k²+1/(2κ²))kk]`, with the **same explicitly defined discrete mode convention**.
Minimum bounded route: retain orthorhombic production support; implement an independent test-only energy oracle
for arbitrary H using fractional charge spread, a small direct DFT, and physical reciprocal vectors H^-T n.
Use its six affine derivatives at orthorhombic H; add real/self/background/exclusion energy independently.
The oracle must reproduce production energy at H0 before its derivatives can certify anything.
Odd meshes first remove Nyquist ambiguity. Even meshes need an explicit derivation of Nyquist/alias representatives
and conjugate symmetry under shear; compare the resulting derivative to production stress for mixed Nyquist modes.
Do not silently choose an extension whose shear derivative differs from the intended mesh Hamiltonian.
If the actual chosen discrete energy cannot be differentiated independently, B07 SPME shear remains unpassed.
An alternative is genuine triclinic production support, requiring fractional spread/gather, Cartesian gradient H^-T,
full-cell cache keys, geometry-safe neighbors and corresponding tests; that is larger scope, not prerequisite by fiat.
Coordinate the oracle and convention with B05; isotropic scaling tests or trace comparisons alone do not close this gate.

## 3. Pair impulse ledger with exact stage geometry

For constraint e=(i,j), choose `r_e=q_j-q_i` in a validated bonded image. Let J_e be impulse on j, -J_e on i.
Record `I_e = r_e⊙J_e` (energy*time). Divide the accumulated step ledger by the OUTER dt only once.
A position correction over drift τ means `J_i=m_i Δq_i/τ`, not `m_i Δq_i/τ²`.
The latter has force units and lacks the operator half-kick factor used by RATTLE.
Accumulate every Newton/GS update, including repeated solves for coupled constraints; do not overwrite the last multiplier.
The loose SHAKE code has `Δq_i=+g r_ref/m_i`, so `J_e=-g r_ref/τ`, `I_e=-g r_ref⊗r_ref/τ`.
Loose velocity RATTLE has `Δv_i=+μ r/m_i`, hence `J_e=-μ r`, `I_e=-μ r⊗r`.
Water's current signs differ: its code uses `Δq=M^-1 G_ref^T lam` with G row (-r,+r).
Thus water `J_e=lam_e r_ref/τ` and velocity `J_e=mu_e r`; add positive lam/τ or mu times `rr`.
This sign mapping follows the actual updates, not variable names. Test reconstructed per-atom impulses against mass*delta-v.
Standalone `settle` initialization is a geometry repair, with no physical elapsed time: never invent a pressure impulse for it.

Freeze a consistent periodic bond image over a stage, selected before prediction, and carry its continuous vector through projection.
Wrap only afterward. Do not re-minimum-image a large iterative correction into another branch.
Reject ambiguous half-cell bonds, singular/duplicate constraints, inconsistent cycles, and stage motion that loses the chosen image.
For water, verify `r_H1H2=r_OH2-r_OH1`; three independent minimum-image calls can otherwise select an inconsistent triangle.
For general contractible graphs, build image shifts by graph traversal and validate cycle closure.
Pair-edge tallies avoid origin dependence; an unwrapped molecule sum is an independent test, not a replacement for edge tracking.
Networks with nonzero periodic winding need a separately derived image convention; reject unsupported topology instead of promising it.
Keep momentum impulses and full 3×3 torque checks in test diagnostics before symmetric compression.

## 4. Two pressure observables; do not combine their time conventions

### 4a. Interval operator estimator, actual impulses

For each constrained drift, reinterpret the implemented position projection exactly as a pre-drift impulse
at q_start, then a straight drift at `v_A = v_before + Δq/τ`; optional tangent correction is a post-drift impulse.
This is an algebraic equivalence of `project_drift`, since `q_end=q_start+τv_A` before wrapping.
For `G=Σ q_i⊙p_i`, a fixed-q kick gives `ΔG=Σq_i⊙J_i`; straight drift gives `ΔG=τK(v_A)`.
Consequently the following ledger satisfies the discrete virial identity for an isolated unwrapped cluster.
For periodic interactions use the force/constraint image tallies above; do not assert global wrapped G telescopes.

VV, outer dt=h:
- B0: conservative impulse hF0/2 at q0; model contribution hS0/2.
- SHAKE drift τ=h: position impulse at q0, then drift K_A at the corrected velocity.
- B1: conservative impulse hF1/2 at q1; model contribution hS1/2.
- final RATTLE: velocity impulse at q1.
- `P_interval V = K_A +(S0+S1)/2 +(I_pos+I_final)/h`.
  If expressed as effective endpoint constraint forces, each is `J/(h/2)` and receives weight 1/2.
  Thus position effective force is `2mΔq/h²`, not `mΔq/h²`; final is `2mΔv/h`.

BAOAB, outer dt=h, τ=h/2:
1. B0 at q0 then tangent projection: tally `hS0/2` and `I_B0`.
2. A1: tally position impulse `I_A1pre` at q0, drift `τ K_A1`, tangent `I_A1post` at q_mid.
3. O at q_mid: store unconstrained OU impulse and following reaction `I_Oprojection` separately.
4. A2: tally `I_A2pre` at q_mid, drift `τ K_A2`, tangent `I_A2post` at q1.
5. B1 at q1 then tangent projection: tally `hS1/2` and `I_B1`.
`P_interval V = (K_A1+K_A2)/2 +(S0+S1)/2`
`                  +(I_B0+I_A1pre+I_A1post+I_A2pre+I_A2post+I_B1)/h`.
For each A pair, effective impulse force is J/(τ/2) and its outer-step weight is τ/(2h)=1/4.
Each B projection has force J/(h/2), weight 1/2. These weights are not all identical.
Never divide every projection by τ and then sum as full forces.

O is a thermostat exchange, not an interatomic mechanical pressure term. Retain its complete impulse ledger for the
identity `ΔG = h V P_interval + Σq_mid⊙J_Ototal` (isolated cluster), but exclude it from P_interval.
At fixed q the momentum tangent projector is `Π=I-G^T(GM^-1G^T)^-1 GM^-1`.
Its projected Gaussian covariance is `kT[ M-G^T(GM^-1G^T)^-1 G ]`. With tangent p, OU projection reaction is
`J_Oprojection=(Π-I)[(c-1)p+sqrt(1-c²)M^1/2 sqrt(kT) ξ]`.
Its deterministic part is zero and its conditional mean is zero. Including it would add a noisy bath reaction estimator,
not endpoint mechanical pressure. This statement relies on scalar uniform γ and the current mass-weighted OU covariance.
A general thermostat requires a fresh derivation. Verify both conditional zero mean and projected covariance numerically later.
The operator identity proves bookkeeping, not the invariant distribution or finite-dt thermodynamic accuracy.
A weak-error order claim for this precise projected BAOAB implementation remains a validation gate, not a citation shortcut.

### 4b. Endpoint mechanical estimator, independently solved constraint acceleration

Use the final tangent velocities at the same q as the cached final conservative forces.
Define `C_e=(|r_e|²-d_e²)/2`, G row (-r_e,+r_e), M diagonal masses, `A=GM^-1G^T`.
Set `F_c=-G^T λ`. Differentiating `C=0` twice gives
`0=|v_j-v_i|² + [GM^-1F]_e - [Aλ]_e`.
Thus solve `Aλ=b`, `b_e=|v_j-v_i|² + r_e·(F_j/m_j-F_i/m_i)`;
then `S_c_endpoint=-Σ_e λ_e r_e⊗r_e` and `P_endpoint=[K(v_end)+S_model(q_end)+S_c_endpoint]/V`.
Use a connected-block rank-revealing solve with finite/rank/residual checks; disjoint water gets a 3×3 block.
This is not an inferred force from geometry residuals, and it does not replace the impulse ledger.
For Langevin, this is the deterministic mechanical reaction; friction is tangent and projected noise is bath exchange.
It supplies an unambiguous endpoint quantity aligned with `Energies::from`; at finite dt sampled states still have integration bias.
If the runner applies a subsequent thermostat or velocity edit, recompute K and this acceleration solve after that edit.
Never add interval S_c to endpoint K, and never mix S0/S1 quadrature with a single-endpoint kinetic tensor.
Keep external barostat use restricted until its specific estimator/coupling scheme is validated; replacing W/3 alone is insufficient.

## 5. Analytic oracles that can falsify the implementation

### Unequal-mass free rigid rotor: continuum AND finite-dt identity

Set masses m1=1,m2=3 (also test 2,5), total M, reduced μ=m1*m2/M, length d=1.3.
Choose orthonormal e,t oblique to box axes, relative position r=d e, relative velocity u=dω t;
`q1=R-m2*r/M`, `q2=R+m1*r/M`, `v1=V-m2*u/M`, `v2=V+m1*u/M`.
Use ω=0.7, nonzero COM V, no external forces, large box; translate and split the pair across each face/corner.
Continuum: λ=μ|u|²/d², `S_c=-μω² r⊗r`, `K=M V⊗V+μu⊗u`.
Hence `tr(K+S_c)=M|V|²` exactly; the rotational tensor need not vanish instantaneously.
Orientation/time average of rotational K+S_c is zero; finite-window tensor oracle is ΔG_internal/T.

Exact single free RATTLE drift of duration h, `a=ωh`, `c=sqrt(1-a²)`, |a|<1:
`e1=c e+a t`, `r1=d e1`, `u_A=d[(c-1)e+a t]/h`, `u1=dω[-a e+c t]`.
`I_pos=μd²(c-1) e⊗e/h`, `I_final=-μd²(1-c)e1⊗e1/h`.
Thus `tr(K_A,internal)=2μd²(1-c)/h²` and `tr((I_pos+I_final)/h)` is its exact negative.
The correct interval scalar has no spurious rotational pressure at any solvable h.
Mixing endpoint K instead produces the explicit numerator error
`μd²ω² - 2μd²(1-c)/h² = -μd²ω²(1-c)/(1+c) = O(h²)`.
Require this negative-control mismatch to appear when deliberately using the wrong estimator in the test oracle.
For γ=0 BAOAB, two A substeps use h/2, not one h rotor map; test each substep then the composition.
For forced rotor add a weak smooth pair/environment force and compare the independent Aλ=b oracle at every endpoint.

### Rigid water / SETTLE-comparable oracle

Use O,H1,H2 masses 15.9994,1.008,1.008; OH=1, HH=2 sin(109.47°/2), then isotope masses 16,1,2.
Build the planar triangle, subtract its true mass COM, call internal coordinates ρ_i and `C=Σm_iρ_iρ_i^T`.
Rotate the whole fixture by a generic orthogonal matrix; choose Ω along its plane normal (a principal inertia axis).
With no external force, `v_i=V+Ω×ρ_i`, `a_i=Ω×(Ω×ρ_i)=-ω²ρ_i`.
Analytic `F_c,i=-m_iω²ρ_i`, `S_c=-ω² C`, and `K_internal=Σm_i(Ω×ρ_i)(Ω×ρ_i)^T`.
Their traces cancel exactly; endpoint scalar numerator is again M|V|².
Solve the three pair tensions independently from these analytic atom forces; compare production Aλ=b and all six components.
For a general nonprincipal Ω, use Euler's equation `I Ωdot=-(Ω×IΩ)` and
`a_i=Ωdot×ρ_i+Ω×(Ω×ρ_i)` as a further independent force oracle.
Compare the existing `waters` Newton path with three explicit bond constraints: geometry, per-atom impulses, tensor, COM, torque.
Do not call matching these two numerical implementations independent physical validation: the analytic force oracle is essential.
For finite-step Newton water, check exact ΔG ledger and dt convergence toward analytic rigid-body motion;
there is no promise that current RATTLE water follows an exact analytical SETTLE trajectory at finite dt.


## 6. Exact core contracts and file map

New `crates/modena-core/src/stress.rs` owns tensor algebra, availability, and pressure records; export as `pub mod stress` in `lib.rs`.
All public record fields below are exact. Stress DTOs derive Debug/Clone/PartialEq/Serialize/Deserialize; derive Copy only for fixed scalar records.
Use serde snake_case adjacently tagged enums (`tag="kind",content="value"`), `deny_unknown_fields` on records, and validate all deserialized floats as finite.
Runtime-only ConstraintGeometry and StepEvaluation need Debug/Clone only; existing Energies stays unchanged. Add transparent serde derives to
Vec3 for impulse DTOs. Do not serialize SimBox caches; Checkpoint remains the existing validated mechanical serializer. Tensors are energy
numerators until divided by positive finite volume; pressure records use reduced pressure units. `Vec3` is the existing `[f64;3]` wrapper.

```rust
pub struct SymTensor6(pub [f64; 6]); // xx yy zz xy xz yz
impl SymTensor6 {
    pub const ZERO: Self = Self([0.0; 6]);
    pub fn dyad(a: Vec3, b: Vec3) -> Self; // symmetric a⊙b
    pub fn trace(self) -> f64;
    pub fn scaled(self, x: f64) -> Self;
    pub fn plus(self, rhs: Self) -> Self;
    pub fn finite(self) -> bool;
}
pub enum MissingTensorReason {
    LegacyForceModel, UnsupportedBackend, MissingInitialEvaluation,
    PostStepTransform, LegacyRestart, ProviderPropertyUnavailable,
}
pub struct MissingTensor { pub component: String, pub reason: MissingTensorReason }
pub enum TensorAvailability { Complete(SymTensor6), Missing(Vec<MissingTensor>) }
pub struct ModelEvaluation {
    pub energy: f64, pub scalar_virial: f64, pub tensor: TensorAvailability,
}
impl ModelEvaluation {
    pub fn zero() -> Self;
    pub fn combine(self, other: Self) -> Self;
}
pub enum ValidationStatus { Diagnostic, Validated { pub_record: String } }
pub struct PressureObservation {
    pub kinetic: SymTensor6,
    pub model: TensorAvailability,
    pub constraint: SymTensor6,
    pub total: TensorAvailability, // (kinetic+model+constraint)/volume
    pub scalar: Option<f64>,      // trace(total)/3; None on Missing
    pub volume: f64,
    pub validation: ValidationStatus,
}
pub struct IntervalPressure {
    pub kinetic_average: SymTensor6,
    pub model_average: TensorAvailability,
    pub constraint_average: SymTensor6,
    pub total: TensorAvailability,
    pub scalar: Option<f64>,
    pub volume: f64, pub dt: f64,
    pub validation: ValidationStatus,
}
pub struct ConstraintSolverOptions { pub tolerance:f64, pub max_iterations:usize }
impl ConstraintSolverOptions {
    pub const DEFAULT:Self=Self { tolerance:1e-10,max_iterations:100 };
    pub const TIGHT:Self=Self { tolerance:1e-13,max_iterations:100 };
}
pub struct StepStress {
    pub solver:ConstraintSolverOptions,
    pub final_model: ModelEvaluation,
    pub endpoint: PressureObservation,
    pub interval: IntervalPressure,
    pub stages: Vec<ProjectionReport>,
    pub bath: Option<BathImpulse>,
}
pub struct StepEvaluation { pub energies: Energies, pub stress: Option<StepStress> }
```

`combine` adds energy/scalar independently, adds complete tensors, propagates all missing reasons with component provenance, and never
constructs a tensor from a scalar. `zero()` is Complete(ZERO). An unknown legacy provider on an empty system is still Missing;
a built-in known empty contribution is Complete(ZERO). Current composite is pair + optional bonded/long-range, not an arbitrary boxed list;
extend that existing structure only, and unit-test the general combine rule separately.

In `force_model.rs`, retain `Virial { pub w:f64 }` and `compute` unchanged. Add a nonbreaking trait default:

```rust
fn evaluate(&self, system: &mut System, neighbors: &mut NeighborList)
    -> Result<ModelEvaluation> {
    let mut w = Virial::default();
    let energy = self.compute(system, &mut w, neighbors)?;
    Ok(ModelEvaluation { energy, scalar_virial:w.w,
        tensor:TensorAvailability::Missing(vec![MissingTensor {
            component:std::any::type_name::<Self>().into(),
            reason:MissingTensorReason::LegacyForceModel }]) })
}
```

`evaluate` adds only this contribution to existing forces and returns only this contribution's E/W/S; never zero caller buffers.
Each shipped CPU model implements one shared private tensor pass; `compute` adds returned `scalar_virial` to caller W, and `evaluate`
returns the record. Keep legacy scalar arithmetic/order where practical and require scalar parity tolerance, not bitwise tensor reduction.
Pair partitions gain one SymTensor6 (48 bytes/partition), not extra 3N scratch buffers. Preserve B02 lease/reuse and B03 packed APIs.
`BondedForceModel`, `EwaldForceModel`, `SpmeForceModel`, `LongRangeForceModel`, `CompositeForceModel`, `AseEmtForceModel` override.
GPU/MPI providers remain explicitly tensor-unsupported via the legacy default `Missing(LegacyForceModel)` (or `Missing(UnsupportedBackend)` when overridden) until independent tensor implementations exist; their scalar behavior is retained.
Do not claim validated pressure for those backends. The separately implemented NoseHooverChain integrator remains stress None/Missing/Diagnostic: the derived supported integrator domain is VV and BAOAB. Its thermostat ledger needs a separate derivation and is not implied by this work. No unused new generic composite container.

External provider: preserve inherent `AseEmtForceModel::evaluate(&System)->Result<ProviderEvaluation>` and qualify the new trait method
as `ForceModel::evaluate(&model,&mut system,&mut neighbors)` to avoid the same-name overload ambiguity. Keep fixed ASE protocol wire
shape and pin checks; stress is required configurational tensile-positive. Missing/malformed/nonfinite/asymmetric provider stress remains
an existing protocol error; the generic missing-capability variant is for providers that actually declare optional stress. `converted`
returns `(f64,Vec<f64>,SymTensor6)`; all six values map with `-V_A^3*EV_TO_KCAL_MOL/epsilon`, with trace as legacy W. No kinetic subtraction.
Keep staged additive overflow validation and session invalidation. Test a synthetic symmetric stress with unequal diagonals/nonzero shears.

New `constraint_stress.rs` owns periodic stage geometry, reports, and endpoint solve. Existing constraint updates remain in `constraints.rs`.

```rust
pub enum ProjectionStage { VvPosition, VvFinal, BaoabB0, BaoabA1Pre, BaoabA1Post,
    BaoabO, BaoabA2Pre, BaoabA2Post, BaoabB1 }
pub struct ConstraintEdge { pub i:usize, pub j:usize, pub distance:f64,
    pub image_shift:[i64;3], pub reference:Vec3 }
pub struct ConstraintGeometry { pub edges:Vec<ConstraintEdge>, pub cell:SimBox }
pub struct ConstraintImpulse {
    pub edge:usize, pub lever:Vec3, pub impulse_on_j:Vec3,
}
pub struct ProjectionReport {
    pub stage:ProjectionStage, pub impulses:Vec<ConstraintImpulse>,
    pub tensor_impulse:SymTensor6, pub iterations:usize,
    pub max_position_residual:f64, pub max_velocity_residual:f64,
}
pub struct BathImpulse {
    pub unprojected:Vec<Vec3>, pub reaction:Vec<Vec3>,
    pub unprojected_moment:SymTensor6, pub reaction_moment:SymTensor6,
}
pub struct EndpointReaction {
    pub multipliers:Vec<f64>, pub forces:Vec<Vec3>, pub tensor:SymTensor6,
    pub relative_linear_residual:f64,
}
pub fn constraint_geometry(system:&System, constraints:&Constraints)
    -> Result<ConstraintGeometry>;
pub fn endpoint_reaction(system:&System, geometry:&ConstraintGeometry)
    -> Result<EndpointReaction>;
pub fn kinetic_tensor(system:&System) -> SymTensor6;
pub fn endpoint_pressure(system:&System, model:&ModelEvaluation,
    constraints:Option<&Constraints>) -> Result<PressureObservation>;
```

Cell vectors follow `SimBox`'s column convention. `image_shift` means r=q_j-q_i+H*shift. `constraint_geometry` enumerates bonds then
water OH1/OH2/H1H2 (stable indices), validates undirected duplicate edges and cycle consistency, and chooses images before prediction.
Use B04 validated closest-image result; reject ties at squared-distance difference <=64 eps*max(1,d²,|H|²), noncontractible winding,
invalid distance/mass, overlapping singular graphs, and incompatible water cycles. Reference vectors are carried continuously during a
stage; updates never re-minimum-image them. Keep unwrapped prediction/correction until projection ends, then wrap once. Assert positive
reference/current dot and norm of relative stage displacement < half the shortest lattice vector; reject otherwise as image/branch loss.
Pure rigid translation of the whole component is allowed; apply the displacement criterion to edges, not atom COM travel.

Add exact methods in Constraints, leaving original signatures as wrappers around the same solver with recording disabled:

```rust
pub fn shake_report(&self, system:&mut System, geometry:&ConstraintGeometry,
    interval:f64, stage:ProjectionStage, tol:f64, max_iter:usize)
    -> Result<ProjectionReport>;
pub fn rattle_report(&self, system:&mut System, geometry:&ConstraintGeometry,
    stage:ProjectionStage, tol:f64, max_iter:usize) -> Result<ProjectionReport>;
```

`shake_report` is position-only, called after an unwrapped drift, uses geometry.reference from the start, and divides each mass-correction
by interval for impulse units; the integrator applies total Δq/interval to velocity exactly once. `rattle_report` changes only velocity. Build its geometry at the actual fixed-q stage (after the relevant drift), while a position report retains the pre-drift geometry. A stage cannot reuse a prior stage's lever accidentally.
Water `settle` initialization has no elapsed time and produces no pressure report. Its dynamics path is Newton/M-SHAKE with reference
normals; report the actual per-edge increments on every iteration. Reports retain each update (not just net multiplier); they are local
to a step, not indefinitely appended. Store full lever and impulse vectors so tests can reconstruct unsymmetrized torque. Both report
APIs validate before mutation and restore positions/velocities on solver failure; discard incomplete reports. Wrapper initialization keeps
its documented semantics. Step failure is terminal and emits no completed pressure/checkpoint; no claim of full integrator rollback.

Endpoint solve: build connected graph blocks; A_ef=Σ_i G_ei·G_fi/m_i, b as §4b. Use symmetric pivoted Cholesky on row/column-scaled A;
reject pivot <=1e-12 of maximum initial diagonal, nonfinite arithmetic, or normalized infinity residual >1e-11. Require geometry residual
<=1e-10 and tangent residual <=1e-10 for normal output; deterministic tight-gate fixtures use <=1e-13. Return constraint_error/Physics
with block edge IDs and residual. Solve once per block; disjoint waters are 3x3. Zero constraints returns exact zero without a solve.
Never append Fc to system.forces. Dense general block cost O(E_block³), memory O(E_block²), checked allocations; record the cost.

Integrator extension (default preserves external trait implementers):

```rust
fn step_with_stress(&mut self, system:&mut System, model:&dyn ForceModel,
    neighbors:&mut NeighborList, dt:f64, initial:&ModelEvaluation)
    -> Result<StepEvaluation> {
    Ok(StepEvaluation { energies:self.step(system,model,neighbors,dt)?, stress:None })
}
```

Add `Integrator::step_with_stress_options(&mut self,system:&mut System,model:&dyn ForceModel,neighbors:&mut NeighborList,dt:f64,initial:&ModelEvaluation,solver:ConstraintSolverOptions)->Result<StepEvaluation>`. Its default accepts DEFAULT and delegates to step_with_stress; a nondefault solver request on an unsupported integrator returns Config. Built-in VV/BAOAB override through one private shared step engine; legacy `step` and ordinary step_with_stress use DEFAULT, discard/report respectively and retain old Energies fields. Campaigns call the options method with TIGHT. Validate positive finite tolerance and max_iterations>0 before mutation. This keeps existing default trajectories and exposes tighter gates without changing public integrator struct layouts. Store actual options in StepStress and provenance; checkpoint options must match upon continuation.
`initial` must correspond to the cached input F/q/cell (runner-owned invariant); no extra force call to reconstruct S0 in the hot loop.
The public method requires a matching initial record; it cannot infer caller cache correctness. Runner stores a private position/cell generation with its cache and validates before stepping.
Use exact stage ordering and weights in §4. K_A uses corrected straight-drift velocity before post-A tangent projection, never final K.
OU report includes all per-atom raw bath impulses and O-projection edge impulses separately; mechanical interval excludes both.
At gamma=0 the bath is zero; preserve current RNG consumption even when its amplitude is zero for continuation compatibility.
`dt<=0`, nonfinite dt/volume, invalid stage/interval, or stale initial cache returns Config/Physics before publishing a record.
`StepStress` is pre-controller; runner stores/labels any later transformation separately.

`Energies` stays six scalar fields and `pressure` remains the old configurational-only diagnostic for source compatibility.
Do not silently change the meaning of existing CSV columns. Add `Simulation::pressure_observation(&self)->Option<&PressureObservation>`
and `Simulation::interval_pressure(&self)->Option<&IntervalPressure>`, new `RunEvent::Pressure { step:usize,time:f64,endpoint:PressureObservation,
interval:Option<IntervalPressure> }`, and optional companion output `pressure.jsonl` when `RunSpec.pressure_jsonl:Option<String>` is set.
Document the additive RunEvent variant's exhaustive-match migration and RunSpec struct-literal migration (new serde-defaulted field is None). Endpoint is recomputed after any thermostat velocity edit;
interval remains explicitly pre-controller and does not become the thermodynamic pressure of the post-controller state. On box/controller
move mark interval Missing(PostStepTransform); reevaluate final model on the updated cell using the same model, keeping full cache.
An external custom integrator returning stress None produces Missing diagnostics, never a fabricated tensor.

## 7. Even-Nyquist SPME ruling and independent oracle

The current khat maps Nyquist to +N/2, which does not define a unique continuous-cell shear Hamiltonian. Freeze this **new explicit extension**:
for each discrete index m, let R(m) contain all independent choices ±N_a/2 on every even-N Nyquist axis and the ordinary signed index
on other axes. Give its 2^r representatives equal weight. At arbitrary H define

`G_m(H)=2π/[det(H)*b_m²] * mean_{n in R(m)} exp(-|2πH^-Tn|²/(4κ²))/|2πH^-Tn|²`.

Mask the all-zero index. `b_m=product_a (4+2cos(2πm_a/N_a))/6`. Let Q_m be the unnormalized direct DFT of independently spread real
order-4 B-spline charges; `U_rec=Σ_m G_m |Q_m|²` (normalization matches current alpha=2π/V). This is a declared Nyquist-even product-cosine
extension of the mesh kernel, not a claim of a unique physical continuation. G_m=G_{-m} for arbitrary H by the lift bijection n→-n;
therefore a real charge grid has real energy/potential and conjugate symmetry. At orthorhombic H all lift |k|² agree, so energy, Cartesian
coordinate forces at fixed H, and isotropic trace agree with the old model exactly apart from reduction roundoff. Shear terms touching
any Nyquist axis average to zero; cross terms between two non-Nyquist axes survive, even if the third axis is Nyquist.

Production remains orthorhombic. Retain the old scalar vfac to preserve its arithmetic and add48bytes per mesh mode for the six strain factors (plus allocator metadata); cache six strain factors: same scalar diagonal formula; for ab off diagonal use
`-(2/k²+1/(2κ²))*k_a*k_b` unless a or b is Nyquist, then zero. Accumulate `U_m*factor`.
Do not zero whole modes, all shears on a mode, mesh forces, or non-Nyquist mixed terms. Reuse `InfluenceKey` (κ and three edges plus
existing workspace mesh identity); arbitrary-H test oracle has no shared production cache.

New test/example support `stress_validation/mesh_oracle.rs` uses explicit cardinal B-spline polynomials and direct complex DFT loops,
small fixed grids only. It may use Vec3/basic math but may not call production spread, influence, FFT, reciprocal helpers or their tensor
formula. It obtains derivatives by energy only; full oracle adds independent real/image/exclusion/self/background energy. Freeze the grid,
κ and mode set for all ±strains. Odd grids [5,5,5] and [5,9,5]; even [4,6,8] and [6,6,6]; max225 grid points for the full matrix. Per deformed energy direct DFT <=50,625 complex terms, ±6 strains*4 eps plus H0 <=2,480,625 terms per fixture.
Add public explicit-grid `spme_reciprocal_tensor(kappa:f64,box_size:&SimBox,positions:&[f64],charges:&[f64],forces:&mut[f64],mesh:[usize;3])->Result<(f64,SymTensor6)>` alongside the existing standalone helpers; it uses the same kernel and still rejects triclinic production. Existing helpers discard tensor/return trace. Do not change `mesh_for_box` or automatic mesh behavior.

Mixed-Nyquist synthetic grid has Q only at (Nx/2,Ny/2,1) and its conjugate: expected Sxy=Sxz=Syz=0. A second mode (Nx/2,1,1) has
Sxy=Sxz=0 and Syz nonzero per the formula. Compare actual DFT energy derivatives, not the convention function against itself.
Verify production H0 U/F/W first (coordinate FD for F), six independent affine derivatives next, and force/energy/conjugacy for arbitrary-H
oracle itself. Full particle meshes also cover net charge, fractional pair scales, cold/warm cache, integer wraps, non-grid translations,
and fixed-grid convergence toward independently converged Ewald. No exact continuous-translation invariance claim for SPME.

Strict independent H0 scalar verification supplements the fixed-epsilon derivative tests: evaluate the independent reciprocal energy under uniform complex cell scale lambda=1+i*1e-20 and obtain W=-Im(U)/1e-20. Real lambda=.97/1.03 probes must match the independent real-cell energy at1e-11. Full-model nonmesh contributions use first-order dual differentiation of their independent scalar energies, including elementary erfc differentiation, and match real energies at lambda=.97,1,1.03. Neither path calls production virial/tensor factors. Independent full and reciprocal U/W1e-11, coordinate FD2e-7, all6 fixed-epsilon strain derivatives2e-6, and legacy scalar/tensor-trace5e-12 are distinct gates; finite differences are not falsely required to achieve the stricter algebraic identity tolerance.

The auxiliary small automatic-grid full-model derivative fixture uses kspace_cutoff1.2 (mesh[4,4,4]), fixed by geometry-only knot checks before energy outcomes. The initial1.5/[4,5,5] preparation is retained as rejected because outside y=.4 lies on a spline knot. This auxiliary low-resolution identity control changes neither the four mandatory explicit grids nor the8³/16³/32³ accuracy ladder, and it grants no precision acceptance.

B05 compatibility: preserve `SpmeForceModel::compute_profiled(&self,&mut System,&mut Virial,&mut NeighborList)->Result<(f64,SpmeProfile)>`.
Twelve Duration fields remain exactly `neighbor,real_space,pair_correction,workspace_setup,influence_setup,spread,fft_forward,influence,
fft_inverse,gather,self_background,total`. Add `evaluate_profiled(&self,&mut System,&mut NeighborList)->Result<(ModelEvaluation,SpmeProfile)>`
using the same compile-time observer path. Tensor work is charged to existing influence/real/correction phases; no thirteenth phase,
no observer timers on ordinary calls, no recomputation, no fake zero durations. Final B05 synchronization is required before implementation.

## 8. Controllers and compatibility migrations

Maintain runner configuration rejection of anisotropic, legacy piston, and constrained box scaling. B07 does not validate B08 or NPT coupling.
Isotropic unconstrained existing behavior remains; diagnostic pressure status alone never grants controller permission.

Explicit source-breaking exception: change `Barostat::apply(&mut self,&mut System,&Virial,f64)` return type from `()` to `Result<()>`.
Isotropic paths return Ok after their same arithmetic; anisotropic scalar calls return Config("anisotropic coupling requires tensor pressure")
BEFORE any cell, positions, velocities or piston mutation. No panic, silent no-op, or W/3 path remains. Update `NoneBarostat`, `BerendsenBarostat`,
`ParrinelloRahman`, `MonteCarloBarostat`; default `apply_with_model` propagates `?`; in-module calls currently at barostat.rs:549,567 unwrap.
Runner builders at runner.rs:1508,1520 may still construct experimental objects, with production validation unchanged. Search all final-base
callers/implementers at adoption. Downstream implementers add Result/Ok; callers handle errors. This is documented in release notes.

Add trait `apply_with_pressure(&mut self,system:&mut System,pressure:&PressureObservation,dt:f64)->Result<()>` default Config("tensor coupling
unsupported by this barostat"). Only Berendsen/legacy piston implement experimental unconstrained complete endpoint-diagonal updates;
require complete finite total tensor, same positive volume, no constraints, Diagnostic permitted only for this direct explicitly experimental
API (runner remains closed). All other cases fail before mutation. This API's kernel tests validate tensor use, not canonical NPT sampling.
MC remains energy-based and unconstrained; no molecular Jacobian change. Post-move full ModelEvaluation cache refresh is mandatory.

Restart: keep `Checkpoint`, `ContinuationState`, `save`, `load` source-compatible and exact v2 14-field continuation shape. Introduce separate
public `StressCheckpoint { pub mechanical:Checkpoint, pub stress:Option<RestartStress> }`, `save_with_stress<P:AsRef<Path>>(path:P,cp:&StressCheckpoint)
->Result<()>`, `load_with_stress<P:AsRef<Path>>(path:P)->Result<StressCheckpoint>`. New writer uses MODENA_RESTART 3 with mechanical v2 payload
nested under `mechanical` and stress under `stress`; `save` still writes v2, `load` accepts v1/v2 and rejects v3 with an explicit use-load_with_stress
message. Runner migrates to new functions. New reader accepts v1/v2 via existing parser and gives stress=None. Old exact v2 reader fixture stays.

`RestartStress { estimator_version:u32, solver:ConstraintSolverOptions, model:ModelEvaluation, endpoint:PressureObservation, interval:Option<IntervalPressure>,
last_stages:Vec<ProjectionReport>, last_bath:Option<BathImpulse>, step:usize, dt:f64, physical_signature:serde_json::Value }`.
Estimator version=1, output timing endpoint_after_controllers / interval_pre_controllers included in the physical signature. v3 validates step,
volume, dt, scalar traces, mechanical vector lengths, provenance and completeness before restore. No partial-step save. Full exact continuation
keeps all cached forces and RNG/controller states unchanged. Legacy v2 input has unavailable initial interval; reconstruct endpoint model only
on scratch System and cloned/pinned evaluation context if evaluation is pure and cannot mutate worker counters/cache state; otherwise mark
Missing(LegacyRestart) until first normal force call. Do not evaluate a live external session as a supposedly pure reconstruction. The first
interval missing S0 remains Missing; next complete step restores it. Controllers requiring pressure must reject missing state before advancing
or require explicit state-only restart; current unsupported controllers remain rejected. v3 reader rejects unsupported estimator versions.

## 9. Frozen validation protocol and budgets

All values below are fixed before numerical outcomes; changing them creates a new protocol hash and invalidates prior acceptance.
Deterministic f64 tests use reduced units, one Rayon worker then four for pair/cache parity. No randomized fixture generation without listed seed.
Six-component tolerance uses `abs(error)<=tol*(1+abs(expected))` componentwise, not norm-only. FD requires at least two adjacent eps values
passing and their derivatives differing within twice the stated tolerance; archive all four values, not just the best.

| Gate | Fixed protocol and pass requirement |
|---|---|
| Kernel six strain | eps=1e-4,3e-5,1e-5,3e-6; pair/bonded/direct Ewald tol2e-7; mesh oracle tol2e-6; trace5e-12*(1+sum abs terms). |
| Mesh H0 | U/W1e-11*(1+abs), each force coordinate energy-FD2e-7; integer-cell shifts1e-10; translation error decreases under mesh [8,8,8]→[16,16,16]→[32,32,32] against independently converged Ewald (finest<=2e-4 absolute/component, charged/scaled fixtures). These larger meshes use FFT, not direct DFT. |
| Impulses | projection tol1e-13; momentum/torque/per-atom reconstruction/ΔG error<=5e-10*(1+sum abs ledger terms); residuals reported separately. |
| Endpoint | rotor and principal/general water force and six tensor components<=2e-9*(1+abs expected); linear residual<=1e-11. |
| State/time convergence | hω=.08,.04,.02,.01, Tω=8, same initial conditions and time; last two non-roundoff state-error ratios3–5 for second-order claim. Exact free-rotor finite-step identity mandatory regardless of fitted order. |
| OU conditional | 65,536 independent draws, seeds0..7 split equally, gamma*h=.2,kT=1; means <=4SE plus1e-3; each normalized projected covariance element error<=.02 and<=4SE; raw-reaction mean zero, total covariance matches formula §4. |
| Free rigid gas | 32 rotors or32 waters, cube edge12, kT=1,gamma=1; no COM removal. Seeds1001..1008; h=.04,.02,.01,.005; burn time50, measurement time400 initially, samples each.1. If precision target fails, extend same frozen seeds to time800 then1600 only; do not change tolerances or seed set. |
| Gas distribution | 2K_com/kT~chi²(3M); rotor 2K_rot/kT~chi²(2M); water~chi²(3M). Per independent trajectory mean/variance and standardized CDF at chi-square .1,.5,.9 quantiles; extrapolated mean/variance within3 combinedSE and1%/3%; CDF error<=.02 and3SE. Orientation E[e_a e_b]=delta_ab/3 within.02 and3SE. Test both endpoint and interval means independently. |
| Gas pressure | analytic p=32/1728=1/54; total tensor mean diag1/54 offdiag0. Target95% CI halfwidth<=.005*p for diagonal/scalar and<=.005*p for shear; extrapolated error<=.01*p and<=3SE. Report finite-dt h rows; at finest h same practical1% target and3SE. |

Block protocol: sample interval .1; discard burn-in; estimate integrated autocorrelation by initial-positive-pair sequence for pressure, COM energy,
rotational energy, orientation. Choose common block length ceil(10*max tau_int) samples per fixture/dt, using the first half only; use the
second half for inference. Require at least20 blocks/seed and N_eff>=1000 pooled per observable. Primary uncertainty is between the eight
independent seed means (Student-t, df7), backed by block bootstrap with fixed resampling seed7001, 2000 resamples. Insufficient precision at
max duration is Inconclusive, not Passed. Weighted quadratic-in-h extrapolation uses all four steps and reports linear-in-h² lack-of-fit;
compare dropping the coarsest value, intercept shift<=combined2SE and.005*p. A second-order weak-error claim requires bias resolved above2SE
at two steps and ratios consistent with order2 (95% interval inside[1,3]); unresolved slope means no order claim, while distribution/pressure
agreement may pass if their independent practical and statistical targets pass. Do not impose monotonicity on noisy unresolved estimates.

Analytic fixtures: rotor masses(1,3) and(2,5), d1.3,omega.7, e=(1,2,2)/3,t=(2,-1,0)/sqrt5, COM V=(.2,-.1,.3), R=(5,6,7), cube20.
Test each face/corner split and a common translation(.31,-.27,.19). First rotor exact endpoint trace M|V|²=.56, internal trace cancellation;
mu*d²*omega²=.621075. Negative-control interval/endpoint mixing error `-.621075*(1-c)/(1+c)`, c=sqrt(1-(.7h)²), must be observed.
Water masses15.9994,1.008,1.008 then16,1,2; OH1,HH=2sin(109.47deg/2), omega.7. Rotate by Rodrigues axis(1,2,3)/sqrt14 angle.41;
principal axis is rotated plane normal; nonprincipal Ω=(.2,.3,.4). Use analytic Euler acceleration, not production matrix inverse, as oracle.
Forced rotor uses external harmonic trap U=.05*Σ|q_i-(5,6,7)|², exact F=-.1*(q_i-center), plus an oblique orientation; confirm b-force term.

Pair fixture harmonic U=r²/2 with d=(1,2,3), expected S=[-1,-4,-9,-2,-3,-6], W=-14 and U=7. Bond/angle/dihedral fixtures are oblique,
nonplanar where required, and include an angle whose trace is zero but shear nonzero. Electrostatics cube edges(8,9,10), fractional positions
(.13,.27,.39),(.61,.44,.73),(.32,.81,.58), charges(1,-.7,.2), κ=.8,rc3.5, fixed direct modes radius8; scale(0,1)=.4 and exclusion(1,2)=0.
One additional excluded pair beyond rc uses positions(.1,.1,.1),(.59,.4,.3) and fixed image branch. The original second x=.6 is retained only as a rejected exact-half-cell tie negative control; the replacement was fixed by pure geometry proof before its FD energies. No cutoff seams in FD. Ewald reciprocal
modes fixed under strain; scalar background standalone U_bg=-πQ²/(2Vκ²) has tensor U_bg I, independently tested with net charge.
Supplementary triclinic direct-Ewald cell columns are (8,0,0),(.8,9,0),(.3,.1,10). This was the first all-branch/all-48-membership-valid cell in a predeclared eight-candidate geometry-only roster; earlier crossing cells and protocol revisions are retained. The choice preceded accepted derivatives and changes no kappa, rc, kcut, charges, modes, eps or tolerance.

Work caps: deterministic matrix<=60s CPU/test binary and<=512MiB RSS; total standard suite<=10min one worker. Scientific campaign<=2h wall,
<=4 CPU workers (sequential jobs), <=512MiB Rust RSS / <=1GiB LAMMPS RSS, <=5GiB raw artifacts. ≤64 interacting waters, ≤32 ideal molecules,
≤8 seeds, frozen dt/length ladders above; ≤400million total molecule-steps for gas ensemble (counts all dt,seeds and initial+extension duration;
continuation extends rather than restarts). Before execution calculate exact upper work and reject any campaign configuration over the cap.
Conditional OU run counts against a separate cap1million projected draws. External runs capped24 individual jobs, total1h within campaign2h.
Measure wall/RSS with host time wrapper, retain actual values and physical-step counters, not static counters. Stop at budget with retained
partial rows. No extended run while B03/B05 holds timing resources. These are budget limits, not speedup targets or performance claims.

## 10. External matched validation — real execution required

Create `bench/b07_pressure_reference.py`, `bench/fixtures/b07/`, example `stress_validation.rs`, modules `stress_validation/{mod,fixtures,
mesh_oracle,statistics,protocol}.rs`, and `docs/scientific-hardening/B07-PRESSURE-VALIDATION.md`. Use pinned LAMMPS CPU double with RIGID+KSPACE.
Preparation resolves one installed or authorized official release into an immutable source revision and binary SHA256, captures `-h`, package
list, compile options, units and exact executable path before numerical evaluations. No guessed command feature: fail preflight if fix rattle,
compute pressure contribution switches or required output are unsupported. An unavailable executable is a pending required gate; no synthetic
success. Freeze the exact version and all generated input/snapshot hashes in protocol.json before the first physics run. Reference preparation
may install only with existing user authorization; this draft does not perform installations.

Use `units lj`, `atom_style full`, explicit masses/charges/bonds/angles, periodic orthorhombic cells, no removal of free COM. Reference pressure
Voigt order(xx,yy,zz,xy,xz,yz), compressive sign, kB=1. For rotor and free water use zero nonbonded force, explicit constraints with fix rattle
1e-13 100, `fix nve`, pressure split into kinetic/configurational/fix terms. Match snapshots, masses, distances, dt, initial velocities,
no thermostat and no COM removal. Export endpoint q,v and conservative F at every retained sample. Preserve actual LAMMPS kinetic timing
and its impulse virial convention in metadata. Check same-snapshot conservative tensors plus analytic rotor/water endpoint oracle, and compare
reference impulse pressure to MODENA interval and endpoint only via common-time-window dt extrapolation; no assumed same finite-step scheme.
For nonprincipal rigid body compare rotational state convergence and Euler-acceleration endpoint force oracle at reference snapshots.

Interacting fixture:32 waters on 4x4x2 lattice, cube12, listed water geometry/masses, all O LJ epsilon=.1 sigma1 cutoff2.5 shifted-force,
H LJ zero, charges(-.8,.4,.4), intrawater Coulomb/LJ excluded. Initialize orientations by quaternion from seed901 and COM/rotation Maxwell
draws seed902; freeze full coordinate/velocity arrays, not generator-only provenance. kT1 for NVT comparison, fixed κ=.8,real cutoff3.5;
MODENA converged direct Ewald is the reference branch first, then SPME. Native LAMMPS LJ/coulomb cutoff conventions must match the exact
MODENA Hamiltonian: implement the small explicit shifted-force LJ table in the reference input with independently generated values/derivatives,
retain it, and verify pair table U/F to2e-7 before trajectories. Use reference Ewald, not PPPM force-error parameter as a tensor certificate.
Converge reciprocal reference across tolerances1e-7,1e-9,1e-11 on frozen snapshots; require final tensor changes<=2e-6*(1+abs), F<=2e-6.
Freeze actual reciprocal settings after this predefined convergence phase, never fit to MODENA outcomes. Verify real-space/exclusion branches
and energy shifts on pair snapshots before dynamics; inability to match them blocks the interacting reference gate.

External trajectory roster: rotor and water NVE hω=.08,.04,.02,.01,Tω8 (8 jobs); interacting32 waters NVE h=.01,.005,.0025,.00125,T=2
(4 jobs); interacting NVT h=.01,.005, four seeds2001..2004, burn20,measure100 (8 jobs), plus at most4 preflight/convergence reference jobs.
Reference NVT uses LAMMPS `fix rigid/nvt` with molecule IDs, target1, damping1 and explicitly configured chain length3 (`tparam 3 1 3`). The installed pinned-version documentation has default10; explicit configuration preserves this protocol's required3. It replaces fix nve/rattle for these NVT jobs; do not apply both integrators. Its ensemble comparison is independent of MODENA BAOAB, so compare distributions and dt limits rather than random sequences or finite-step trajectories. Two of the four preflight jobs validate this reference on32 free rotors and32 free waters at h=.005,burn50,measure1000,seed2001, no COM removal, against χ² kinetic and ideal pressure targets using frozen block analysis/practical tolerances. The other two preflight jobs perform all frozen snapshot reciprocal/table checks, reusing the process for the finite tolerance ladder. If this pinned reference cannot satisfy its free-gas gates or support its declared thermostat settings, fail the dependent NVT gate instead of swapping schemes. Independent Python dense projector/RATTLE code may cross-check splitting but is not required and does not replace the external gate.

Matched gates: identical-snapshot U/F/conservative tensor<=2e-6 scaled; analytic free rotor/water endpoints as §9; common-window continuum
mean all six total pressures difference<=3 combinedSE and1% of max(kT*M/V,abs(reference diagonal)); shear uses ideal pressure scale, not
relative-to-zero. Deterministic NVE extrapolated average tolerance2e-5*(1+norm pressure), finest state ratios3–5. Interacting ensemble CI
halfwidth<=1% ideal scale and N_eff>=500; if the frozen roster cannot reach precision, mark Inconclusive and report the extra measured cost
needed, without calling B07 accepted or widening errors. Archive input/stdout/stderr/binary hash/raw snapshots and both estimators per row.

## 11. Acceptance and final report

Default pressure validation remains Diagnostic. Scientific report has a row for each provider/backend/integrator/constraint/controller tuple,
with Passed, Failed, Inconclusive, Unsupported or Pending. Validated records identify the immutable report hash and validated domain; serialized
self-asserted status never overrides local compatibility. Scalar-only, backend-incomplete, unsupported periodic winding, or experimental controller
cases cannot be promoted. NVE unconstrained tensors require six strain/trace gates; constrained endpoint and interval each require analytic,
impulse, dt and external matches; BAOAB adds OU/distribution/uncertainty gates. Completing code is not acceptance. All required supported cases
must finish with actual outcomes; a blocking numerical/reference problem is reported explicitly, not hidden behind a fallback.

Costs must report: six-component pair/bond/reciprocal accumulation and workspace deltas, stage ledger entries/bytes, rollback scratch, dense
reaction block sizes/solve time, output/checkpoint sizes, optional profile instrumentation, independent DFT counts, campaign wall/RSS/work,
and barostat source-break migration. There is no speedup requirement. Compare old/new whole-step cost on frozen tiny and32-water fixtures
only after exclusive timing lease release; do not compete with B03/B05 or fabricate counter values. Retain all failed physical rows.

## 12. Validation artifact records

Example-only records use serde JSON with deny_unknown_fields; they are not a production input format. Units are reduced throughout.
All artifact records below are binding; protocol/provenance schemas are shared by the example and launcher:

```rust
pub enum CaseStatus { Prepared,Passed,Failed,Inconclusive,Unsupported,Pending,BudgetExceeded }
pub struct GateResult { pub name:String,pub observed:Vec<f64>,pub expected:Vec<f64>,
    pub absolute_error:Vec<f64>,pub tolerance:Vec<f64>,pub status:CaseStatus }
pub struct Sample { pub seed:u64,pub step:u64,pub time:f64,pub dt:f64,
    pub endpoint:[f64;6],pub interval:[f64;6],pub com_energy:f64,
    pub rotation_energy:f64,pub orientation_second:[f64;6] }
pub struct Estimate { pub mean:f64,pub se:f64,pub ci95:[f64;2],pub n_eff:f64,
    pub block_samples:usize,pub retained_blocks:usize }
pub struct ValidationFixture {
    pub id:String, pub positions:Vec<[f64;3]>, pub velocities:Vec<[f64;3]>,
    pub masses:Vec<f64>, pub charges:Vec<f64>, pub cell:[[f64;3];3],
    pub constraints:Constraints, pub kappa:Option<f64>, pub real_cutoff:Option<f64>,
    pub reciprocal_cutoff:Option<f64>, pub pair_scales:Vec<(usize,usize,f64)>,
}
pub struct WorkLimits {
    pub wall_seconds:u64, pub max_workers:usize, pub rust_rss_bytes:u64,
    pub reference_rss_bytes:u64, pub artifact_bytes:u64, pub molecule_steps:u64,
    pub projected_draws:u64, pub reference_jobs:usize,
}
pub struct ValidationProtocol {
    pub version:u32, pub solver:ConstraintSolverOptions,
    pub strain_eps:Vec<f64>, pub gas_seeds:Vec<u64>, pub gas_dt:Vec<f64>,
    pub burn_time:f64, pub measure_times:Vec<f64>, pub sample_interval:f64,
    pub fixtures:Vec<ValidationFixture>, pub limits:WorkLimits,
    pub specification_sha256:String, pub source_sha256:std::collections::BTreeMap<String,String>,
    pub snapshot_sha256:std::collections::BTreeMap<String,String>,
}
pub struct ExecutionProvenance {
    pub source_head:String, pub compiler:String, pub build_profile:String,
    pub features:Vec<String>, pub command:Vec<String>, pub protocol_sha256:String,
    pub binary_sha256:String, pub external_source_revision:Option<String>,
    pub external_packages:Vec<String>, pub threads:usize,
    pub started_utc:String, pub ended_utc:String, pub elapsed_seconds:f64,
    pub peak_rss_bytes:u64, pub molecule_steps:u64, pub projected_draws:u64,
    pub reference_jobs:usize,
}
pub struct ValidationReport {
    pub protocol:ValidationProtocol, pub provenance:ExecutionProvenance,
    pub gates:Vec<GateResult>, pub raw_files:std::collections::BTreeMap<String,String>,
    pub status:CaseStatus,
}
```

Validation rejects protocol-version!=1, nonfinite arrays/parameters, inconsistent atom-vector lengths, duplicate fixture IDs,
invalid boxes/masses/constraints, unsupported seeds/dts/limits or spec hash mismatch. Exact native constraints preserve their serde schema.
Statistical tolerances are normative spec §9 constants bound by specification_sha256, not undocumented mutable CLI arguments.
Record actual elapsed work counts as each phase executes; continued phases load consumed budgets. No counter may be initialized to a claimed
measured total. Output directories are new on prepare, resumed only with matching immutable protocol/source hashes. Status reflects worst required
gate (Failed/BudgetExceeded/Inconclusive/Pending before Passed). Raw file map is path→SHA256. External identity is frozen before physical runs.

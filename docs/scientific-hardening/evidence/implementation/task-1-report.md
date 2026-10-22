# Task 1 report: sampling, constraints, DOF, continuation state and replicas

Status: complete within Task 1 scope, with the explicitly permitted constraint-pressure limitation. Scientific support is bounded by the tests below; this is not a proof of arbitrary constrained molecular canonical sampling.

Implementation commit: `c6fef76ae58500bd545930224cda0980c607fed0` — `:bug: Fix: canonical sampling constraints and replica state`.

Worktree: `/Users/ricardoleal20/.codex/worktrees/modena-scientific-hardening/modena`; branch `ricardoleal20/modena-scientific-hardening`; Rust `1.97.1 (8bab26f4f 2026-07-14)`. Original audit notes/probes were read at `/Users/ricardoleal20/Documents/Code/Misc/modena/docs/audits/2026-10-03-modena-vs-lammps` and preserved. No runner/config production file was edited. No external publication or subagent delegation occurred.

## Delivered behavior

- CSVR samples `chi_square(f-1) = 2*Gamma((f-1)/2,1)`, uses the canonical signed factor, and evaluates the transition as a sum of squares to avoid cancellation. `f=1` has a zero chi-square contribution and never invokes a shape-zero Gamma sampler. The previous implementation-derived unit expectation was replaced with canonical moments.
- Constrained BAOAB projects each drift independently, applies only that drift's position-correction impulse to velocity, projects OU noise onto the mass-weighted tangent space, and projects after kicks. It never replaces the OU velocity with the full step displacement. Velocity-Verlet and NHC use the same correction-only drift helper, avoiding the former displacement-based periodic aliasing. The free atom absent from a constrained dimer retains variance `T/m`.
- SHAKE/RATTLE/SETTLE-class Newton projections return `Result<()>`. Nonfinite or invalid geometry parameters, invalid indices, nonpositive masses, mismatched buffers/references, overlapping water triples, invalid tolerance, singular solves, and residuals above tolerance are errors. Residual diagnostics include all water distances and tangent components. Selected water triples may coexist with unrelated atoms; the former assumption that all position chunks correspond to waters is gone. All production callers propagate errors and existing successful tests explicitly unwrap them.
- System carries explicit independent constraint count and conserved COM DOF metadata, defaulting to zero reductions for manually constructed systems. Integrators record constraint count before temperature/thermostat operations. Langevin records zero conserved COM reduction because its COM modes are stochastic. Resizing resets the metadata; cloning preserves it.
- Triclinic minimum imaging enumerates the lattice translations that could improve the rounded fractional candidate. If its Cartesian norm is `R`, every closer image satisfies `|s_k-n_k| <= ||row_k(h^-1)|| R`, which gives finite bounds containing all candidate improvements. This covers translations outside a fixed adjacent-cell stencil. Uniform and Cartesian-axis box scaling rebuild basis/inverse/bounding size consistently.
- Integrator and Thermostat expose checkpoint/restore hooks. CSVR and Langevin use explicitly named ChaCha12 via the existing rand-0.10-compatible `chacha20` crate, rather than the replaceable StdRng implementation. The snapshot contains the key, stream ID and exact 68-bit word offset (encoded as decimal text), with an algorithm/version marker. NHC snapshots include every chain momentum and position. JSON uses `float_roundtrip`. Restoring validates component kind/version, parameters, seed, and constraints for coupled integrators; invalid restoration is checked before mutating live state.
- Replica exchange constructs validated positive-temperature, positive-dt, nonnegative-gamma replicas with full Maxwell velocities including COM. `run` obtains its real-space cutoff from the model, keeps PE returned by every integrator step, exchanges entire System and matching neighbor caches, and rescales incoming velocities by the fixed destination/source temperature ratio. The first two-replica interval now attempts its pair. Acceptance is expressed as `log(weight_after/weight_before) = (beta_i-beta_j)*(Ui-Uj)`. The old executable acceptance sign was already equivalent; the incorrect module equation was corrected, not blindly inverted.

## Public interfaces / integration handoff

```rust
// System: setters validate reductions against 3N; COM count must fit 0..=3.
fn constraint_count(&self) -> usize;
fn set_constraint_count(&mut self, count: usize) -> Result<()>;
fn conserved_com_dof(&self) -> usize;
fn set_conserved_com_dof(&mut self, count: usize) -> Result<()>;
// dof() now returns 3N - constraint_count - conserved_com_dof.

// Constraints: count assumes independent/nonredundant distances.
fn count(&self) -> usize; // bonds.len() + 3*waters.len()
fn validate(&self, system: &System) -> Result<()>;
fn shake(&self, system: &mut System, x_ref: &[f64], tol: f64, max_iter: usize) -> Result<()>;
fn rattle(&self, system: &mut System, tol: f64, max_iter: usize) -> Result<()>;
fn settle(&self, system: &mut System, tol: f64) -> Result<()>;
fn settle_velocities(&self, system: &mut System, tol: f64) -> Result<()>;
fn max_velocity_residual(&self, system: &System) -> f64;

// Integrator and Thermostat: stateless defaults encode null and reject non-null.
fn checkpoint_state(&self) -> crate::error::Result<serde_json::Value>;
fn restore_state(&mut self, state: &serde_json::Value) -> crate::error::Result<()>;

// ForceModel default is None; concrete providers are Task 2 ownership.
fn neighbor_cutoff(&self) -> Option<f64>;

// ReplicaExchange::new now returns Result<Self> with its existing arguments.
// run's dt must equal the validated construction dt.
```

Task 2 must implement cutoff providers; until then direct replica use of a provider retaining the trait default fails clearly. `tests/ensembles.rs` temporarily supplies the known LJ cutoff through a test-only delegating adapter so the Task 1 commit builds and tests independently. Task 2 can remove that adapter after PairForceModel implements the seam. Triclinic neighbor stencils still need Task 2's safe enumeration/fallback: the repaired min-image method alone cannot prove Cartesian bounding-box cell-list correctness. Audited fixture: `a=(10,0,0), b=(4.9,10,0), c=(0,0,10)`, displacement `(6.37,3,0)`, correct image `(-3.63,3,0)` with distance `4.709235606762524 < 4.8`. A second test uses tilt 49 and a nearest translation of `2*a`, ruling out a fixed +/-1 fractional search.

Task 3 must set reductions before velocity initialization and output, project initial positions/velocities, and set conserved COM to three only for initialized zero-momentum dynamics that preserve it. Integrator recording during a step is not a substitute for this initialization wiring. Task 3 must serialize System metadata and physical identity with the mechanical state and use these component state hooks for full continuation. Timestep belongs to the runner's physical signature because the direct Integrator trait receives it per step.

## Observed RED evidence

All commands ran in the worktree and began with `rtk`.

1. Before production sampling/constraint/geometry fixes:

```text
rtk proxy cargo test -p modena-core --test sampling_hardening -- --nocapture
exit 101: 0 passed, 7 failed
CSVR conditional mean K': 3.609192169134449, canonical 6.045571576502769
CSVR stationary mean K (seed 7): 3.25273591836196, canonical 6
Signed decorrelated CSVR: 0 negative draws out of 10000; expected approximately half
BAOAB free atom with unrelated dimer, mass 2, seed 7 axis 0:
  variance 0.1980496380752639, independent T/m expectation 0.5
Singular SHAKE: integrator silently returned success
Skew image: (6.370000000000001,3.0000000000000004,0), expected (-3.63,3,0)
Triclinic uniform scale 2: volume remained 1000, expected 8000
```

2. Replica tests before the fixes:

```text
rtk proxy cargo test -p modena-core replica_exchange::hardening_tests --lib -- --nocapture
exit 101: 1 passed, 5 failed
Accepted swap cached force: -2 instead of incoming configuration's -8
Final PE without swaps: retained 2.5 vs final-geometry energy 202.4145231817191
First two-replica swap interval: 0 attempts instead of 1
Unknown model cutoff: silently guessed and returned success
Nonfinite dt: silently integrated and returned success
Independent strongly unfavorable Boltzmann test already passed.
```

3. New state and DOF tests were first written against missing interfaces (compile RED). After adding the required null defaults and metadata API, the state behavior produced a real failing regression before its stateful implementations:

```text
rtk proxy cargo test -p modena-core --test sampling_state -- --nocapture
exit 101: 3 passed, 5 failed
CSVR split continuation velocity: 0.7738631017771412 vs uninterrupted 0.898211326757774
Langevin split continuation: differing velocity arrays
NHC split continuation velocity: 0.7520106853551479 vs uninterrupted 0.8123053029625782
Different CSVR parameters wrongly accepted null state
Constrained Langevin DOF: 9 instead of required 11 (four atoms, one distance, stochastic COM)
```

4. Additional observed checkpoint RED before its defensive fix:

```text
rtk proxy cargo test -p modena-core --test sampling_state checkpoint_rejects_nonfinite -- --nocapture
exit 101: 0 passed, 1 failed
NaN constraint length was silently serialized instead of returning an error.
```

The Result-oriented constraint error cases were written before the signature change; the existing singular-integrator regression above provided the behavioral RED for the former silent projection failure. The final cases cover iteration exhaustion, invalid mass/index input, singular/invalid water, and selected water/free-atom behavior.

## Independent expectations and statistical tolerances

- Exact CSVR transition with `c=exp(-dt/tau)`: conditional mean `c*K+(1-c)*f*T/2`; conditional variance `(1-c)^2*f*T^2/2 + 2*c*(1-c)*K*T`. 80000 draws, absolute tolerances 0.045 on mean and 0.14 on variance.
- Canonical stationary kinetic energy is `Gamma(f/2,T)`, independent of the implementation. At `f=12,T=1`, mean and variance are both 6. Three fixed seeds `[7,42,91]`, 1000 burn-in and 40000 retained samples each; mean tolerance 0.085, variance tolerance 0.25. At `f=1`, mean and variance are 0.5, with 40000 retained samples and tolerances 0.015/0.04.
- An unrelated unconstrained mass-2 atom at `T=1` must have every Cartesian velocity variance 0.5. Three seeds, constrained and unconstrained variants, 1000 burn-in and 25000 retained samples each; per-component tolerance 0.035. This tests the audited defect; it does not certify water's full orientation/rotation distribution.
- CSVR signed decorrelation uses 10000 independent resets with 4500..5500 negative draws accepted. Restart tests compare exact f64 arrays after 13 uninterrupted steps, a JSON serialization/deserialization seam, and 17 further steps; all three stochastic/extended state paths agree bit-for-bit.

## GREEN verification

```text
rtk proxy cargo test -p modena-core --test sampling_state --test sampling_hardening --test constraints_hardening
exit 0: 10 state + 9 sampling/geometry/pressure assessment + 5 constraint tests passed.

rtk proxy cargo test -p modena-core
exit 0: 166 passed, 0 failed, 1 existing doc test ignored.
Harness counts: lib80; bonded6; config_builder3; conservation3; constraints5;
decomp5; electrostatics9; ensembles3; gpu_pair0(default feature off);
io_and_lattice17; mpi_decomp0(default feature off); property_consistency9;
restart3; sampling9; state10; water4; doc0 passed/1 ignored.

rtk proxy cargo clippy --all-targets -- -D warnings
exit 0, workspace default targets checked.

rtk proxy cargo clippy --all-targets --features gpu -- -D warnings
exit 0, workspace GPU-feature targets checked; actual GPU execution not claimed.

rtk proxy cargo fmt --all -- --check
exit 0.

rtk git diff --check
exit 0.
```

Earlier covering checks independently passed NVE conservation, both constrained bonded conservation cases, NHC extended Hamiltonian, replica integration, rigid-water Newton geometry, full rigid-water/SPME NVE, and existing v1/NVE/protocol-switch restart tests. No scientific performance claim is inferred from verification runtimes. MPI execution was not performed in Task 1.

## Constraint-pressure assessment and remaining limits

An equal-mass centered dimer of length 1 rotating with endpoint speed 0.5 has `K_rot=0.25`. Its physical centripetal constraint virial must be `W_constraint=-2*K_rot=-0.5`; a nontranslating isolated rotor has no ideal translational pressure. The dedicated assessment computes the velocity-change impulse at the midpoint of the independently known rigid motion:

```text
rtk proxy cargo test -p modena-core --test sampling_hardening rotating_rigid -- --nocapture
exit 0
rotor K=0.25
independent impulse W=-0.5000000000023244
reported W=0
reported P=0.00000016666666666662456
complete rotor P=0
```

**Constrained pressure remains diagnostic and incomplete.** This single-rotor identity does not establish the correct impulse stress for overlapping constraints, the Newton water projection, BAOAB's different substeps, thermostat impulses, or box coupling. No generalized scalar or tensor constraint virial was invented. Task 3 must expose this limitation and reject constrained NPT/box rescaling, as allowed by the accepted brief/spec. Fixed-volume constrained dynamics remain supported with the tested geometry/conservation properties.

Further practical limits: constraint counts assume nonredundant distances; a general rank detector was not added. Task 3 validates topology/contradictory duplicates at the input boundary. CSVR needs a nonzero initial kinetic direction and intentionally cannot heat an exactly resting system by global rescaling; zero/empty DOF safely do nothing. Arbitrarily ill-conditioned triclinic bases can make the rigorous finite search expensive; singular/nonfinite cells require input rejection. Replica exchange currently supports unconstrained systems only and fails explicitly for constraint metadata. The driver remains sequential; no new parallel-execution claim is made. Complete runner restart identity/atomic files, concrete cutoff implementations, triclinic neighbor fallback, and constrained-pressure reporting remain Tasks 2/3 responsibilities.


## Fix round 1: rigid-water reference-direction drift projection

Independent review: `.superpowers/sdd/2026-10-03-scientific-hardening/task-1-review.md`.
Immutable FIX_BASE: `c6fef76ae58500bd545930224cda0980c607fed0`.
Fix commit: `dd22f471a40781085c9eb375e0afcc1f7be239a9` — `:bug: Fix: reference-direction rigid-water RATTLE projection`.
Fix-round status: implementation and verification complete; independent scoped re-review pending. All covering tests, including the full water trajectory, formatting, and default workspace Clippy passed.

### Findings addressed and implementation scope

The review correctly identified first-order rotational dissipation from the current-normal coupled water position projection. Adding the position-correction impulse and final tangent projection forms a RATTLE drift only when the position solve holds the **saved pre-drift constraint gradients** fixed. Repeatedly reorienting those gradients at the predicted positions cools a force-free rotor.

Changed files:

- `crates/modena-core/src/constraints.rs`: dynamical SHAKE now passes each selected water's saved reference positions to the coupled Newton solve. Residuals use the predicted bond vectors; Jacobian columns and the mass-weighted correction use the fixed reference bond vectors. This Jacobian is generally nonsymmetric. Existing general 3x3 Cramer solving, final residual/error handling, selected-water indexing, and velocity tangent projection remain applicable. Standalone `settle` remains a geometric initialization projection using current directions; its documentation now distinguishes initialization from dynamics. Removed the obsolete parallel-contiguous-layout claim.
- `crates/modena-core/tests/rigid_sampling.rs`: independently expected free-rotor conservation and nonlinear rigid-body rotational-energy equilibrium tests.

No public API, checkpoint format, runner production file, parameter, random generator, or pressure calculation changed. The already-correct correction-only `project_drift` logic serves BAOAB, velocity Verlet, and NHC; correcting the water SHAKE directions fixes its dynamical impulse. Root's uncommitted plan edits were preserved and are excluded from the fix commit.

### Observed RED before the production fix

After fixing an initial test-only ambiguous float literal, the new regressions produced actual behavioral failures against FIX_BASE before any production change:

```text
rtk proxy cargo test -p modena-core --test rigid_sampling -- --nocapture
exit 101: 0 passed, 2 failed
rotor baoab=true  dt=0.01  initial_K=1.5 final_K=1.428571428501615 max_relative_drift=4.761904766559022e-2
rotor baoab=false dt=0.01  initial_K=1.5 final_K=1.363636363665614 max_relative_drift=9.090909088959043e-2
rotor baoab=true  dt=0.005 initial_K=1.5 final_K=1.463414634099744 max_relative_drift=2.4390243933503974e-2
rotor baoab=false dt=0.005 initial_K=1.5 final_K=1.428571428501615 max_relative_drift=4.761904766559022e-2
rotor baoab=true  dt=0.001 initial_K=1.5 final_K=1.492537314974684 max_relative_drift=4.975123350210886e-3
rotor baoab=false dt=0.001 initial_K=1.5 final_K=1.485148514448276 max_relative_drift=9.900990367815887e-3
rigid rotation gamma=0.1 dt=0.04 samples=300000
  mean=1.2658767246023268 variance=0.9590147535354732
  canonical_mean=1.5 canonical_variance=1.5
  mean_tolerance=0.18677716670002661 variance_tolerance=0.5403315001000798
```

Both tests failed. The conservation assertion collected all six rotor cases before failing. The rotational test failed the weak-friction mean assertion before reaching its second parameter case.

### Independent physical expectations and tolerances

Three unit masses form an equilateral triangle of circumradius 1 and constrained distance sqrt(3), centered at (50,50,50) in a cubic box of side 100. Initial angular velocity is 1 about its z principal axis. Zero physical force and zero Langevin friction imply constant kinetic energy `K=3/2` at every time. Both BAOAB and velocity Verlet integrate to `t=10` at `dt=[0.01,0.005,0.001]`; the maximum relative energy drift over every step must remain below `1e-7`, and final squared-distance residuals must be below `1e-10`. These thresholds reject the measured first-order cooling by orders of magnitude.

For equilibrium, translational kinetic energy is removed by subtracting the instantaneous mass-weighted COM velocity. A nonlinear rigid body has three rotational quadratic modes; canonical tangent momenta therefore give `K_rot ~ Gamma(3/2,T)`, with mean and variance both `1.5` at `T=1`. Cases `(gamma,dt)=(0.1,0.04)` and `(2,0.01)` each use three seeds `[7,42,91]`, twenty friction times of burn-in, and 100000 retained samples per seed (300000 pooled). Statistical tolerances use five standard errors with the OU energy correlation estimate `rho=exp(-2*gamma*dt)` and inflation `(1+rho)/(1-rho)`, plus a `0.01` finite-timestep margin. The canonical variance of the centered squared-energy estimator is `13.5`. This accounts explicitly for the weak-friction samples' correlation; it is not an independence assumption or a proof of arbitrary molecular canonical sampling.

### GREEN observed with the reference-direction solve

```text
rtk proxy cargo test -p modena-core --test rigid_sampling -- --nocapture
exit 0: 2 passed, 0 failed
rotor baoab=true  dt=0.01  final_K=1.500000000009644 max_relative_drift=2.1194823673909013e-11
rotor baoab=false dt=0.01  final_K=1.499999999985798 max_relative_drift=1.2918222047630934e-11
rotor baoab=true  dt=0.005 final_K=1.499999999814251 max_relative_drift=1.4215795207661586e-10
rotor baoab=false dt=0.005 final_K=1.500000000009644 max_relative_drift=2.1194823673909013e-11
rotor baoab=true  dt=0.001 final_K=1.499999997961690 max_relative_drift=1.358873458201515e-9
rotor baoab=false dt=0.001 final_K=1.500000000456273 max_relative_drift=3.905971102113881e-10
rigid rotation gamma=0.1 dt=0.04 samples=300000
  mean=1.5256781848314105 variance=1.5663987421237122
  mean_tolerance=0.18677716670002661 variance_tolerance=0.5403315001000798
rigid rotation gamma=2 dt=0.01 samples=300000
  mean=1.5056806796694442 variance=1.482278587162741
  mean_tolerance=0.08906221165077505 variance_tolerance=0.2471866349523252

rtk proxy cargo test -p modena-core constraints::tests --lib
exit 0: 4 passed, 0 failed (76 unrelated filtered)

rtk proxy cargo fmt --all -- --check
exit 0, no output

rtk proxy cargo clippy --workspace --all-targets -- -D warnings
exit 0: core, MCP, CLI, bench default targets checked
```

The retained spectator and exact direct JSON restart seam regressions pass in the covering run. No unsupported-water guard was added: fixed-volume rigid-water dynamics remain supported. Constrained pressure is still diagnostic/incomplete because constraint impulse virial is absent; constrained NPT/box-rescaling rejection and pressure metadata remain required Task 3 work. The equilateral free-rotor and energy statistics test specific independent physical requirements, not all orientation distributions or arbitrary constrained configurations. Existing singular/unconverged projection errors remain necessary for excessive timesteps or bad reference geometry.


Final covering verification after formatting:

```text
rtk proxy cargo test -p modena-core --test rigid_sampling --test sampling_hardening --test sampling_state --test constraints_hardening --test bonded --test conservation --test water
exit 0: 39 passed, 0 failed
Harness counts: bonded6; conservation3; constraints5; rigid_sampling2;
sampling9; state10; water4.
Rigid-water/SPME stack NVE trajectory: passed (water harness 66.35 s).
Together with constraints::tests --lib: 43 covering tests passed.
```

The unchanged whole-core/GPU-feature checks from the initial implementation were not rerun in this focused fix round. Default workspace Clippy covers every default target; no actual GPU or MPI execution is claimed. No timing/performance claim is made from the test runtimes.

`rtk git diff --check` exited 0 before the fix commit. After commit, `rtk git status --short` shows only the root-owned plan edit; both owned source/test files are clean.

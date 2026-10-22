# Task 3 report: validated runner, isotropic MC NPT and v2 continuation

Status: implemented and covering checks pass. Independent review remains with the orchestrator.

Implementation commit: `c5810faca4324b5ab3b6de8eeddff028c8122eb8` — `:bug: Fix: validate ensembles and preserve complete checkpoints`. Immutable task BASE: `9fa3b9e6fefa326a16a473ccade05709d7e76585`.
Worktree: `/Users/ricardoleal20/.codex/worktrees/modena-scientific-hardening/modena`; branch `ricardoleal20/modena-scientific-hardening`. No subagents, remote publication, or edits to the original checkout/audit evidence. Read the exact brief, full spec/plan, runner handoff and Task 1/2 public interfaces. Applied test-driven-development and verification-before-completion.

## Delivered behavior

- All configuration entry points go through `InputConfig::normalized()`: raw real-unit TOML/JSON/Rust configs convert once; populated `unit_factors` mark already-normalized real configs. Validation runs before and after conversion. Real neighbor skin also converts as a length. Signatures use normalized parameters and unit factors. Output thermodynamics remain kcal/mol, K, atm for real configs (ps and Angstrom for time/coordinates); reduced configs report engine units.
- Strict serde records/enums reject unknown keys, including nested keys. The common boundary checks finite fields, positive dt/tau/mass/box/lengths, nonnegative T/gamma/skin, nonzero cadences, structure composition, generator arithmetic, canonical coupling compatibility and electrostatic parameters. Live-system checks validate atom masses, geometry/kinetic finiteness and topology/constraint indices before any force evaluation or output opening. Normalized duplicate constraints are rejected to avoid counting the same constraint twice. Contradictory pair overrides and bonded pairs are rejected. Valid zero LJ parameters for noninteracting water H types remain accepted by their existing potential constructors.
- Anisotropic scalar-virial coupling, constrained box rescaling, and the legacy piston return explanatory production-config errors. MC NPT requires a positive-temperature canonical Langevin (gamma > 0), NHC, or CSVR controller with matching temperature. Berendsen remains equilibration only. Positive-target CSVR cannot initialize or restore zero kinetic energy silently.
- Initial constraint positions and tangent velocities are projected before final temperature scaling. The runner sets constraint count and COM reductions before initialization: VV/NHC/global rescaling preserve an initialized zero COM and normally subtract three DOF; Langevin draws and retains stochastic COM modes. State-only protocols and later full continuations preserve COM modes already present in their loaded state. Manual System defaults are unchanged.
- The model-aware barostat hook refreshes forces, PE, virial and neighbors after box changes before observables or the next kick. MC accepted moves retain freshly computed caches; rejected/failed moves restore the full System, scalar virial and PE and rebuild neighbors for the restored configuration. Config construction uses the actual ForceModel cutoff and the fallible `NeighborList::with_pair_scales`; existing binary exclusions are translated to 0/0 scales. Graph/fractional topology wiring remains Task 4.
- Isotropic MC uses a uniform symmetric proposal in log V and `log_accept = -beta*(delta U + P*delta V) + (N+1)*log(Vnew/Vold)`. It preserves velocities. Seed, interval phase, proposal width, acceptance/attempt counters and the exact ChaCha12 state survive checkpoints. Model computation errors propagate after rollback.
- Version 2 writes full mechanics (including forces, PE and virial), atom identity, physical DOF, normalized physical signature, and integrator/thermostat/barostat state. It validates first, creates a unique sibling temporary file, writes and syncs it, atomically renames it and syncs the parent directory. Failed pre-write validation preserves an existing file. A create-new collision does not delete a file we did not create.
- The v1 parser remains available, with required headers, finite/positive geometry, checked atom-buffer arithmetic, duplicate/trailing-header rejection and trailing-data rejection. V2 JSON rejects unknown/duplicate record fields, truncated/malformed/trailing data, inconsistent lengths, nonfinite/overflow/invalid mechanics and invalid component state at restoration.
- Default continuation checks normalized model/table contents, atom masses/species/charges/elements, topology, constraints, dt, ensemble/controller parameters and seeds/initialization identity; restores mechanical and stochastic/extended state exactly. Steps, cadences, output paths, restore/write paths and neighbor skin may vary. V1/mechanical-only files require explicit `restart_mode = "state_only"`. State-only intentionally retains loaded reduced positions/velocities/box/step with the new config's identity and freshly initialized controller state, projecting the new constraint manifold as required.
- `Simulation::metadata()` and `RunSummary.metadata` expose ensemble, kinetic DOF, constraints, conserved COM modes, output unit system and pressure status. Constrained pressure is explicitly `DiagnosticWithoutConstraintImpulseStress`. CLI labels this limitation and reports NVT/NPT energy changes as energy exchanged instead of a conservation fault.
- CLI `check` performs complete setup/identity validation through an output-disabled constructor. Real binary tests verify existing trajectory/log/restart sentinels are preserved and absent outputs are not created.

## Exact public/config APIs for Task 4/5

```rust
impl InputConfig {
    pub fn normalized(&self) -> modena_core::Result<Self>;
    // from_str/from_path call this; Simulation and builder call it too.
}
impl Simulation {
    pub fn from_config(config: &InputConfig, base_dir: &Path) -> Result<Self>;
    pub fn from_config_without_outputs(config: &InputConfig, base_dir: &Path) -> Result<Self>;
    pub fn metadata(&self) -> RunMetadata;
    pub fn write_checkpoint(&self, path: &Path) -> Result<()>;
}
pub struct RunSummary {
    pub steps: usize, // this invocation, while current_step is global
    pub final_energies: Energies,
    pub conservation: ConservationReport,
    pub metadata: RunMetadata,
}
pub struct RunMetadata {
    pub ensemble: Ensemble,
    pub kinetic_dof: usize,
    pub constraint_count: usize,
    pub conserved_com_dof: usize,
    pub pressure_status: PressureStatus,
    pub output_units: UnitSystem,
}
pub enum Ensemble { Nve, Nvt, EquilibrationNpt, CanonicalNpt }
pub enum PressureStatus {
    UnconstrainedPairVirial,
    DiagnosticWithoutConstraintImpulseStress,
}
// Metadata and enums implement Serialize; serde enum strings use snake_case.
// These types, and RestartMode, are re-exported from modena_core.

pub enum RestartMode { Continue, StateOnly }
pub struct RunSpec { /* existing fields */ pub restart_mode: RestartMode }
// TOML/JSON strings: "continue" (default), "state_only".

pub enum BarostatSpec {
    MonteCarlo {
        target_pressure: f64,
        target_temperature: f64,
        max_log_volume: f64, // default 0.05, dimensionless
        interval: usize,    // default 10
        seed: u64,          // default 42
    },
    /* existing Berendsen and rejected legacy ParrinelloRahman */
}
// TOML/JSON type = "monte_carlo". Other numeric fields use selected units.

pub trait Barostat: Send {
    fn apply(&mut self, system: &mut System, virial: &Virial, dt: f64);
    fn apply_with_model(&mut self, system: &mut System, model: &dyn ForceModel,
        neighbors: &mut NeighborList, pe: &mut f64, virial: &mut Virial,
        dt: f64) -> Result<bool>; // true means changed/accepted volume
    fn checkpoint_state(&self) -> Result<serde_json::Value>;
    fn restore_state(&mut self, state: &serde_json::Value) -> Result<()>;
}
impl MonteCarloBarostat {
    pub fn new(P: f64, T: f64, max_logv: f64, seed: u64) -> Self;
    pub fn with_interval(self, interval: usize) -> Result<Self>;
    // Public accepted/attempted counters; direct new defaults to interval 1.
    // Legacy energy callback attempt remains; production uses model hook.
}
```

`io::restart::Checkpoint` retains atoms/step/time/box_size/positions/velocities and adds `continuation: Option<ContinuationState>`. `ContinuationState` is a strict Serialize/Deserialize record containing forces, masses, species, charges, elements, constraint_count, conserved_com_dof, potential_energy, virial, physical_signature, integrator, thermostat and barostat. `None` means mechanical-only/v1 import; it never qualifies for full continuation. V2 starts with `MODENA_RESTART 2\n`, followed by a strict JSON record. `load`/`save` signatures are unchanged. Mechanical floats use serde_json's existing float_roundtrip support.

Keep normalization at the common boundary when adding typed import/special pairs. New config record fields/types must derive Serialize as well as Deserialize because validation/signatures inspect typed configs. JSON input is raw physical-unit config; `unit_factors` is not deserialized. Serializing an already-normalized real InputConfig while dropping its skipped marker is not a raw-input interchange workflow.

Manual `runner::init_velocities` keeps its existing signature/metadata behavior and removes COM. Its doc comment now requires callers to set constraint count and conserved COM DOF first (usually 3 for multi-atom NVE) to match runner initialization. Only the GPU ladder example currently calls this public helper; Task 5 should set physical COM metadata before using it in matched benchmarks. Langevin runner initialization uses its separate full-COM draw.

Minimal integration outside primary ownership: `crates/modena-mcp/src/lib.rs` RunSummary destructure adds `..` to compile with metadata; actual MCP metadata presentation/classifier/schema is Task 4. The intentional water NVT->NVE fixture and `bench/validation/water/modena_spce_nve_restore.toml` now explicitly set state_only. Restart fixture temp dirs include PID+atomic sequence to avoid parallel cleanup collisions.

## Observed RED before corresponding production fixes

```text
rtk proxy cargo test -p modena-core --test runner_hardening -- --nocapture
Initial: exit101, 0 passed / 9 failed.
DOF12 instead of9; initial constraint length/tangent not projected;
Langevin initial COM wrongly zero; stale PE -0.3872942341733775 vs fresh -0.3220323939714096;
raw real coordinates4.08 vs normalized1.2; unknown keys accepted;
changed physical identity accepted; infinite input caused neighbor overflow.

After isolating the unconstrained seam fixture:
rtk proxy cargo test -p modena-core --test runner_hardening stochastic_and_extended -- --nocapture
exit101: split variants1/2/3 differed (Langevin/CSVR/NHC); NVE matched.

rtk proxy cargo test -p modena-core --test barostat_hardening -- --nocapture
exit101: ideal gas volume mean0, independent expected5.

rtk proxy cargo test -p modena-cli --test check -- --nocapture
exit101: check replaced trajectory sentinel with an empty file.
Additional CLI diagnostic RED: NVT constrained run printed a Drifting
conservation verdict and omitted constraint-pressure diagnostics.

Legacy parser RED: accepted NaN time, negative box, trailing corruption,
extra magic-header token; maximum atom count panicked on 3*N overflow.
Invalid checkpoint-write RED: save accepted NaN time and overwrote old file.

Additional boundary RED/GREEN:
finite initial T1e308 accepted but overflowed KE and silently froze;
full-u64 JSON seed rejected by TOML inspection's signed-int limit;
state_only Langevin->VV checkpoint could not subsequently continue its COM DOF;
FCC replicas[65536,65536,1] panicked in u32 lattice multiplication.
Each dedicated regression observed failure before its corrective follow-up.
```

The first constrained FCC seam fixture encountered the expected periodic half-box tie and did not isolate continuation. It was replaced by an unambiguous XYZ dimer/free-atom fixture; no constraint tolerances were loosened. New model-hook rollback/cache/counter/corrupt-payload cases supplement the initial force-cache/MC/state behavioral regressions and exercise real components.

## Independent expectations and tested continuation domain

- Ideal gas uses N=4, P=T=1. Independently, p(V) is Gamma(N+1,kT/P), with mean5 and variance5. Three fixed seeds13/417/9091 run120,000 proposals each, burn20,000 each, pooling300,000 samples. Measured mean `4.997425867311`, variance `4.985568810455`; absolute tolerances0.12 and0.35. These are bounded regression gates including correlation, not a claim of independent samples or a universal mixing-time certificate. Proposal width0.8.
- Cache tests use actual LJ at r1.2, rc2.5 and skin0.7; 100 seeded attempts must contain both accepted and rejected moves. Every move compares forces/PE/virial with a fresh independent list; rejected state compares box/positions/velocities/forces/PE/virial exactly. A deliberate force-model error checks rollback and restored neighbor contents. Interval3 state restoration reproduces 23 subsequent model moves exactly, including final RNG/counters.
- Runner continuous17 vs split7+10 compares positions, velocities and forces with exact equality for NVE, Langevin, CSVR, NHC and MC NPT. Four-atom FCC, lattice2.4, dt0.003, T1; CSVR tau0.1/seed231, Langevin gamma0.7/seed771, NHC tau0.1/chain3, MC P0.1/width0.25/interval3/seed93 with CSVR. MC's interval phase crosses the seam.
- Constraint seams use three XYZ atoms in edge8, dimer distance1.2 plus a free atom; continuous25 vs split11+14 compares exact mechanics for VV, Langevin and NHC. Metadata confirms one constraint and diagnostic pressure. Initial unconstrained n4 has DOF9/KE4.5; constrained n4 has DOF8/T1 and projected tangent velocities; Langevin n4 retains DOF12 and nonzero initial COM.
- Existing 108-atom NVE restart tests retain exact checkpoint-state gates, relative logged seam/energy comparisons1e-12 and continuous-vs-split drift checks. Existing water64/SPME/CSVR1000 -> state_only NVE500 workflow passed; its historical seam tolerance1e-9 and NVE drift bound0.05 remain unchanged, with no stronger water conservation claim inferred.

Neighbor reference positions, pairs and force-model workspaces are rebuilt rather than serialized. Saved mechanics and RNG/extended variables round-trip exactly. The exact trajectory comparisons above use deterministic small-system enumeration under unchanged execution settings. Different cache/rebuild summation order, parallel resources, potential versions or floating-point platforms are not certified universally bitwise-identical; allowed neighbor-skin/output/execution changes can alter rounding. There is no new GPU/MPI execution or timing claim in this task.

## GREEN commands and integration findings

```text
rtk proxy cargo test --workspace
exit0: 228 passed, 0 failed, 1 ignored doc test.
Raw output: task-3-workspace.log in this evidence directory.
This full run includes runner18 before the final generator-capacity guard.

rtk proxy cargo test -p modena-core --test runner_hardening --test barostat_hardening --test config_builder --test io_and_lattice -- --nocapture
exit0: 43 passed (runner19, MC4, builder3, io/lattice17).
Final guard/mismatch checks covered here; task-3-focused.log.

rtk proxy cargo test -p modena-cli --test check -- --nocapture
exit0: 2 passed; existing and absent outputs plus real constrained/NVT output.
Raw output: task-3-cli.log.

rtk proxy cargo test -p modena-core --test restart water_two_phase -- --nocapture
exit0: 1 passed, 32.49s isolated.

rtk proxy cargo clippy --workspace --all-targets -- -D warnings
exit0, including final production guard and all targets.
rtk proxy cargo fmt --all -- --check
exit0.
rtk git diff --check
exit0.
```

Earlier full-suite integration exposed the old error-message expectation for run.log_every; precise field naming was restored and all17 io/lattice checks passed. Another parallel run lost nvt.csv while checkpoint/nve files survived; the exact isolated water workflow passed and the timestamp-only shared fixture naming was hardened with PID+sequence. The subsequent whole workspace passed. Clippy found the MSRV1.85 incompatibility of usize::is_multiple_of; modulo preserves interval behavior and the final default lint is clean.

## Remaining limitations / later ownership

Constraint impulse stress/full stress tensor remain absent; constrained pressure is diagnostic and constrained/anisotropic NPT stays rejected. MC supports unconstrained isotropic particles and requires consistent canonical thermal control; positivity of P/T/width/interval is enforced. Legacy direct piston/anisotropic APIs are experimental and unavailable through production config. Constraint counts assume physically independent distance constraints as documented by Task 1; exact duplicate pairs are now rejected, but there is no new general numerical constraint-rank analysis.

State-only deliberately changes protocol and resets thermostat/integrator/barostat state; stored mechanics are always reduced-unit values, so callers must keep an appropriate reference-unit interpretation when changing models. Full restart signatures deliberately include initial seed/temperature identity and model table contents; initialization cannot silently change during continuation. No arbitrary cross-version/platform bitwise continuation promise is made.

Typed molecular import, graph special-pair policy, Python interface breadth, MCP propagation/classifier and new examples remain Task4. Fresh whole-branch optional backend checks, release build, campaigns and matched timings remain root/Task5, including physical COM metadata for manual benchmark initialization. Original audit evidence remains untouched.

## Fix round 1 — four Important review findings

Base: `c5810faca4324b5ab3b6de8eeddff028c8122eb8`.
Follow-up commit: `17a6a9a6614ecb78eca743ee0106b2d8b5c61c6b`
(`:bug: Fix: bind restart identity and physical ensemble boundaries`).
Post-commit working tree is clean. Ownership stayed within
`src/runner.rs`, `src/config/validation.rs` and `tests/runner_hardening.rs` in
modena-core; no other production files or earlier-task APIs changed.

1. Full continuation now includes the normalized `config.electrostatics` value
   in `continuation.physical_signature.electrostatics`: `null` when absent, or
   the serialized `EwaldSpec` with `method`, `kappa`, `real_cutoff` and
   `kspace_cutoff`. This prevents old mechanical/force caches being restored
   across electrostatic Hamiltonian or evaluation-method changes. The existing
   common normalization boundary supplies this field. All six regressions cover
   removal, Ewald-to-SPME, each of the three parameters, and addition to a
   formerly pair-only checkpoint. An identical Ewald specification continues
   successfully. A v2 signature produced before this fix lacks the field and
   cannot certify full identity; it requires an explicit state-only import.
2. Runner construction checks normalized unordered atom pairs after explicit
   distance constraints and generated rigid-water triples have been combined,
   before setting constraint count/DOF. A repeated pair is a Config error even
   when target lengths agree. Six real-unit eight-water cases cover O–H1,
   O–H2 and H1–H2 overlap in either order, including same O–H lengths and
   contradictory lengths. The valid unmodified fixture has 24 constraints.
   General numerical constraint-rank analysis is still outside the supported
   claim; this fix is exact repeated-pair detection.
3. Deterministic state-only imports assess COM reduction from the new masses
   and loaded velocities, then reassess after any new constraint projection.
   The existing scale-relative zero-momentum criterion is
   `|sum(m_i v_i)| <= 64*EPSILON*sum(m_i |v_i|)` for more than one atom.
   Langevin keeps its stochastic COM modes. Validated full continuation retains
   saved COM metadata instead of reclassifying it. The eight-atom two-species
   rocksalt fixture starts with masses1/1 and zero momentum (<1e-13), then imports
   unchanged velocities with masses2/1, producing momentum >1e-5 and requiring
   COM reduction0/DOF24. Unchanged masses still use COM reduction3. The new
   deterministic protocol's continuous7 vs split3+4 full continuation preserves
   positions, velocities and forces exactly, including saved COM reduction0.
4. MC canonical pairing now explicitly requires controller T>0 and uses
   `|T_controller-T_MC| <= 1e-12*max(|T_controller|,|T_MC|)` without a unit-sized
   absolute floor. For CSVR, Langevin and NHC, matching T=1e-13 constructs
   successfully; zero controller versus MC1e-13 and controller2e-13 versus
   MC1e-13 reject. NHC already rejected zero on its own parameter boundary;
   the other five combinations were newly exposed RED cases. Public signatures
   and configuration syntax remain unchanged.

Observed TDD evidence (tests written before the production fixes):

```text
rtk proxy cargo test -p modena-core --test runner_hardening -- --nocapture
RED exit101: 19 passed, 4 failed.
- electrostatic identity changes accepted: [0,1,2,3,4,5]
- water overlaps lacked duplicate Config rejection: [0,1,2,3,4,5]
- changed-mass state-only COM metadata: actual3, expected0
- invalid near-zero MC pairings accepted: [(0,0),(0,1),(1,0),(1,1),(2,1)]
Raw: task-3-fix-round-1-red.log.
GREEN exit0: 23 passed, 0 failed; task-3-fix-round-1-green.log.

rtk proxy cargo test -p modena-core --test runner_hardening --test config_builder --test io_and_lattice --test restart -- --nocapture
exit0: 46 passed (runner23, builder3, io/lattice17, restart3).
Water protocol test included; restart suite completed in32.43s.
Raw: task-3-fix-round-1-covering.log.

rtk proxy cargo test -p modena-cli --test check -- --nocapture
exit0: 2 passed; task-3-fix-round-1-cli.log.
rtk proxy cargo clippy --workspace --all-targets -- -D warnings
exit0; task-3-fix-round-1-clippy.log.
rtk proxy cargo fmt --all -- --check
exit0.
rtk git diff --check
exit0.
```

The water RED fixture was corrected before production changes: one molecule
put its H–H distance across the small box's minimum-image limit, and an initial
zero-sigma trial hit the existing LJ constructor domain. Eight molecules with
the fixture's standard H sigma and zero epsilon isolate the reviewed overlap
behavior without changing either unrelated scientific domain.

The earlier full-workspace and MC moment evidence remains the evidence for
unchanged broader behavior. This review round reruns covering tests and default
lint/format checks only; final fresh whole-branch/backend suites remain root's
later gate. Neighbor/workspace rebuild and exact-seam limits documented above
are unchanged.

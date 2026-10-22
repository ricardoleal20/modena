# Task 4 report: typed chemical workflows and accurate client interfaces

Status: implemented; covering checks pass. Root owns independent review and later fresh full-branch/backend suites. No publication requested/performed. Original research/audit evidence and original checkout preserved. No subagents.

Task BASE: `17a6a9a6614ecb78eca743ee0106b2d8b5c61c6b`.
Implementation commit: `3bf7f54c683115aabc366d59f3d7d33b36669297` — `:sparkles: Add: typed chemical workflows and reliable clients`. Post-commit working tree clean.
Worktree: `/Users/ricardoleal20/.codex/worktrees/modena-scientific-hardening/modena`.
Branch: `ricardoleal20/modena-scientific-hardening`.
Read exact Task4 brief, scientific scope spec, chemical handoff and Task3 API/handoff. Applied test-driven-development and verification-before-completion. All shell invocations use RTK. Root explicitly expanded ownership to units.rs and charge conversion docs/tests after an independently detected physical conversion error.

## Delivered inventory

- `crates/modena-core/src/io/xyz.rs`: strict documented typed Properties records, per-atom numeric model type/mass/charge, positive/finite validation, coordinate finiteness, atom-capacity guard and incremental parsing without trusting file count for allocation. Existing `read_xyz` delegates to identity-aware reader and returns the same System type. Trajectory writer stays unchanged.
- `src/io/input.rs`: Serialize/Deserialize topology LJ/Coulomb three-distance arrays and explicit pair exceptions; charge conversion documentation corrected.
- `src/config/mod.rs`: shortest bond graph BFS limited to distance3, indexed normalized policies, explicit exception precedence, legacy endpoint exclusions, finite/nonnegative weights and conflicts/index validation. Sparse visited state avoids clearing N-sized storage for every BFS.
- `src/config/builder.rs`: raw real/reduced units selection, XYZ convenience method, explicit restart mode; shared normalization retained.
- `src/runner.rs`: converts typed file masses/charges and positions exactly once after normalized config, rejects typed global identity overrides, installs the single indexed neighbor special-pair policy shared by pair and long-range calculations. Rigid-water exclusions merge through indexed maps; incompatible explicit water intrapair exceptions reject. All Task3 validation/DOF/restart/pressure contracts preserved.
- `src/units.rs`: correct physical charge factor sqrt(C/(sigma*epsilon)), including its derivation/doc comments.
- `python/modena`: all production structure/model/controller/topology/constraint/electrostatic/barostat/restart dataclasses; robust TOML strings and namespaced recursive tables; absolute source/table/restore/output paths before temporary TOML relocation; numeric CSV values; subprocess errors retain exit/stdout/stderr; release/debug binary fallback. Two original parity tests ported to stdlib unittest with actual numeric steps/observable gates.
- `modena-mcp`: all supported production config schemas in embedded config_schema.json; capabilities/limitations/backend-source-presence tool; RunMetadata propagation, actual thermo unit labels, explicit reduced conservation units, barostat-aware NVE status, real accumulated post-initial finite PE means/count/windows, FIFO-independent aggregation, transition evidence uses actual means and counts. MC-coupled target-temperature sweeps explicitly reject.
- `crates/modena-core/tests/chemical_workflows.rs`: 14 real-component tests, independent graph/electrostatic energy/force expectations, ring shortest paths, malformed/unknown numeric types, global override rejection, capacity error, default legacy policies, water exception compatibility, typed/special policy continuation identity, old-convention checkpoint rejection, Rust builder normalization.
- `crates/modena-mcp/tests/chemical_workflows.rs`: 5 actual call_tool/request-parser workflows covering metadata/units/means, NVE classification, discovery/probe, MC sweep rejection, FIFO accumulation/transition.
- `python/modena/tests/test_workflows.py`: 8 new stdlib workflows; together with 2 legacy cases, 10 total. No dependency installation.
- `examples/chemical`: checked/runnable typed four-atom molecular NVT, v2 continuation, real ionic SPME, unconstrained MC NPT, fixed-volume real water; README states parameter conventions, units, APIs and scientific limitations. chain.xyz is intentionally tracked source; generated results and Python bytecode ignored by local .gitignore files.

## Exact public/config APIs and compatibility

```rust
pub fn io::xyz::read_xyz_with_identity<P: AsRef<Path>>(
    path: P, box_size: SimBox) -> Result<(System, bool)>;
// bool=true only for the supported explicit typed Properties declaration.
// Existing read_xyz<P>(path, box_size) -> Result<System> unchanged.

pub struct BondedTopologySpec {
    pub lj_scales: Option<[f64; 3]>,
    pub coulomb_scales: Option<[f64; 3]>,
    pub exceptions: Vec<PairExceptionSpec>,
    // Existing bonds/angles/dihedrals/exclude_bonded retained.
}
pub struct PairExceptionSpec {
    pub atoms: [usize; 2],
    pub lj: f64,
    pub coulomb: f64,
}
pub fn config::topology_pair_scales(
    spec: &BondedTopologySpec, n: usize) -> Result<Vec<PairScale>>;

impl SimulationBuilder {
    pub fn units(self, units: UnitSystem) -> Self;
    pub fn xyz(self, file: impl Into<String>, box_size: [f64; 3]) -> Self;
    pub fn restart_mode(self, mode: RestartMode) -> Self;
}
```

Typed XYZ uses the exact whitespace-separated comment token:
`Properties=species:S:1:pos:R:3:type:I:1:mass:R:1:charge:R:1`.
Following rows: `element x y z type mass charge`, exactly7 fields. Numeric type is zero-based and must refer to an actual model type; element strings are presentation labels, not an implicit chemical parameter database. The reader rejects unsupported Properties layouts, malformed indices/masses/charges, nonfinite coordinates/charges and nonpositive masses. Capacity exceeds u32 indexed atoms (or 3*N arithmetic capacity) return Input errors; truncated large files allocate only parsed rows.

Plain XYZ retains type0/mass1/charge0 defaults, explicit global structure.mass/charge override behavior, and X-label replacement by structure.element. Typed XYZ rejects both global mass and global charge, avoiding silent erasure of heterogeneous identity. For real typed XYZ the reference mass is explicitly1 amu. Positions/masses/charges use declared input units and convert once after common InputConfig::normalized; box already normalized. Absolute/relative source paths use the config base_dir through existing CLI runner plumbing.

TOML arrays `topology.lj_scales` and `topology.coulomb_scales` index shortest bond distances1,2,3, independently. Each omitted array defaults to all1 when exclude_bonded=false, or all0 when true. Bond graph derives all1-2/1-3/1-4 pairs without requiring explicit angles/dihedrals. Legacy explicit angle/dihedral endpoint exclusions remain when exclude_bonded=true and both arrays omitted. Explicit unordered exceptions override graph/legacy entries. Duplicate equal exceptions are compatible; contradictory duplicates, invalid indices/self pairs, nonfinite/negative weights reject. Weights may exceed1 (not artificially restricted). Water's three intramolecular exclusions remain mandatory; explicit incompatible water exceptions reject and compatible zeros succeed. One shared indexed neighbor policy supplies LJ and Coulomb paths, with full Coulomb bare-pair correction beyond the real-space cutoff; no duplicate model-local scaling.

Existing normalized restart signature naturally includes added serialized topology fields and typed masses/types/charges/elements; no identity weakening. Source Rust literals of BondedTopologySpec need new fields or ..Default::default(). Serializing a topology now includes the new default fields, so pre-change v2 topology signature shapes cannot establish default continuation equivalence; deliberate state_only remains available. Old real-charge-convention checkpoints are incompatible at sigma!=1 and reject normal continuation through unit factor/charge identity checks. State-only imports reduced mechanics into a new protocol, and do not migrate or certify old force/charge physics.

Python exports: Config, Fcc, Xyz, Water, RockSalt, RandomAlloy; LennardJones, Buckingham, Morse, SoftSphere, Tabulated; VelocityVerlet, Langevin, NoseHooverChain; Run, Csvr, Berendsen, NoThermostat; Ewald; Bond, Angle, Dihedral, DihedralTerm, PairException, Topology; Constraint, Constraints; MonteCarlo, BerendsenBarostat. Type dictionaries for model types/overrides remain supported and render [[model.types]]/[[model.overrides]] correctly. Topology dataclass lists render recursive [[topology.dihedrals.terms]] tables. Config positional structure/model/integrator/run order remains; added units defaults reduced, all optional sections default None. Config.to_toml performs no Python validation/unit conversion; Rust remains authoritative.

```python
check(cfg: Config, binary: Optional[str] = None, *, base_dir: Optional[Path] = None) -> str
run(cfg: Config, binary: Optional[str] = None, *, base_dir: Optional[Path] = None) -> RunResult
```

base_dir defaults current directory. Transport deep-copies and absolutizes structure/model file, restore_file, restart_file, trajectory_file, log_file before temporary config writing. check preserves output paths without creating/truncating them. run intentionally captures thermo in a temporary CSV (as before), preserving declared restart/trajectory outputs; caller config is not mutated. RunResult keeps steps_rows/conservation/perf/toml. step values are integers; time and observables are floats. Nonzero subprocess results raise RuntimeError including status/stdout/stderr. Missing binary/timeouts retain native subprocess exceptions. NVT conservation parsing can remain unknown because CLI correctly prints energy exchanged; use ensemble metadata through MCP when needed.

MCP get_capabilities() exposes CPU/Rayon production config availability, GPU/MPI source presence and optional feature names, unavailable production config selection and runtime_available=null/runtime_probe=not_performed. This does not claim a compiled/enabled device or distributed production run. Discovery includes all5 structure/model families, real/reduced units, NHC/CSVR/Langevin, Ewald/SPME, topology/scales/exceptions/constraints, canonical MC/equilibration Berendsen and restart modes. Physical coupling rules remain enforced by normalized Rust config and runner, not duplicated loosely in JSON schema. Legacy piston/anisotropic production support is absent from discovery and explicitly rejected by Rust.

RunStatus adds metadata:Option<RunMetadata>, thermo_units:UnitSystem, conservation_units:"reduced", mean_potential_energy:Option<f64>, sample_count:usize, sample_window JSON. ProbeResult adds metadata, thermo_units and conservation_units. ThermoSample adds thermo_units. Both raw and normalized real configs are labeled correctly: RunEvent::Thermo and RunSummary.final_energies already pass through runner.display_energies; they are real physical output values, not reduced, and are not converted a second time. System/checkpoint/conservation still use internal reduced mechanics/energy. Finished status uses initial global step + invocation steps for continuation.

Mean policy: arithmetic online mean of all finite potential energies from emitted post-initial samples in this invocation; excludes fresh initial step0 or restored seam sample, includes equilibration, removes no burn-in. Cadence is run.log_every; no extra final off-cadence point is inserted into the mean. sample_count counts finite PE samples; first/last step exposed. Latest events FIFO capped4096, mean retains the entire invocation independently. No latest-point fallback. NVE classifier includes absence of barostat and is refreshed from shared ensemble metadata. Target-temperature MC-coupled sweeps reject synchronously; initial_temperature changes only initial draw, not canonical targets. NHC target sweep not implemented and stated; thermostat target supports CSVR/Berendsen, Langevin target its matching integrator.

## Physical charge correction and independent expectations

At r*=r/sigma, U*=U/epsilon, requiring epsilon*q*²/(r/sigma)=C*q_e²/r gives q*²=C*q_e²/(sigma*epsilon), hence q*=q_e sqrt(332.06371/(sigma*epsilon)). Prior code multiplied by sigma, yielding physical energy/force error sigma². Reduced mode and sigma=1 unaffected; every real electrostatic system with sigma!=1 changes physically. Original research/baselines intentionally remain preserved; historical real-water/ion energies cannot be treated as corrected validation evidence.

Analytic fixture: r4 Angstrom, q_i=+1e, q_j=-2e gives energy -166.031855 kcal/mol, radial force -41.50796375 kcal/(mol Angstrom). Three reference choices (sigma,epsilon,mass)=(2,0.5,1),(3.4,0.238,39.948),(1,1,23) reproduce both within1e-10. Actual Ewald and SPME runner test subtracts identical electrostatic systems with Coulomb pair weights1 versus0, leaving independently expected bare removal despite mesh/self/background terms; physical energy and force agree within1e-10. Real typed Na/Cl fixture checks converted coordinates (positions1,2,3 and3,2,3), masses23/35 with 1amu reference, species0/1 and independently derived charge factor.

Four-atom chain: reduced positions x2,5,8,11, box24, charges+1,-1,+1,-1; LJ sigma2/epsilon0.5/cutoff8. Zero-stiffness bonds isolate nonbonded policy. LJ weights[0,0.2,0.5], Coulomb[0.1,0.4,0.8], explicit 0-3 weights0.7/0.6. Oracle adds sum[(lj-1)U_shifted_LJ(r)+(c-1)qiqj/r] to unscaled total; r9 pair is beyond both real electro cutoff4 and LJ cutoff8, so it contributes only the independent Coulomb correction. Both Ewald/SPME totals agree within1e-10. Ring0-1-2-3-0 asserts direct0-3 shortest distance1, cross0-2 distance2. Legacy chain with only bonds and exclude_bonded=true has zero nonbonded energy.

## Observed RED and GREEN evidence

All raw logs below are in this same ignored evidence directory. RED failures were observed before corresponding production fixes; initial compile/test-fixture API errors were corrected before recording behavioral RED. The first LJ oracle erroneously called low-level energy_force beyond its cutoff; it was corrected to apply the physical cutoff explicitly. The initial constrained FCC probe was changed from one-cell half-box ties to replicas2 to isolate metadata; no projection tolerances or physics checks were weakened.

1. `rtk proxy cargo test -p modena-core --test chemical_workflows -- --nocapture`
   - task-4-rust-red.log: exit101,0/4 pass; typed identities [0,0] vs[0,1], negative typed mass accepted, lj_scales config rejected, graph chain omitted distance2/3 exclusions.
   - task-4-rust-green.log: exit0,4/4 initial pass (before expanded physical-charge correction).
2. `rtk proxy cargo test -p modena-core --test chemical_workflows real_coulomb -- --nocapture`
   - task-4-charge-red.log: exit101 at sigma2, incorrect physical energy independently expected -166.031855.
   - task-4-charge-green.log: exit0,5/5 after correction, including initial typed/graph tests with corrected independent charge expectation.
3. Supported Python invocation below, initial task-4-python-red.log:6 cases,2 failures+2 behavioral errors (typed classes missing, numeric CSV strings, unescaped quote/relative source TOML). task-4-python-green.log:6/6 pass. Later combined paths/restart/water/table and DEL coverage culminate in10/10 final pass.
4. `rtk proxy cargo test -p modena-mcp --test chemical_workflows -- --nocapture`
   - task-4-mcp-red.log:exit101,0/3; missing capabilities, barostat true NVE, null mean/unit metadata.
   - task-4-mcp-sweep-red.log:exit101; coupled MC target sweep incorrectly returned success/run ids instead of explicit synchronous limitation.
   - task-4-mcp-green.log:exit0,5/5 after fixes and FIFO evidence.
5. task-4-compatibility-red.log:exit101,10/12; legacy angle endpoint missing and explicit nonzero water exception silently discarded. Final chemical log passes both with graph/precedence safety.
6. task-4-xyz-count-red.log:exit101; usize::MAX atom count produced capacity panic. Final14-case chemical log returns error without panic.
7. task-4-quoting-red.log: Python10 cases,1 error; TOML parser rejects literal DEL. Final Python10 cases passes with unicode escape.

## Fresh covering commands/results

```text
rtk proxy cargo test -p modena-core --test chemical_workflows --test pair_scales --test config_builder --test io_and_lattice --test electrostatics --test runner_hardening --test water --test restart -- --nocapture
exit0:77 passed,0 failed. Chemical12 at this run; remaining65 existing covering tests.
Raw task-4-rust-final.log. Includes water NVT->state_only NVE restart26.87s and rigid-water stack68.05s.

rtk proxy cargo test -p modena-core --test chemical_workflows -- --nocapture
exit0:14 passed,0 failed, latest parser/compatibility/charge/restart/builder checks.
Raw task-4-chemical-final.log.

rtk proxy cargo test -p modena-mcp -- --nocapture
exit0:22 passed (16 unit,5 chemical request-path,1 reproducible sweep).
Raw task-4-mcp-final.log.

rtk proxy cargo test -p modena-cli --test check -- --nocapture
exit0:2 passed, output sentinels and constraint/NVT presentation preserved.
Raw task-4-cli-final.log.

rtk proxy cargo build -p modena-cli
exit0, refreshed debug binary for Python/end-to-end paths.
rtk proxy /Users/ricardoleal20/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/bin/python3 --version
Python3.12.14, supported runtime; no pytest installation.
rtk proxy /Users/ricardoleal20/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/bin/python3 -m unittest discover -s python/modena/tests -v
exit0:10 passed,0 failed (2 legacy+8 new), task-4-python-final.log.

rtk proxy cargo clippy --workspace --all-targets -- -D warnings
exit0, task-4-clippy.log.
rtk proxy cargo fmt --all -- --check
exit0.
rtk proxy git diff --check
exit0.
```

Additional earlier narrower raw logs: task-4-rust-covering.log45 pass; task-4-python-covering.log9 pass; task-4-mcp-covering.log22 pass. These are retained transparently; final evidence above is authoritative for latest corresponding paths. No complete workspace/backend campaign claimed here; root handles that final gate.

## Actual examples and measured sample domains

`task-4-examples.log` retains every argv/stdout/stderr/exit for sequential CLI check/run of molecular_nvt, continue_nvt, ionic_spme, mc_npt, water_nvt: all10 actions exit0. Inputs are checked-in source; generated results stay ignored. These are smoke/workflow demonstrations, not equilibrated material studies or validated force-field parameter sets.

| Deck | Atoms | Invocation/global steps | Units | Last T | Last PE | Domain |
|---|---:|---|---|---:|---:|---|
| molecular_nvt |4|20 /20|real|78.286805 K|-0.658291 kcal/mol|typed charged chain, graph-weighted LJ+Coulomb, NHC|
| continue_nvt |4|10 /30|real|61.539265 K|-0.508333 kcal/mol|full same-identity v2 continuation|
| ionic_spme |64|10 /10|real|100.619687 K|-6596.068128 kcal/mol|neutral rocksalt ions, corrected charges, Langevin+SPME|
| mc_npt |32|20 /20|reduced|1.027592|-19.422753|unconstrained isotropic canonical MC+CSVR matching T|
| water_nvt |24|5 /5|real|97.930771 K|0.770521 kcal/mol|rigid fixed-volume SPC/E, CSVR+SPME, diagnostic pressure|

MCP real-unit request sample:8 post-initial samples, step1..8, independent arithmetic event mean60.958701026164 kcal/mol; reported mean60.95870102616357 kcal/mol. Shared metadata output_units=real, thermo_units=real; mean excludes initial and includes equilibration. This sample intentionally exercises real-unit conversion at reference sigma2/epsilon0.5; it is not a physical equilibrium study.

MCP two reduced NVE sweep runs, initial T0.7 and1.1,4200 steps with log_every1: counts4200 each; FIFO lengths4096 each. Accumulated PE means -2.9332380452019446 and -3.079379048255322, while latest PE values are -3.779749856400813 and -4.28518366204284. The transition heuristic slope equals abs(mean1-mean0)/0.4, not latest-point slope. Original floating-point NVE energy drift diagnostics remain visible for these longer samples; tests assert the estimator mechanics, not a new conservation/transition claim.

## Remaining limits and handoff

- Constrained pressure is diagnostic without validated impulse stress. No constrained/anisotropic NPT, piston promotion, general numerical constraint-rank analysis, electronic structure/reactive/many-body force fields, or GPU/MPI production input wiring added.
- Typed extended XYZ supports one exact documented Properties layout. Labels do not lookup element-specific parameter tables; numeric model type is explicit. Single-frame importer remains distinct from the plain trajectory writer, which does not round-trip typed metadata by itself. Type records require source import/checkpoint identity for molecular restart.
- Real electrostatic historical results with sigma!=1 require new simulations under corrected conversion. Reduced mechanics do not automatically migrate old Hamiltonians. Signature checks deliberately prevent silent old-convention continuation; state_only imports require appropriate user interpretation of reduced reference units.
- Added default topology signature fields may reject older v2 topology signatures; state_only is explicit available protocol change. Typed identities/policies participate naturally in full identity checks. Neighbor/correction caches rebuild, and no universal cross-platform/rounding bitwise trajectory promise added.
- MCP means are unweighted emitted-sample arithmetic means including equilibration; no confidence interval/autocorrelation/burn-in inference. Cadence, invocation boundaries and finite-sample counts are explicit. Transition detector is a largest finite-difference heuristic only. Event FIFO excludes initial/seam sample and retains at most4096 post-initial events.
- MCP backend device/runtime presence is not probed; GPU/MPI source presence does not mean active hardware/runtime/production capability. Cost estimates remain rough heuristics.
- Python run still captures its log in a temporary file rather than publishing the configured log path (existing transport behavior). It preserves intentional restart/trajectory destinations and check sentinels, returns numeric captured rows, and exposes all current production configuration sections. Standalone Python source path defaults cwd or explicit base_dir, unlike config-file CLI paths.
- Normalized Rust InputConfig serialization drops its conversion marker and is not a raw real-config interchange mechanism (Task3 contract); MCP JSON callers should supply raw physical config values. NHC target sweep and coupled MC target-temperature sweep are explicitly limited.

Root should review the physical charge compatibility impact, sample-unit evidence and policy precedence, then run fresh full branch/default/optional backend tests as planned. No root-final benchmark or acceptance campaign results are inferred from these smoke examples.

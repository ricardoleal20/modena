# Final whole-branch review — scientific hardening

**Verdict: changes required before local branch delivery.** Two Important findings remain; neither requires a new engine or a broader performance campaign. One Minor discovery inconsistency is also recorded. No Critical findings. This is a local-delivery review, not approval to merge or publish.

## Scope and review basis

- Worktree: `/Users/ricardoleal20/.codex/worktrees/modena-scientific-hardening/modena`.
- Base: `b68abd112339aa6970261f470a572da6633bc617`; reviewed HEAD: `10bd9ee3077509ee472875f6697fd6d1a849753e`. Final read-only Git checks confirmed this HEAD and an empty tracked/untracked-short status.
- Binding spec and five-task plan: `docs/superpowers/specs/2026-10-03-scientific-hardening.md` and `docs/superpowers/plans/2026-10-03-scientific-hardening.md`; ledger, original findings, task reviews/rereviews and focused implementation reports were inspected as evidence, not treated as proof of correctness.
- Read the complete source/docs-only review diff in passes, including all current implementation, tests, tools and documentation hunks outside the archived evidence directory. Used the complete package/index for scope. Historical byte-preserved outputs were not regraded. Named outside-diff checks covered the existing lattice population constructor and MCP lifecycle/discovery calls reached through the expanded input surface.
- No production source, index, HEAD, branch or tracked documentation changed. No helpers, subagents, Cargo builds, suites or timing campaigns run by this reviewer. The coordinator explicitly authorized two small standalone Rust probes and their scratch evidence below.

## Spec alignment and strengths

The implementation substantially follows the accepted scope: classical pairs, ions and existing fixed-volume molecules; constrained dynamics with explicit pressure limits; canonical unconstrained isotropic MC NPT; strict production rejection of unsupported pressure coupling. It does not substitute universal LAMMPS parity or new reactive/electronic/HPC features for the requested hardening.

The strongest changes are:

1. Sampling corrections have analytical expectations, multiple seeds and real component tests. CSVR uses the corrected chi-square factor and signed rescaling; constrained BAOAB separates stochastic velocity from drift corrections, and rigid reference-direction projections address the previously observed rotor cooling. Constraint failure now propagates.
2. Force/energy/state contracts are more coherent: joint Hermite tables, strain-derived electrostatic virials, separate indexed LJ/Coulomb policies, shared reciprocal corrections, owned SPME workspaces, explicit unsupported-backend errors, and refreshed mechanical caches after box moves and replica swaps.
3. Full-state restart seams include stochastic/extended state and physical identity. Intentional protocol changes are explicit. Typed atom import, dimensional charge correction, graph-distance policy and common normalization connect through Python and MCP without a second client-side physics implementation.
4. The real-charge ruling is physically justified: `q*=q*sqrt(C/(sigma*epsilon))` follows from `r*=r/sigma` and `U*=U/epsilon`. Independent physical energy/force tests and fail-closed identity checks support the correction. Documentation clearly retracts the earlier water/unit claim and does not promise historical migration.
5. Current evidence is described with its actual limitations: failed coarse-mesh validation is preserved; corrected production units are distinguished from historical results; matched CPU timing scopes, resources, RSS and the remaining performance gap are explicit. The GPU wrapper fix retains failure metadata. Criterion teardown scope and Plotters fallback are disclosed.

The remaining COM issue is a cross-task exception to the otherwise improved sampling/DOF contract. The alloy issue is a reachable older constructor defect that the newly promised shared input boundary does not yet contain.

## Findings

### Important I1 — SPME is assigned a conserved COM manifold it does not preserve

**Primary location:** `crates/modena-core/src/runner.rs:304–309`.
**Other required fix paths:** `runner.rs:317–322`, `runner.rs:400–403`; misleading force comment at `crates/modena-core/src/spme.rs:639`. Thermostat consumers include `thermostat.rs:182–186` and `nose_hoover.rs:150–155`.

The runner subtracts three COM degrees of freedom whenever the integrator is not Langevin. Its SPME model analytically differentiates the finite B-spline mesh energy, whose value changes under a continuous common translation relative to the grid. That model therefore has a nonzero total force even in exact arithmetic. The sum of the four derivative weights being zero establishes partition-of-unity behavior, not continuous translation invariance of this mesh energy.

The focused probe reused the already accepted neutral three-charge validation fixture (box9; charges1,-.7,-.3; kappa.8, real cutoff4.2). On the exact recorded default release rlib:

| Model | Summed force (x,y,z) | Independent common-x translation derivative `-dU/dx` |
|---|---|---|
| SPME kcut4.8 | (-.002196221738, -.001628763985, -.001085552043) | -.002196221729 |
| SPME kcut6.4 | (.000145507097, .000916990311, -.000144074430) | .000145507101 |
| Direct Ewald kcut6.4 | components below3.1e-16 | reference momentum-conserving model |

The same fixture through `InputConfig` and `Simulation`, VV dt.001,100 steps, initially zero total momentum, reports `kinetic_dof=6, conserved_com_dof=3` throughout. Its final SPME momentum is (-.000185287727,-.000159710518,-.000123094149); the Ewald control remains below2.2e-15. Thus the removed momentum manifold is not invariant, including at the mesh resolutions accepted by the numerical suite. This is an actual mismatch between the force model and the new DOF policy, not evidence that the mesh force is inconsistent with its energy.

**Impact:** incorrect temperature/kinetic-DOF metadata and CSVR/NHC target DOF for supported charged SPME systems. For this three-atom unconstrained mesh Hamiltonian the kinetic dimension is9, while the runner uses6. This also reaches fixed-volume water and MC NPT with deterministic/canonical controllers. Current continuation/state-only logic can reinstall COM3 even if fresh initialization is fixed.

**Practical fix:** make conserved-COM eligibility depend on both the dynamics/controller and the force model. The minimal consistent policy for the present SPME energy/force is COM0, initialize its full kinetic space, and use `3N-constraints` DOF. Propagate the same eligibility through fresh construction, both restart modes and identity checks; an old unsupported COM3 checkpoint must not silently override the corrected policy. Correct the translation-invariance comment and relevant capability/DOF documentation. Do not merely subtract the summed mesh force while retaining the existing mesh energy, since that would break the energy derivative contract just hardened.

**Regression boundary:** use a nonsymmetric charged fixture, establish the nonzero mesh total force/common-translation derivative, and assert the corresponding fresh VV/CSVR/NHC metadata, target DOF, constrained behavior and restart handling. Preserve the momentum-conserving Ewald/pair COM3 controls and Langevin COM0 behavior. A long new liquid campaign is unnecessary to demonstrate this structural defect.

### Important I2 — Valid alloy concentrations pass validation and panic during construction

**Primary location:** `crates/modena-core/src/lattice.rs:102–107`.
**Integration locations:** `crates/modena-core/src/config/validation.rs:147–153`, `crates/modena-core/src/runner.rs:159`, `crates/modena-mcp/src/lib.rs:698–705`.

A valid `random_alloy` FCC structure with replicas[1,1,1] has four atoms. Supply six positive, normalized concentrations `[.13,.13,.13,.13,.13,.35]`, positive masses and six corresponding LJ model types. The new strict `InputConfig::from_str` accepts it. The constructor independently rounds each of the first five populations: `round(.13*4)=1`, so five labels have already been emitted when it calculates the last population as `4-5`. In release this underflows and `labels.extend` panics with `capacity overflow`; debug arithmetic would also panic. The later `labels.resize(n,...)` guard is unreachable.

The focused actual release `Simulation::from_config_without_outputs` probe catches and records precisely that panic after successful `InputConfig` parsing. This requires no huge allocation or malformed data. It is a supported small chemical configuration through the newly expanded Python/MCP surface. The synchronous path cannot return its advertised validation error; the MCP async worker only handles `Result::Err` here, so this panic also bypasses its Failed-state transition.

**Practical fix:** apportion the integer populations with a bounded deterministic algorithm whose total is exactlyN (for example floor quotas plus largest remainders), then apply the existing seeded shuffle. The last-species remainder must never underflow, and legitimate small mixtures should not need to be rejected just to avoid the panic. Keep actual input/capacity failures as `Result` errors.

**Regression boundary:** the exact six-species/four-atom fixture plus zero fractions and several populations rounding upward; assert exact atom count, valid species identities, deterministic seed and successful common-config construction. Include a real MCP construction/run path so completion/error lifecycle remains observable. The pre-existing lattice helper is explicitly included because the accepted hardening promises an error-safe common configuration boundary; this is not a request for a new model.

### Minor M1 — Potential discovery still lists only Lennard-Jones

**Location:** `crates/modena-mcp/src/lib.rs:364–370`; advertised tool at `crates/modena-mcp/src/tools.rs:23–25`.

`list_potentials` still returns only `lennard_jones`, although the expanded production schema/runner and Python surface support Buckingham, Morse, soft-sphere and tabulated models. An MCP user following the dedicated discovery tool receives an incomplete inventory even though the new schema exposes the other families.

**Practical fix:** update this inventory and parameter names from the supported model families, ideally keeping it tied to the same capability/schema source. A request-path assertion for all five families is sufficient. This is a nonblocking discovery inconsistency by itself; it can be fixed with the two Important items or explicitly tracked.

## Focused verification and retained reproduction

Artifacts are under `.superpowers/sdd/2026-10-03-scientific-hardening/final-probe/`:

- `spme_com.rs`, `charges.xyz`, `run.stdout`, `run.stderr` preserve the mesh/COM reproduction and Ewald control.
- `alloy.rs`, `alloy.toml`, `alloy-run.stdout`, `alloy-run.stderr` preserve accepted configuration and caught capacity-overflow panic.
- Compilation logs, `provenance.txt` and `manifest.json` preserve commands, file hashes and artifact identity.
- Both probes used `target/release/libmodena_core.rlib`, SHA256 `794cca85de5b6eae1dd3ec9ef80387b4b5f905254991a758f8460765d107b295`, verified before/after and matching `docs/scientific-hardening/evidence/raw/original-probes-current.json`. The coordinator confirmed unchanged production source since that artifact. Compiled with `rtk proxy rustc --edition=2021 -O -C lto=thin -C codegen-units=1 ... --extern modena_core=target/release/libmodena_core.rlib -L dependency=target/release/deps`; run with `RAYON_NUM_THREADS=1`. No timing claim is made.

Prior reviews resolved their concrete findings, including rigid-rotor cooling, restart electrostatic identity/COM import/water duplicate checks/temperature comparison, and GPU launch-error metadata. Their code changes were inspected in the whole-branch pass; they do not cover the two reproductions above.

## Cannot-verify and deferred triage

- **Pending coordinator gates:** the coordinator reports fresh workspace252 pass/0fail/1existing ignored, format/default Clippy,10 tooling tests and7 GPU adapter tests without skips at this frozen HEAD. These are coordinator results, not suites independently rerun here. Remaining final feature/release/Python/numerical/MPI/example checks were held for the combined fix to avoid redundant campaigns. Do not claim the complete final gate set has passed yet; after fixes, preserve fresh results for the new HEAD and rerun affected verification.
- **Evidence campaign review:** the Task5 independent reviewer recomputed24 successful matched processes, medians/RSS/deck hashes,67 numerical observations and46 archived hashes. I inspected the tooling, contracts and relevant review reports rather than repeating that campaign or regrading historical raw bytes. These reports support provenance review; the new focused findings show why passing existing numerical rows is not a universal correctness claim.
- **Accepted scientific deferrals:** constraint impulse pressure remains diagnostic; constrained/anisotropic NPT and the legacy production piston fail closed. No validated converged bulk water/LAMMPS parity, broad molecular parameter matrix, corrected real ionic/water cross-engine campaign, or historical charge-checkpoint migration is claimed. General constraint-rank analysis remains an explicitly recorded limit. None is silently upgraded to validated production here.
- **Accepted execution/performance limits:** optional GPU/MPI remain limited library backends, with unsupported pair-policy/box contracts rejected. Resident GPU/distributed production is outside scope. Current CPU gap1.95–2.65x and N20000 RSS427MB versus22MB are observed open limits, not blockers demanding speculative optimization. Criterion destruction overhead and plotting fallback are disclosed; no source warning is attached to intentional historical-log whitespace.

## Local delivery readiness

After I1 and I2 are fixed with failing-before/passing-after regressions, request a focused independent rereview of the fixes and their restart/client interactions, then finish the coordinator's final gates on the resulting HEAD. M1 is minor and should be fixed or tracked explicitly. Preserve the original/current probe evidence and update DOF claims if the correction changes existing numerical baselines. No PR, merge, publication or additional engine work is required for this acceptance decision.

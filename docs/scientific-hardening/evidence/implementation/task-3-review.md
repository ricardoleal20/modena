### Spec Compliance

- ❌ Issues found: full continuation does not bind electrostatic configuration (`crates/modena-core/src/runner.rs:941`); implicit water constraints can be counted twice (`crates/modena-core/src/runner.rs:285`); state-only imports can retain invalid COM reductions after changing masses (`crates/modena-core/src/runner.rs:307`); the MC consistency check accepts a zero-temperature controller for a positive MC temperature (`crates/modena-core/src/config/validation.rs:306`). These violate physical identity, physical DOF, and canonical-temperature requirements.
- ⚠️ Cannot verify from this task diff: unchanged integrator/thermostat checkpoint implementations and replica-exchange cutoff behavior belong to earlier tasks; final whole-branch review should retain their earlier evidence. The full workspace log precedes the final generator guard, as disclosed; the focused final runner/config/io tests cover that guard. No whole-branch/backend certification is implied.

### Strengths

- `crates/modena-core/src/barostat.rs:377` implements the requested symmetric log-volume acceptance with the N+1 Jacobian. The real model hook refreshes caches and restores mechanics/neighbors on rejection or error (`:455`, `:490`). The analytical Gamma mean/variance regression uses three fixed seeds (`crates/modena-core/tests/barostat_hardening.rs:7`).
- `crates/modena-core/src/io/restart.rs:73` validates before saving, uses a unique create-new sibling, syncs contents, renames atomically and syncs the directory (`:113`). Strict record parsing and legacy checks accompany the format transition.
- `crates/modena-core/src/runner.rs:151` shares scientific construction with the output-disabled path; `crates/modena-cli/tests/check.rs:2` exercises the real executable against existing and absent outputs. Pressure limitations are explicit in shared metadata (`crates/modena-core/src/runner.rs:477`).
- `crates/modena-core/tests/runner_hardening.rs:163` and `:197` compare actual continuous/split mechanics across stochastic, extended-state, MC and constrained integrations. Configuration tests exercise the common boundary (`:98`).

### Issues

#### Critical (Must Fix)

- None found within this task scope.

#### Important (Should Fix)

1. **Include electrostatics in restart physical identity.** `crates/modena-core/src/runner.rs:941` serializes the pair model, integrator, thermostat, barostat, topology and constraints but omits `config.electrostatics`. Consequently adding/removing electrostatics, switching Ewald/SPME, or changing kappa/real cutoff/reciprocal cutoff passes the identity comparison at `:370` when all other fields match. Restoration then overwrites freshly evaluated forces/PE/virial with values from the old Hamiltonian (`:383`), contaminating the first kick and silently continuing a different model. Include the normalized electrostatic specification and add mismatch regressions for method, parameters, and presence/absence. Existing identity cases at `crates/modena-core/tests/runner_hardening.rs:431` do not cover this.

2. **Reject duplicate constraints introduced by generated rigid water.** `crates/modena-core/src/runner.rs:285` combines explicit bonds with generated water triples, then uses the combined count at `:290`. Duplicate validation at `crates/modena-core/src/config/validation.rs:321` sees only explicit bonds. A water deck adding an explicit O–H constraint already supplied by its water triple (with the same length) therefore counts that distance twice, lowers kinetic DOF, and rescales/thermostats at the wrong physical temperature. The focused dependency check confirms `Constraints::count` simply sums bonds plus three per water (`crates/modena-core/src/constraints.rs:79`) and parameter/system validation does not detect bond/water overlap (`:83`, `:109`). Validate the expanded distance-pair set before setting DOF; reject both same-length duplicates and contradictory lengths, with a generated-water regression. This is exact duplicate detection, not the deferred general constraint-rank analysis.

3. **Reassess COM reduction when state-only changes atom masses.** `crates/modena-core/src/runner.rs:307` copies `conserved_com_dof` from every v2 checkpoint, including state-only imports using the new configuration’s masses. With an old zero-momentum state and nonuniform changes to per-species masses, retained velocities have nonzero momentum under the new masses, yet the code still subtracts three DOF. Reported temperature and thermostat kinetic targets then incorrectly include COM kinetic energy in a reduced mode count. For state-only imports, establish COM metadata using the resulting masses and velocities (or explicitly reject incompatible reuse); preserve stored metadata for validated full continuation. Add a multicomponent mass-change test and its subsequent full-continuation seam.

4. **Require a strictly positive canonical controller temperature and compare on its scale.** `crates/modena-core/src/config/validation.rs:306` uses tolerance `1e-12 * max(abs(T_controller), abs(T_MC), 1)`. A VV/CSVR controller at T=0 and MC T=1e-13 passes this test, even though the MC constructor accepts the positive temperature and the controller continuously cools toward zero. This violates the stated positive-temperature canonical pairing. Explicitly require controller T>0 and use a relative comparison without the unit-sized absolute floor (or exact equality for identically normalized inputs). Add zero-versus-small-positive and unequal-small-positive regression cases.

#### Minor (Nice to Have)

- None added.

### Checks and Evidence

- Read the supplied complete diff once in consecutive manageable ranges, together with the task brief, implementer report, reviewer method and scientific scope spec. No source/index/HEAD changes or broad codebase crawl.
- Named focused dependency check: determine whether the existing constraint validator already prevents explicit-bond/generated-water overlap. Read only `Constraints::count`, parameter validation and system validation in `crates/modena-core/src/constraints.rs`; they do not.
- Inspected provided workspace/focused/CLI logs for test outcomes and warnings/errors. Workspace results are recorded at `.superpowers/sdd/2026-10-03-scientific-hardening/task-3-workspace.log:103` and subsequent test-result lines; final focused results show 4+3+17+19 passing at `task-3-focused.log:12`, `:21`, `:44`, `:69`; CLI shows two passing at `task-3-cli.log:11`. No warning/error lines were found in these logs. The ideal-gas measured moments appear at `task-3-focused.log:9`.
- No completed suites were rerun. Findings follow directly from the reviewed branches and signature fields; the new regression cases above are recommended for fixes.

### Assessment

**Task quality: Needs fixes.** The main mechanics, serialization and test structure are substantial and coherent, but the restart-identity and DOF gaps can silently produce scientifically inconsistent results. Correct the four boundary cases before accepting this task gate.

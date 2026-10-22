### Finding Verdicts

1. **Include electrostatics in restart physical identity — ADDRESSED.** `crates/modena-core/src/runner.rs:938` includes the complete normalized optional electrostatic specification in the saved physical signature. Existing equality validation therefore rejects presence, method and parameter changes before restoring old caches. `crates/modena-core/tests/runner_hardening.rs:460` covers all six mismatch cases and successful continuation with identical settings. Older signatures lacking the field fail closed, with state-only import retained as the explicit escape hatch.

2. **Reject duplicate generated-water constraints — ADDRESSED.** `crates/modena-core/src/runner.rs:291` checks the fully combined constraint set before setting DOF. The helper at `:949` expands all three rigid-water distances, canonicalizes unordered atom pairs and rejects any repetition independently of target length. `crates/modena-core/tests/runner_hardening.rs:502` covers both orientations of O–H1, O–H2 and H1–H2 overlaps, including matching and conflicting lengths, while checking the valid water constraint count.

3. **Reassess COM after state-only mass changes — ADDRESSED.** `crates/modena-core/src/runner.rs:310` retains stored COM metadata only for full continuation; state-only imports evaluate momentum using their resulting masses and velocities. The helper at `:973` applies the existing relative zero-momentum criterion, and `:395` reassesses after constraint projection. Langevin retains stochastic COM modes. `crates/modena-core/tests/runner_hardening.rs:552` verifies changed-mass nonzero momentum produces DOF24/COM reduction0, unchanged masses retain reduction3, and the new protocol subsequently continues with exact positions, velocities and forces.

4. **Require positive, consistently scaled MC/controller temperatures — ADDRESSED.** `crates/modena-core/src/config/validation.rs:307` explicitly rejects nonpositive controller temperature and compares relative to the temperature magnitudes without the unit-sized floor. `crates/modena-core/tests/runner_hardening.rs:597` checks matching small-positive temperatures and rejects zero or unequal small-positive controller targets across CSVR, Langevin and NHC.

### New Breakage in the Fix Diff

None found. The changes are confined to the four reviewed scientific boundaries and their behavioral regressions.

### Out-of-Scope Observations

None added. Earlier task-scope limitations and final whole-branch verification remain as recorded in the first review.

### Checks

- Read the complete supplied fix diff once, the scoped re-review method, previous four findings and appended fix report. Reviewed only these findings and fix-induced risks; no source/index/HEAD mutation, subagents or broad codebase reads.
- Inspected the supplied GREEN and covering logs. All four named regressions passed (`task-3-fix-round-1-covering.log:56`, `:57`, `:63`, `:66`), with 46 total covering tests passing (`:10`, `:33`, `:42`, `:71`). The CLI log records two passing tests (`task-3-fix-round-1-cli.log:10`), and the clippy log records successful completion (`task-3-fix-round-1-clippy.log:6`). No warning/error lines were found in those logs. These paths are under `.superpowers/sdd/2026-10-03-scientific-hardening/`.
- No completed suites were rerun; no new focused probe was needed because the fixes and reported behavioral coverage resolve the specific doubts.

### Verdict

**Fix round: All findings addressed, no new Critical/Important breakage.**

**Task-scoped spec compliance: Approved. Task quality: Approved.**

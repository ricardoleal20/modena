# Final combined integration fix — binding requirements

Worktree /Users/ricardoleal20/.codex/worktrees/modena-scientific-hardening/modena.
BASE10bd9ee3077509ee472875f6697fd6d1a849753e. Read final-review.md first: its
Important I1/I2 and Minor M1 are the full findings. Binding spec/plan are under
docs/superpowers; completed-task reports and final-probe evidence remain here.

Ownership: runner.rs, model/COM contract helpers strictly needed (ForceModel,
SPME/Ewald/composite as appropriate), init_velocities/DOF/checkpoint identity
paths, lattice.rs and necessary validated construction callers, MCP
list_potentials/discovery, covering core/MCP tests, relevant current docs and
durable evidence. No new engine or speculative performance rewrite. Do not
edit historical archived artifacts, unrelated optional backends, or
broaden chemical parameterization. Not alone; preserve all earlier changes.
No subagents. All shell commands begin rtk. No timing campaign/build concurrency
with timings: root has paused heavy work; no timing is active.

## I1 — model-aware conserved COM

Finite-mesh analytical SPME forces are derivatives of its translation-dependent
mesh energy. Preserve that derivative contract; don't subtract net force while
keeping old energy. Charged SPME must not claim COM3/3N-3-constraints: use COM0,
3N-constraints, initialize the full kinetic space and propagate eligibility
through fresh VV/CSVR/NHC, Langevin, fixed-volume constraints, continue and
state_only. Momentum-conserving pair/direct-Ewald controls retain COM3 where
initialized zero momentum is actually conserved. Handle zero-charge/no-mesh-
interaction cases consistently rather than inventing a removed manifold.
Current restart metadata must not reinstall unsupported COM3; fail closed for
incompatible continuation and reassess state_only using the actual new force
model plus existing momentum/mass logic. Update physical signature/metadata
validation as needed; report any checkpoint compatibility consequences.

Observe meaningful RED regressions against BASE, using the exact nonsymmetric
fixture in final-probe/spme_com.rs and reviewer table. Establish nonzero total
mesh force/common-translation energy derivative independently, and test
fresh/controller target DOF, constrained behavior and restart paths. Add proper
initial velocity/COM tests. Preserve Ewald/pair and Langevin controls and all
previous sampling/cache/restart tests; no weakened energy/force tolerances.
Correct spme.rs translation comment and DOF/capability documentation.

## I2 — bounded alloy populations

Exact case: FCC[1,1,1] N4, concentrations [.13,.13,.13,.13,.13,.35], six positive
masses and six LJ types. Shared parser accepts; current constructor panics.
Use deterministic bounded integer apportionment (e.g. floor plus largest
remainders) with total exactlyN, then existing seeded shuffle. Legitimate small
mixtures and zero fractions should work; no remainder underflow. Actual input/
capacity failures remain Result errors on common-config construction. Observe
RED then GREEN for exact atom counts, valid species, seed repeatability,
round-up/zero fractions and real Simulation + MCP lifecycle path.

## M1 — accurate potential discovery

Dedicated MCP list_potentials must list all five production pair families:
Lennard-Jones, Buckingham, Morse, soft-sphere, tabulated, with actual input
names/parameters consistent with schema. Add a real request-path regression.

## Verification/report

Run focused core/MCP covering tests, fmt and default Clippy as justified. Root
owns final fresh workspace/feature/Python/release/GPU/MPI/example/numerical gates;
don't duplicate broad campaigns. Original matched LJ24 processes and67 numeric
rows remain historical campaign evidence from clean11d7683; don't modify or
relabel them. Explain whether changes affect those workloads/scopes. Preserve
reviewer RED probe sources/raw/hashes and new RED/GREEN proof durably under
docs/scientific-hardening/evidence, update RESULTS limits/compatibility as needed.
Commit cohesive gitmoji+imperative Action, no AI attribution, no publish/merge.
Write final-fix-report.md here and durable copy, with exact commands/counts,
interfaces, compatibility, self-review and remaining concerns. Return only
status/HEAD/checks/concerns. One scoped independent re-review follows; this is
the single final combined fix dispatch, so resolve all three findings.

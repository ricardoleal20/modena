# Final combined integration fix

## Scope and status

Addressed Important I1 (charged SPME COM manifold), Important I2 (tiny alloy
population underflow) and Minor M1 (incomplete MCP discovery). Based on
`10bd9ee3077509ee472875f6697fd6d1a849753e`; the cohesive commit containing this
report identifies the final change. Exact validated source SHA256s and commands
are in [final-fix-manifest.json](../final-fix/final-fix-manifest.json).

Ownership was limited to the runner/model COM seam, public velocity helper,
alloy construction, MCP potential discovery, covering tests, current docs and
new evidence. No subagents, engine expansion, optional backend changes,
publishing, merges, pushes or timing campaigns.

## Behavior and interfaces

- Added `ForceModel::conserves_total_momentum(&System) -> bool`, defaulting to
  false so new providers must opt in. Pair and bonded forces and direct Ewald
  opt in. Long-range/composite providers propagate constituent eligibility.
  SPME opts in only when every atom charge is zero: any nonzero charge can
  couple to the fixed mesh, including mesh self energy. Empty/single-atom
  systems never acquire a removed three-mode manifold from runner setup.
- Runner COM eligibility now requires both a momentum-conserving model and a
  non-Langevin controller. Charged SPME uses COM0 and `3N-constraints` for VV,
  CSVR, NHC, Langevin and unconstrained isotropic MC NPT. Initialization draws
  the full kinetic space before tangent projection and final temperature
  scaling. Rigid charged water uses `6Nw` modes. Existing pair/direct-Ewald and
  neutral-mesh zero-COM controls keep COM3. No SPME force/energy formula or
  tolerance changed; no summed-force subtraction was introduced.
- Full continuation imports saved COM metadata only when eligible. Otherwise
  COM0 remains in force and the existing physical identity/DOF comparison
  rejects unsupported COM3. `state_only` checks actual new model/controller
  eligibility, then existing mass-weighted momentum logic, and rechecks after
  constraint projection. Physics signatures already bind electrostatics,
  charges, masses, controller, topology, seeds and initialization; no new
  signature key or blanket version bump was necessary.
- Public `runner::init_velocities` respects the caller's declared COM DOF.
  COM0 callers now draw COM motion, including a single atom. COM3 callers
  retain the previous seeded removal/rescaling arithmetic. Manual drivers
  remain responsible for setting valid constraint/COM metadata first.
- `random_alloy` retains its `Result<System>` signature. It validates finite
  fractions, positive masses, finite charges, valid lattice dimensions and
  bounded atom/type counts. It uses normalized floor quotas plus largest
  remainders, with input index resolving ties, then the existing seeded
  shuffle. All label/population/system buffer reservations are fallible.
  FCC coordinate order and the public infallible `fcc` signature are unchanged;
  the Result-based alloy path uses its validated/fallible construction helper.
- MCP `list_potentials` retains its request/response shape and derives names
  and actual input fields from the existing bundled model schema: LJ,
  Buckingham, Morse, soft-sphere and tabulated. Parameter order is schema map
  order; tests compare names/fields without depending on that order.

## Compatibility

Checkpoint format remains `MODENA_RESTART 2` and physical signature remains
`modena-physical-v2`. Old charged SPME COM3 checkpoints now fail normal
continuation; users can intentionally import mechanics with `state_only`, which
uses corrected DOF and fresh controller/RNG state. Pair/direct-Ewald checkpoints
and prior eligible SPME COM0 protocols (including Langevin) continue with the
same identity checks. Current SPME VV/CSVR/NHC/Langevin split continuations are
bit-identical in positions, velocities and force caches. Existing mass-change
state-only and electrostatic identity regressions still pass.

Alloy rounding can change species counts/seeded assignments for fractional
quotas relative to historical independent rounding. Seeds remain deterministic
under the corrected algorithm. For the exact N4 fixture, the six populations
are `[1,1,1,0,0,1]`; zero-fraction fixture gives `[0,2,2,0]`; the seven-species
multiple-round-up case gives `[1,1,1,0,0,0,1]`. Zero fractions receive no atoms.

Python/CLI/MCP physics continues to use the common core construction path.
MCP async invalid scientific input remains observable as Failed with an error;
valid alloy build/probe/async requests complete successfully.

## RED proof and independent expectations

The initial test-first run occurred before any production change: core
2 passed/8 failed; MCP 0 passed/2 failed. It reproduced COM3/6-mode charged mesh
setup, constrained DOF, unsupported COM3 continuation, single-atom/public
initializer loss of kinetic space, exact alloy underflow and incomplete
inventory. Logs are preserved in the new final-fix evidence directory.

After expanding controller/water coverage, production was temporarily replaced
with exact Git BASE bytes, all expanded regressions were run, and the final
production bytes were restored in a `finally` block. Source/test hashes and
commands are in [expanded RED provenance](../final-fix/final-fix-expanded-red-provenance.json).
Expanded BASE core: 2 passed/10 failed; MCP: 0 passed/2 failed. This proves the
additional real CSVR/NHC target consumption and charged rigid-water cases fail
without the fix. Baseline control tests pass intentionally.

The nonsymmetric neutral three-charge fixture exactly matches the reviewer:
box9, positions `(1.1,2.3,3.4)`, `(2.7,3.2,3.8)`, `(6.2,5.8,7.0)`, charges
`1,-.7,-.3`, kappa.8, real cutoff4.2. Independent common-x translation finite
difference (h1e-5) matches the summed analytical mesh force within2e-10, with
literal net-x references `-.002196221738` at kcut4.8 and `.000145507097` at6.4
within1e-11. Direct Ewald summed force stays below1e-14. These establish the
translation-dependent Hamiltonian; they do not require changing its derivatives.

Expected kinetic dimensions are independently `3*3=9`, constrained dimer8,
rigid two-water12, and `K=f*T/2`. The one-step CSVR/NHC consumption test invokes
the actual controller outside runner setup with the independently declared
nine-mode dimension and compares runner velocities. Proper COM0/COM3 velocity
sampling, zero-charge/direct-Ewald controls, rigid/distance tangency, old COM3
rejection, state-only force-model switches and exact continuation are covered.
The alloy fixture parses through real InputConfig and constructs/runs a real
Simulation; MCP tests call actual tools and join actual worker lifecycles.

Intermediate harness failures remain archived and are not accepted gates:
water initially omitted required real units; a probe assertion initially used
an absent atoms field; invalid async input was initially asserted synchronous
rather than checking Failed; timestamp-only fixture directories collided under
concurrent tests and were corrected with an atomic sequence. The confirmed
covering logs below are the accepted final proof. No physical assertion was
weakened to address these harness issues.

## Accepted focused verification

All commands start with rtk. Full stdout/stderr remains in
[final-fix evidence](../final-fix/). All tests used default debug features.

| Exact command | Result | Log |
|---|---|---|
| `rtk proxy cargo test -p modena-core --test final_integration --test runner_hardening --test sampling_hardening --test sampling_state --test rigid_sampling --test force_hardening --test chemical_workflows --test io_and_lattice --test restart -- --nocapture` | 101 pass / 0 fail / 0 ignored | final-fix-core-covering-confirmed.log |
| `rtk proxy cargo test -p modena-core --lib lattice:: -- --nocapture` | 4 pass / 0 fail | final-fix-lattice-lib.log |
| `rtk proxy cargo test -p modena-core --test force_cache -- --nocapture` | 1 pass / 0 fail; warm serial/parallel mesh allocations0 | final-fix-force-cache.log |
| `rtk proxy cargo test -p modena-mcp --test final_integration --test chemical_workflows --lib -- --nocapture` | 23 pass / 0 fail / 0 ignored | final-fix-mcp-covering.log |
| `rtk proxy cargo clippy -p modena-core -p modena-mcp --all-targets -- -D warnings` | exit0 | final-fix-clippy.log |
| `rtk proxy cargo fmt --all --check` | exit0 | final-fix-fmt-check.log |
| `rtk proxy git diff 10bd9ee3077509ee472875f6697fd6d1a849753e HEAD --check -- . ":(exclude)docs/scientific-hardening/evidence/**"` | exit0 | source/current docs only |

Total:129 distinct passing tests, including12 new core and2 new MCP tests.
Full staged/committed diff whitespace checks report exit2 solely for preserved
raw test-output trailing whitespace and blank EOF lines. Those bytes remain
intact as evidence; the source/current-doc check excluding archived evidence
passes. This is not a production source whitespace failure.
The coordinator owns fresh full-workspace, feature/GPU/MPI Clippy, release,
Python, tooling/GPU-adapter, MPI2, numerical and example gates; they were not
repeated or claimed complete by this fix.

## Durable provenance and historical boundaries

The original [whole-branch review](final-review.md) and
[reviewer probe text/raw/manifests](../final-review-probe/) are copied unchanged.
All12 original probe file hashes were verified against the reviewer's manifest.
Its release rlib SHA256 remains recorded as
`794cca85de5b6eae1dd3ec9ef80387b4b5f905254991a758f8460765d107b295`; executable
probe binaries/rlibs were not committed. Original scratch probe files were not
modified. All218 pre-existing tracked scientific-evidence files were compared
byte-for-byte against BASE; [preserved-evidence hashes](../final-fix/final-fix-preserved-evidence.json)
record this check.

The historical67 numerical rows and24 successful matched-LJ processes remain
attached to their original commits/source hashes. Pair/FCC COM3 initialization,
FCC positions, pair force evaluation and timestep scopes used by that LJ
campaign are unchanged. The SPME static energy/force/virial formulas are
unchanged, so those archived static rows still describe that finite-mesh
Hamiltonian. Charged SPME initialization/temperature/controller interpretations
change; historical outputs are not relabeled as a fresh corrected trajectory
campaign. No fresh performance claim is made here.

## Self-review and remaining concerns

All three review findings are addressed with real RED/GREEN proof. Verified
that no mesh-force projection breaks the energy derivative contract; unsupported
restart metadata cannot reinstall COM3; constraint projection cannot bypass
new-model eligibility; zero fractions do not acquire populations; deterministic
ties cannot underflow label counts; schema-derived inventory reaches real tool
requests. Earlier constraint-rank/pressure, rigid/anisotropic NPT, converged
liquid/real ionic-water parity, optional production-backend and CPU memory/
performance limitations are unchanged. New force providers should explicitly
opt into conserved momentum only when justified. Scoped independent rereview
and coordinator final gates remain pending at this immutable handoff.

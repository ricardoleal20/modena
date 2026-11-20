# B07 pressure diagnostics: measured results

Measured on 2026-11-20 at revision `3345c6f09ce6a61f6c92c3273eadad6c35ab6f09`. Production validation remains **Diagnostic**. Aggregate scientific acceptance is **Failed**. Tensor availability, software-test success and correct classification of failed experiments do not establish a validated physical domain.

B07 adds conservative virial tensors, independent endpoint constraint reactions, stage-timed impulse pressure, separate BAOAB bath reactions, optional runner observations and stress-aware restart siblings. Legacy scalar observables and mechanical v2 checkpoints remain separate. See the [API migration](../releases/b07-pressure-api.md), [protocol](B07-PRESSURE-VALIDATION.md) and [design](../superpowers/specs/2026-11-06-b07-tensor-pressure-design.md). Production constrained/anisotropic NPT and GPU/MPI pressure validation are not established.

## Frozen measurements

| Family | Passed | Failed | Inconclusive | Pending |
| --- | ---: | ---: | ---: | ---: |
| Internal deterministic | 163 | 0 | 0 | 0 |
| Fresh conditional OU | 108 | 0 | 0 | 0 |
| External reference/provider | 175 | 24 | 46 | 135 |
| Final rotor, duration 1600 | 129 | 4 | 45 | 0 |
| Final water, duration 1600 | 123 | 3 | 52 | 0 |

The shared ledger has 1,145 rows, including execution and missing-action records: 765 Passed, 33 Failed, 143 Inconclusive and 204 Pending. The initial unlaunched clock denial remains separately recorded. Independent reviews approved evidence integrity and reporting after a separate attribution corrigendum. Their approval does not promote scientific status.

Both free-body NVE families pass their fixed four-rung native/MODENA continuum and state comparisons. Those results apply to the frozen zero-interaction rotor and nonprincipal-water fixtures and finite common windows. The real oblique Cu4 ASE 3.26.0/EMT family passes eight gates: all six virials are nonzero and the retained signed-strain derivatives agree. Its 51 request IDs comprise one info request and 50 evaluations: wire base plus 49 production snapshots. No DFT calculation, learned model or model-weight download was involved.

## Domain acceptance

Each outcome applies to the listed estimator and fixture on measurement revision `3345c6f`. VV denotes velocity Verlet. Snapshot rows cover conservative evaluations before integration. Grouped providers or controllers share the indicated outcome, subject to the configuration exclusions below. The production column describes API availability. All pressure observations retain Diagnostic status.

| Provider / backend | Integrator · constraints · controller | Estimator | Evidence or blocker | Outcome | Production availability |
| --- | --- | --- | --- | --- | --- |
| CPU known-empty contributions | snapshot · none · none | conservative virial | `stress_tensor` and `stress_bonded`, empty pair/bonded/composite identities | Passed | exact zero tensor |
| CPU pair (serial and Rayon) | snapshot · none · none | conservative virial | `stress_tensor`, harmonic and LJ periodic/scaled fixtures, 1 and 4 threads | Passed | six components |
| CPU Buckingham, Morse, soft sphere, tabulated and typed pair models outside the tested fixtures | snapshot · none · none | conservative virial | shared pair accumulator, separate fixture-specific strain evidence absent | Pending | six components |
| CPU bonded and complete composite | snapshot · unprojected coordinates · none | conservative virial | `stress_bonded`, bond/angle/dihedral six strains, connected images and composition | Passed | six components |
| CPU direct Ewald | snapshot · unprojected coordinates · none | conservative virial | `stress_electrostatics`, neutral/charged, scaled and excluded branches, fixed reciprocal modes | Passed | six components |
| LAMMPS CPU double, four pair and five interacting-water snapshots | run-zero · as registered · none | native reference convergence and table checks | `converge-pair-*`, `converge-interacting32*` and `outside-table`, pinned CPU double RIGID + KSPACE | Passed | external reference only |
| CPU direct Ewald, O–O probe | snapshot · none · none | energy, force and virial | `evaluate-pair-oo-r0-ewald` at 8/10/12, all nine companion gates pass | Passed | six components |
| CPU direct Ewald, charged pair and inside/outside excluded probes | snapshot · none · none | energy, force and virial | `evaluate-pair-{charge,excluded-in,excluded-out}-r2-ewald` at 8/10/12, nine tensor failures, U/F pass | Failed | six components |
| CPU LJ + direct Ewald, neutral interacting32 | snapshot · 32 rigid waters · none | energy, force and conservative virial | `evaluate-interacting32*-r2-ewald` at 8/10/12, fifteen energy failures, F/S pass | Failed | six components |
| CPU SPME, orthorhombic | snapshot · unprojected coordinates · none | discrete virial and trace | `stress_mesh` six strains at the reference cell and finite-grid energy oracle | Passed | six components |
| CPU SPME, charged/charged-scaled/inside-scaled | snapshot · unprojected coordinates · none | E/F/S and translation | `stress_mesh` 8/16/32 reference ladders, three fixture classifications | Passed | six components |
| CPU SPME, outside-excluded | snapshot · unprojected coordinates · none | E/F/S and translation | `stress_mesh` 8/16/32, maximum 32³ error and translation exceed 2e-4 | Failed | six components |
| CPU SPME, external pair probes and interacting32 | snapshot · none or rigid water as registered · none | matched energy, force and virial | `evaluate-*-spme`, required direct-Ewald matching did not pass | Pending | six components |
| CPU SPME, triclinic | any · any · any | conservative virial and pressure | production reciprocal cell guard, affine oracle is test-only | Unsupported | rejected |
| ASE 3.26.0 / EMT, periodic oblique Cu4 | snapshot · none · none | conservative virial | `ase-six-strain`, eight external gates, 49 production snapshots | Passed | six components |
| CPU zero / harmonic trap on the registered rotor | VV or BAOAB deterministic stages · distance bond · none | endpoint reaction and timed mechanical interval | `deterministic-gates` and `stress_vv`/`stress_baoab` analytic stage identities, trap tests use their specified forces | Passed | endpoint and interval |
| CPU zero, principal and nonprincipal water | VV · rigid water · none | endpoint and interval | `deterministic-gates`, two mass sets and refined independent Euler solutions for general rotations | Passed | endpoint and interval |
| CPU zero, free rotor and light nonprincipal water | VV · distance bond or rigid water · none | endpoint reaction, interval and state | `nve-rotor`, `nve-water-nonprincipal` and `evaluate-*-nve-*-zero`, four fixed common-time rungs | Passed | endpoint and interval |
| CPU zero, fixed water geometry | BAOAB O operator · rigid water · Langevin bath | projected OU law and impulse | formal conditional OU, 108 gates, 65,536 draws, operator only | Passed | separate bath ledger |
| CPU zero, free rotor gas | BAOAB · distance bonds · Langevin bath | endpoint, interval and ensemble law | `gas-rotor-1600-gates`, 129 Passed, 4 Failed, 45 Inconclusive | Failed | endpoint and interval |
| CPU zero, free water gas | BAOAB · rigid waters · Langevin bath | endpoint, interval and ensemble law | `gas-water-1600-gates`, 123 Passed, 3 Failed, 52 Inconclusive | Failed | endpoint and interval |
| LAMMPS CPU double, free rotor and water gas | native rigid/nvt · rigid bodies · native thermostat | native pressure and ensemble law | `preflight-gas-rotor` and `preflight-gas-water`, each has 23 inconclusive precision/ESS gates | Inconclusive | external reference only |
| CPU LJ + direct Ewald or SPME, interacting32 | VV · 32 rigid waters · none | endpoint, interval, state and native comparison | `nve-ewald` / `nve-spme`, matching direct-Ewald energy and charged tensors failed | Pending | endpoint and interval |
| CPU LJ + direct Ewald or SPME, interacting32 | BAOAB · 32 rigid waters · Langevin bath | endpoint, interval and native ensemble comparison | `nvt-ewald` / `nvt-spme`, matching failures and both inconclusive native thermostat preflights | Pending | endpoint and interval |
| CPU complete providers, including ASE within its restrictions | VV or BAOAB · none or supported independent bonds/waters · fixed cell, VV may add Berendsen or velocity rescaling | endpoint after velocity changes, interval before them | `stress_runner` timing checks, other trajectory/ensemble tuples lack all required physical gates | Pending | endpoint and interval for allowed combinations |
| CPU complete internal providers | VV or BAOAB · none · isotropic Berendsen or energy-based MC under its canonical-input guard | post-move endpoint and pre-move interval | `stress_runner` / `stress_restart` validate timing, NPT sampling outside B07 | Pending | endpoint refreshed, interval total missing after a box move |
| Legacy scalar models, GPU pair or MPI decomposition, composites containing an unavailable contribution | VV or BAOAB · any supported · any allowed | model tensor and total pressure | ForceModel scalar fallback, no independently implemented complete tensor for these domains | Unsupported | scalar mechanics available, total tensor missing |
| Internal CPU or eligible external provider | Nosé–Hoover chain or custom integrator without step stress · any supported · any allowed | step endpoint and interval | integrator stress hook absent, NHC lacks a derived thermostat ledger | Unsupported | step endpoint missing, interval absent |
| Any provider | any · periodic winding, ambiguous image or singular graph · any | constraint endpoint and interval | `constraint_geometry` image, closure and rank checks | Unsupported | rejected |
| Internal CPU | any · nonempty · any box rescaling | box control using constrained pressure | runner configuration rejects constrained box moves | Unsupported | rejected |
| Internal CPU | any · none · anisotropic Berendsen or legacy Parrinello–Rahman piston | tensor-driven box control | production configuration guard, direct orthorhombic diagonal update has engineering coverage only | Unsupported | experimental direct API, production rejected |
| ASE EMT | any · any explicit constraints · any | model and pressure | ASE total-provider input rejects added topology, electrostatics, constraints and barostats | Unsupported | rejected |
| ASE EMT | any · none · any barostat | box control | ASE total-provider input rejects barostats | Unsupported | rejected |
| ASE EMT, missing/invalid provider stress | snapshot · none · none | conservative virial | required six-component, finite, symmetric provider property, `external_provider` contract | Unsupported | evaluation fails |
| CPU complete providers | VV or BAOAB at a full v2 continuation seam · any supported · any allowed | initial endpoint and first interval | `stress_restart`, mechanical v2 lacks the previous conservative tensor | Unsupported | missing legacy contribution until ordinary evaluations refill history |

## Acceptance blockers

The 24 external failures comprise two distinct groups:

- **Nine charged-pair conservative-tensor comparisons** fail at direct cutoffs 8/10/12. Their U/F companions pass. A post-observation source/raw comparison supports a background-convention explanation limited to these nine rows: MODENA includes `U_bg I`, where `U_bg = -pi*Q^2/(2*kappa^2*V)`, while the pinned native global-virial loop lacks that explicit diagonal term. For Q = -.4, kappa = .8 and V = 12^3, the term is -0.00022725641302009494; retained differences agree to residual <= 1.01e-8. This is no preregistered harmonization or corrected pass. The original gates remain Failed. See the [pinned native source](https://github.com/lammps/lammps/blob/9c5ab448c78a14fd534619622162ba418d6a1fb1/src/KSPACE/ewald.cpp#L452-L470).
- **Fifteen neutral interacting-water energy comparisons** fail. Scaled errors 2.045680175280889e-6 at cutoff 8 and 2.0456911230142135e-6 at 10/12 exceed 2e-6. Native finest U is -0.09897314054281026; MODENA U is -0.09897089238321222, a retained difference of 2.2481595980394786e-6. Their F/S companions pass. This discrepancy is unresolved and is not explained by the charged-background observation.

Failed matching prerequisites block interacting Ewald/SPME dynamics. Both native rigid/nvt free-gas preflights completed but retain 23 Inconclusive rows each because precision/ESS requirements are missing; interacting NVT is blocked. Exactly 12 of 24 registered native jobs ran, with no replacements. Every blocked case has an explicit prerequisite link and supplies no physical acceptance.

The independent mesh regression retains the historical outside-excluded 32^3 error 2.4207810466960128e-4 and translation change 4.8298199772208068e-4 against 2e-4. The final four-fixture 8/16/32 physical ladder was executed once. All 12 libtests pass the engineering harness, while outside-excluded and aggregate scientific accuracy remain Failed.

## Gas inference and continuation

| Species | Duration | Passed | Failed | Inconclusive |
| --- | ---: | ---: | ---: | ---: |
| Rotor | 400 | 0 | 0 | 178 |
| Rotor | 800 | 103 | 4 | 71 |
| Rotor | 1600 | 129 | 4 | 45 |
| Water | 400 | 42 | 0 | 136 |
| Water | 800 | 116 | 2 | 60 |
| Water | 1600 | 123 | 3 | 52 |

All 1,024,000 sample rows and 192 continuation records were checked. The 400 -> 800 -> 1600 extensions continue the same states and RNG streams without repeating burn. At h = .04, the nominal .1 cadence uses actual alternating .12/.08 step gaps. Earlier files and failed rows remain retained.

The final seven failures concern rotor endpoint/interval zz drop-coarse shifts and finest xy means, water endpoint/interval xz intercepts and rotation CDF50 intercept. Ninety-seven final criteria remain Inconclusive. All 60 weak-order fits are unresolved; these data support no second-order BAOAB claim.

Endpoint scalar h^2 intercepts show limited agreement, rather than all-component acceptance:

| Species | Intercept | SE | 95% interval |
| --- | ---: | ---: | --- |
| Rotor | 0.018519858493223313 | 1.6437887295501096e-5 | [0.018480989066289, 0.018558727920157628] |
| Water | 0.01851308130195825 | 1.8821448446380922e-5 | [0.0184685756485228, 0.018557586955393703] |

Pressure rows have at least 70 blocks per seed and ESS 4453.77 for rotor, and 78 blocks and ESS 4974.91 for water. Worst 95% halfwidths 0.0002137955241 and 0.0001763008689 nevertheless exceed the 0.00009259259259 target. Some machine display coordinates are clipped to the ESS/block targets. Unrounded means, ESS, blocks and actual CI halfwidths remain in the stats sidecars and determine classification.

Post-measurement inverse-square-root forecasts, assuming unchanged variance/mixing, give work multipliers 5.33144 and 3.62540: about 665.3 million and 403.3 million additional molecule-steps. Pooled measured ensemble time implies approximately 5389.1 and 3266.5 extra seconds. These exceed the campaign limits and were not executed. Species-specific marginal runtime and the work needed to resolve weak order or remove failures remain unknown. More samples do not guarantee that a failed criterion will pass.

## Observed cost

The common driver uses dt = .005 and identical initial states. Final old/new positions and new legacy/reporting q/v agree exactly. A single short batch per row provides no replicated uncertainty or universal performance estimate.

| Batch | Old legacy ns | New legacy ns | New/old | Reporting ns | Reporting/new |
| --- | ---: | ---: | ---: | ---: | ---: |
| Rotor, 128 steps | 176083 | 591708 | 3.36039 | 669541 | 1.13154 |
| 32 waters, 32 steps | 244750 | 2102583 | 8.59074 | 4524375 | 2.15182 |

The release accepts these measured constrained costs for Diagnostic use. The later empty-constraint optimization has no measured bearing on either ratio. No speedup is claimed. Eight alternating SPME observer pairs after one warmup give total ratio 1.00326, with pair ratios .98849–1.02313. The 32^3 cached-factor delta is 1,572,864 bytes. Endpoint reaction including validation costs 916 ns/26,959 ns. Rotor v2 -> v3 checkpoints grow 771 -> 3,545 bytes; water 16,050 -> 63,980 bytes. Stage entries, rollback/reconstruction sizes and exact output bytes remain in `costs.json`.

## Resources and clock

The formal campaign executed 65 sequential children, all reaped and reconciled, with no unknown work suffix, interruption or budget stop. Total primary integration work is 330,258,705 molecule-steps, including 316,800,000 gas steps and 13,443,000 native steps. The deterministic RK4 references add 153,600 reconstructed oracle substeps, separate from the MD counter and included in the deterministic child cost. Fresh formal conditional draws are 65,536; the 196,612 historical-plus-current charged operator calls do not enlarge the fresh statistical sample size.

Child-wall sum is 2,593.183403249 seconds; CPU is 2,548.203794 user + 32.273582 system seconds. Maximum child/worker RSS is 74,825,728 bytes; ASE sampled tree is 79,380,480 bytes. Twelve native jobs use 9.132055791 seconds. Ensemble wall is 2,566.1358605 seconds with peak 65,519,616 bytes. Whole-controller CPU/RSS was not independently measured. Compiler memory is separate from numerical child limits.

The original window remains 2026-11-20 18:59:03.599210 -> 20:59:03.599182 UTC. Python 3.9/macOS process-relative monotonic time caused an initial prelaunch denial. A reviewed append-only Mach-continuous mapping fixed shared deadlines without changing that window. The legacy boot/monotonic prefix is unknown. The 6,031.144813-second report duration is explicitly UTC-derived across it; 4,462.049627542 seconds is shared elapsed from the frozen reference. The last physical child ended at 20:37:27.995339 UTC. Later archival and metadata corrections perform no new physical work.

## Provenance and local archive

| Identity | SHA256 or revision |
| --- | --- |
| Measurement source | 3345c6f09ce6a61f6c92c3273eadad6c35ab6f09 |
| Internal protocol | 89abb8b7a45b5d443bb886fdba0642e78b62fc7f64fe6c8438a06d2886b0759d |
| External epoch | 97d45c8b934c2999dc8f8d77cfcfe563d7bd2a9bc7f1a8a508b9944ddf43381e |
| Rust example | b528493e993d3756f1401daefec6453fa92a26c9e767c7983652fe22a3757eeb |
| Full archive inventory | 7d931afae6961329ecfa53948e9991b5704af45ec8a51f4725b0f5af094d948e |
| Gate ledger | e35d173c6b55fe1d9f3a012715c1554b460e5f4c4b1b3d05b6e9d43590908d56 |
| Corrected readable summary | e324b4b48faeb8121c1082fef28cb7033a16ff7c5044f937f99c120c153a513f |

The build record’s `features: []` lists explicit Cargo CLI feature flags. Cargo enabled the empty default feature set and no optional backend. The frozen summary did not separately inventory resolved artifact features.

The full epoch is retained locally at `target/b07-evidence/combined-campaign-20261120/epoch-v6`: 1,470 files and 847,609,834 bytes. Budget, original-denial and correction metadata are alongside it. Independent reviewers verified 1,475 archive hashes and the relevant receipts, family seals and origins. Registered raw epochs total 861,511,947 logical bytes, below 5 GiB. Original evidence and classifications are immutable; the corrected summary is a separate metadata version. [Compact results](B07-RESULTS.json) are reporting data, not a production input or a validated-domain credential.

## Software verification

**The original 600-second window at the measured revision ended Failed and incomplete.** Format, all-target Clippy, doctests and compile-only GPU/MSRV 1.90/MPI checks passed. Completed Rust binaries reported 513 passed, one failed and three ignored. An interrupted `spme_calibration` retained 57 additional observed passes and one ignored case; its active case and remaining suffix were incomplete. `thermostat_statistics` and `water` were not run after the single 600-second deadline. Python package: 26 passed/nine skipped; bench: 86 passed/one failed/two skipped.

The MCP nonfinite assertion and the benchmark's release-binary location failed in that window. The full mesh regression completed once at the measured revision as described above.

**The separate post-fix software window passed its scoped checks at `a0012846318be454d0be4003806a08229c0602a6`.** The one-worker run took 485.29 of 600 seconds. Cargo compiled 60 workspace test executables. The 26 selected binaries reported 330 passed, 0 failed and 3 ignored. Four passing optional-provider cases explicitly skipped their live ASE paths. Calibration completed with 102 passed and three ignored, thermostat statistics with three passed, and water with four passed. MCP passed all 21, including the actionable kinetic-pressure diagnostic at failed step 1. All 11 restart tests passed, including the zero-constraint VV checkpoint and constrained OU corruption controls.

Python 3.12 package discovery reported 26 passed and 9 skipped. Bench discovery reported 87 passed and 2 skipped. The real benchmark contract verified the path and SHA256 of the isolated release executable, numeric drift, and absence of XYZ/CSV output. Both Python 3.9 and 3.12 passed 35 adapter controls and 13 Darwin clock controls under isolated startup. The package declares Python ≥3.11. The full package/mock suite ran on 3.12 only.

Doctests reported one pass and one ignored. Separate bounded checks covered formatting, all-target warnings-denied Clippy for default/GPU/MPI, Rust 1.90 default/GPU and the standalone validation harness, release workspace and benchmark compilation, and GPU/MPI test compilation. Hardware and Windows execution were outside this local run. The scoped run omitted 34 other corrected-source Rust binaries, including the full mesh and fixed-geometry OU tests. No earlier binary supplied a corrected-source pass. The compact results list the unexecuted artifacts and command scopes.

The measured source closure changed in constraint stress, state checks, restart validation, pressure assembly and two test files. All physical observations and classifications remain attached to `3345c6f`. The corrected tree has engineering evidence and continues to emit Diagnostic pressure.

# MODENA Scientific Hardening Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development to implement this plan continuously.

**Goal:** Finish the accepted audit's scientific correctness, reproducibility, chemical workflow, validation and measured performance priorities for MODENA's classical supported domain.
**Architecture:** Keep the ForceModel/Integrator seams; add explicit state serialization, physical DOF metadata, weighted pair policies and runner guards. Use canonical isotropic Monte Carlo volume proposals rather than claiming the current piston provides validated NPT. Keep acceleration extensions honest and testable.
**Tech Stack:** Rust, serde JSON/TOML, Rayon/rustfft, standard-library Python, existing GPU/MPI optional backends.
**Spec:** docs/superpowers/specs/2026-10-03-scientific-hardening.md

## Global Constraints

- All shell commands start with `rtk`; use `rtk proxy` for unsupported commands.
- Worktree: /Users/ricardoleal20/.codex/worktrees/modena-scientific-hardening/modena. Base b68abd112339aa6970261f470a572da6633bc617. Preserve user files and audit evidence in the original checkout.
- You are not alone in the codebase; do not revert other changes. Implementers must not spawn agents. Task implementations run sequentially, followed by independent spec and quality review.
- Follow test-driven-development: observe meaningful failing regression before its production fix, then run covering checks. Tests use independent expectations and real components.
- Unsupported scientific combinations return clear input errors; never silently produce unvalidated results. Documentation distinguishes supported, experimental, and absent capabilities.
- Checkpoint state methods on Integrator, Thermostat and Barostat: `fn checkpoint_state(&self) -> crate::error::Result<serde_json::Value>` and `fn restore_state(&mut self, state: &serde_json::Value) -> crate::error::Result<()>`. Stateless defaults encode null and reject non-null. Stateful implementations include their kind and all RNG/extended state; preserve exact float round trips.
- `ForceModel::neighbor_cutoff(&self) -> Option<f64>` default None; concrete/composite models return their maximum required cutoff. Replica exchange rejects absent/invalid cutoff rather than guessing.
- No external publication. Commits use gitmoji + Action format; no AI coauthor.

## Task 1: Correct sampling, constraints, DOF and replica state

**Ownership:** core thermostat.rs, integrator.rs, constraints.rs, system.rs, nose_hoover.rs, replica_exchange.rs; sampling/constraint tests; Cargo.toml, core Cargo.toml/Cargo.lock for state dependencies; ForceModel trait cutoff method only. No runner/config edits.

- [x] Read PHYSICS-NOTES.md and probes.rs in original docs/audits/2026-10-03-modena-vs-lammps; reproduce the canonical CSVR mean/variance and free-atom constrained BAOAB variance failures as tests.
- [x] Fix CSVR chi-square scaling (2*Gamma), canonical signed rescaling when necessary, analytic one-step moments and stationary moments. Handle low DOF safely.
- [x] Apply constrained BAOAB projection per drift/OU stage; apply velocity corrections only from each position projection. Do not replace stochastic velocities with the full displacement. Free atom variance must be T/m with or without an unrelated constrained dimer.
- [x] SHAKE/RATTLE/SETTLE errors propagate when residuals remain above tolerance, indices/masses are invalid, or solves fail. Update callers/tests as required for Result, without ignoring failure.
- [x] Add explicit constraint count and conserved COM DOF metadata to System, exposed setters/accessors; constructor defaults preserve 3N for manual systems. Integrators set constraints metadata. Langevin retains COM stochastic DOF. Runner wiring follows in Task 3. Fix triclinic nearest image using a bounded mathematically correct nearest-lattice search and scaling cell basis consistently; add the audited skew-cell regression and neighbor correctness handoff to Task 2.
- [x] Provide checkpoint_state/restore_state on Integrator and Thermostat, using serializable stable RNG state (ChaCha12 with serde or equivalently exact state), and NHC extended variables. Add core serde_json float_roundtrip dependency. Restore validates kind/parameters; stochastic seam tests on direct integrator/thermostat state.
- [x] Fix replica exchange full mechanical state and cached forces, velocity temperature factors, acceptance sign, final PE tracking, validated temperatures/dt and cutoff. Add force/cached-state regression. Add ForceModel cutoff seam; Task 2 implements providers.
- [x] Assess constraint impulse scalar virial using an independently understood rotating rigid dimer test. Implement if correct; otherwise write an explicit limitation for Task 3 to surface and guard constrained NPT (never invent a stress).
- [x] Run core covering tests and format owned files. Commit cohesive changes and write task-1-report.md with RED/GREEN evidence, public interfaces, remaining limits and test commands/results.

## Task 2: Correct force/energy/pressure and remove measured allocation costs

**Ownership:** force_model.rs (remaining implementation), neighbor.rs, potential/tabulated.rs, ewald.rs, spme.rs, composite.rs, decomp.rs, gpu_comm.rs if contract changes; force/pair/electrostatic tests and core lib exports as needed. Do not edit Task 1 files except integration callsites strictly necessary and report these.

- [x] Reproduce tabulated energy gradient and direct Ewald isotropic strain virial failures. Joint Hermite interpolation uses force-derived endpoint slopes; shifted-force cutoff remains consistent. Include tabulated in existing derivative properties.
- [x] Correct reciprocal direct Ewald strain virial and net-charge background virial (W_bg=3 U_bg), with independent finite strain derivative tests. Check SPME background and weighted exclusion corrections.
- [x] Replace linear exclusion scans with an indexed representation. Add `PairScale { i: usize, j: usize, lj: f64, coulomb: f64 }` and NeighborList `with_pair_scales(cutoff, skin, scales)` plus `pair_scale(i,j)`; retain existing constructors as binary 0/0 exclusion compatibility. Do not drop a pair when only one interaction is zero. Apply independent LJ/Coulomb scales, including reciprocal erf corrections proportional to (1-coulomb scale). Electrostatic constructors `with_pair_scales` match this policy.
- [x] Implement neighbor_cutoff for pair, electrostatic, composite and backend models. Triclinic neighbor search must fall back to guaranteed correct pair enumeration if bounding-box cell stencil is unsafe; test the audited skew pair.
- [x] Fix DomainForceModel accumulation contract with an existing nonzero-force fixture. Optional GPU/MPI code must remain compileable; unsupported fractional scaling backends reject it explicitly.
- [x] Cache SPME FFT plans, grids and scratch per instance using a thread-safe owned workspace; reuse buffers and invalidate influence data by box/parameter key. Eliminate global cross-simulation cache thrashing. Cache pair force scratch if profiling/reproducible microbench supports it, retaining numerical tolerances and deterministic serial path.
- [x] Run covering force properties, Ewald/SPME references, exclusion policy tests; provide allocation/reuse or timing evidence. Commit and report RED/GREEN, interfaces and restrictions.

## Task 3: Integrate validated ensembles and complete checkpoints

**Ownership:** runner.rs, config/mod.rs, config/builder.rs, io/input.rs, io/restart.rs, barostat.rs, runner/restart/config tests; CLI check path and regression tests; lib exports and needed serde derives. Use Tasks 1/2 interfaces.

- [x] Reproduce and fix stale force/PE/virial after every accepted box change; recompute caches before observable/output/next kick. Wire physical DOF before velocity initialization/thermostats; project initial constraints and velocities and rescale using correct DOF. Conserved COM: subtract 3 only for dynamics preserving the initialized zero momentum; Langevin includes stochastic COM.
- [x] Surface the assessed constraint-stress limitation in Simulation/RunSummary metadata and CLI result output: constrained pressure is a diagnostic omitting constraint impulse stress. MCP consumes the same metadata in Task 4. Guard constrained NPT; do not mark this observable validated.
- [x] Validate finite positive dt, tau/masses/box, gamma>=0, T>=0, topology/constraint indices, lengths, sampling cadences, electrostatic settings, incompatible thermostats, unknown keys and duplicate contradictory pairs before execution. Explicitly reject anisotropic scalar-virial coupling, constrained box rescaling and unvalidated legacy piston through production config with explanatory errors.
- [x] CLI `check` must not truncate/create trajectory, log or restart outputs: it currently instantiates output writers. Add a real CLI sentinel-file regression and construct validation with outputs disabled, while preserving all scientific setup/identity checks.
- [x] Implement isotropic Monte Carlo barostat with logV symmetric proposals, log acceptance -beta*(delta U+P delta V)+(N+1)*log(Vnew/Vold); seed, interval, log-volume width, counters and RNG checkpoint state. Unconstrained particles only. Trial state rollback includes box, positions, forces and neighbor invalidation. Integrate `apply_with_model` Result trait hook and canonical T consistency validation; tests include ideal gas expected gamma volume mean/variance and rejected/accepted cache integrity.
- [x] Version 2 atomic checkpoint: exact mechanical state plus masses/species/charges/elements, physical signature, DOF and Integrator/Thermostat/Barostat state. Write sibling temp file, sync and atomic rename; reject malformed/nonfinite/overflow/trailing corrupt data. Preserve legacy v1 parser.
- [x] Default restart continuation validates physical identity (model parameters, atom identity, topology, dt, ensemble, seeds/state); execution/output settings may change. Explicit `restart_mode = "state_only"` permits intentional NVT->NVE protocol switch; v1 requires explicit state_only rather than pretending full continuation. Update existing intentional protocol-switch fixtures.
- [x] Tests compare continuous vs split NVE, Langevin, CSVR, NHC, MC NPT and constraints; mismatch rejection tests. Run full workspace suite to catch integration regressions; commit and report.

## Task 4: Deliver typed chemical workflows and practical Python/MCP access

Accepted scope extension (ledger ruling, 2026-10-03): correct the existing real charge conversion to `q*=q*sqrt(C/(sigma*epsilon))`; `sqrt(C*sigma/epsilon)` previously multiplied physical Coulomb energy by sigma². Ownership includes units.rs and independent physical energy/force/reference-scale regression tests. Old real electrostatic trajectories at sigma!=1 use a different Hamiltonian; normal checkpoint signatures include charge factor/charges and reject incompatible continuation. Reduced and sigma1 cases are unaffected. Historical audit data remain preserved and no migration is promised.

**Ownership:** io/xyz.rs, config/builder.rs/mod.rs and runner molecular construction, python/modena, mcp crate, new example inputs and workflow tests. Preserve Task 3 validation and restart wiring.

- [x] Add explicit typed atom import (documented extended XYZ Properties fields or equally unambiguous atom records) carrying per-atom type, mass and charge; plain XYZ legacy behavior remains explicit, typed malformed/unknown species rejected. Resolve paths relative to config file, including Python generated config behavior.
- [x] Graph-derived special pairs use shortest bond connectivity distances 1/2/3 for 1-2/1-3/1-4. Add independently configured LJ and Coulomb scales (default legacy exclusions retained), compatible explicit exceptions with clear precedence; tests on a four-atom chain and long-range Coulomb ensure no double scaling.
- [x] Python typed config supports real/reduced units, electrostatics, topology/scales, constraints, NHC/CSVR, barostat MC, water, restart modes; robust TOML quoting, [[model.types]], absolute source paths, subprocess error propagation and usable numeric outputs. Standard-library tests acceptable; do not require installation just for tests.
- [x] Expand MCP discovery schema to actual supported config, include capability/limitations matrix and backend presence. Correct NVE classifier to account for barostat and transition mean/averaging behavior; behavior tests on real request paths.
- [x] Ship runnable examples for typed molecular NVT, ionic SPME, MC NPT and continuation restart, with units/parameter convention explicit; fixed-volume water example remains supported. Run Rust/Python and end-to-end check/run workflows. Commit/report.

## Task 5: Reproducible validation/benchmark tools and honest release evidence

**Ownership:** bench scripts/decks, criterion bench files, validation harness/tools, README/docs and scientific hardening results. No speculative production rewrites.

- [x] Replace failure-swallowing benchmark behavior: actual configured steps, stdout/stderr logs, return code, thread count/precision/potential/seed/hardware metadata, raw wall/CPU/RSS with unavailable values null, drift, declared warmup, repetitions, same-resource comparison. Fix numeric renderer. Keep engine init and no-output scope explicit. Add small executable fixture tests for failed run and successful metadata generation.
- [x] Criterion separates steady force/integration from neighbor rebuild costs. Add molecular/ionic SPME workloads and scaling ladders; every performance claim cites conditions. No test/build concurrency during final timing.
- [x] Scientific validation runner exports machine-readable multiple-seed canonical mean/variance, constrained/unconstrained sampling, NVE dt convergence, force finite differences, electrostatic grid/strain consistency, restart seams; reference independent analytic or LAMMPS fixtures. Long campaigns optional documented commands; CI smoke is meaningful and bounded.
- [x] Correct README and validation claims: robust supported pairs/ions/fixed-volume molecules vs experimental acceleration/unsupported full anisotropic/rigid NPT/general chemical completeness. Preserve original audit as historical evidence; write docs/scientific-hardening/RESULTS.md with implemented inventory, checks, measured timings and remaining competition gaps.
- [x] Run full workspace and Python checks, release build, optional GPU actual adapter tests and MPI 2-rank driver when available; run original physical probes against new default artifact and new matched benchmarks. Reuse original baseline evidence rather than timing unrelated old binaries.
- [x] Commit documentation and tools; report exact test counts/results and timing artifacts. Independent final whole-branch review follows, fix load-bearing findings, deliver branch and local result links without publishing.

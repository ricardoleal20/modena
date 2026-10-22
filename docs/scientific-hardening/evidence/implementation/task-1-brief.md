## Task 1: Correct sampling, constraints, DOF and replica state

**Ownership:** core thermostat.rs, integrator.rs, constraints.rs, system.rs, nose_hoover.rs, replica_exchange.rs; sampling/constraint tests; Cargo.toml, core Cargo.toml/Cargo.lock for state dependencies; ForceModel trait cutoff method only. No runner/config edits.

- [ ] Read PHYSICS-NOTES.md and probes.rs in original docs/audits/2026-10-03-modena-vs-lammps; reproduce the canonical CSVR mean/variance and free-atom constrained BAOAB variance failures as tests.
- [ ] Fix CSVR chi-square scaling (2*Gamma), canonical signed rescaling when necessary, analytic one-step moments and stationary moments. Handle low DOF safely.
- [ ] Apply constrained BAOAB projection per drift/OU stage; apply velocity corrections only from each position projection. Do not replace stochastic velocities with the full displacement. Free atom variance must be T/m with or without an unrelated constrained dimer.
- [ ] SHAKE/RATTLE/SETTLE errors propagate when residuals remain above tolerance, indices/masses are invalid, or solves fail. Update callers/tests as required for Result, without ignoring failure.
- [ ] Add explicit constraint count and conserved COM DOF metadata to System, exposed setters/accessors; constructor defaults preserve 3N for manual systems. Integrators set constraints metadata. Langevin retains COM stochastic DOF. Runner wiring follows in Task 3. Fix triclinic nearest image using a bounded mathematically correct nearest-lattice search and scaling cell basis consistently; add the audited skew-cell regression and neighbor correctness handoff to Task 2.
- [ ] Provide checkpoint_state/restore_state on Integrator and Thermostat, using serializable stable RNG state (ChaCha12 with serde or equivalently exact state), and NHC extended variables. Add core serde_json float_roundtrip dependency. Restore validates kind/parameters; stochastic seam tests on direct integrator/thermostat state.
- [ ] Fix replica exchange full mechanical state and cached forces, velocity temperature factors, acceptance sign, final PE tracking, validated temperatures/dt and cutoff. Add force/cached-state regression. Add ForceModel cutoff seam; Task 2 implements providers.
- [ ] Assess constraint impulse scalar virial using an independently understood rotating rigid dimer test. Implement if correct; otherwise write an explicit limitation for Task 3 to surface and guard constrained NPT (never invent a stress).
- [ ] Run core covering tests and format owned files. Commit cohesive changes and write task-1-report.md with RED/GREEN evidence, public interfaces, remaining limits and test commands/results.


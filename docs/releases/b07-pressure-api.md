# Tensor pressure API migration

## Reading pressure

Velocity Verlet and Langevin BAOAB report two diagnostic observations. The endpoint uses the current kinetic tensor, conservative model virial and independently solved constraint reaction after velocity and box controllers. The interval uses the kinetic, model and mechanical constraint impulses from the preceding integration step before those controllers. A velocity-only controller leaves that interval intact. A box move gives it a missing total with the reason `post_step_transform` and retains its original volume and timestep.

`Simulation::pressure_observation()` and `Simulation::interval_pressure()` borrow the latest observations in reduced units at every completed step. A fresh run has an endpoint and no preceding interval. A scalar-only force model retains its energy and scalar virial and supplies a missing model tensor with provenance. A step from an integrator without stress support, including the Nosé–Hoover chain, reports a missing endpoint total and no interval. The new scalar is absent when the total tensor is unavailable. Failed simulations expose neither getter.

Set `pressure_jsonl = "results/pressure.jsonl"` under `[run]` to write the initial snapshot and each `log_every` sample. Each JSON line carries the global `step`, display `time`, `endpoint` and optional `interval`. The matching `RunEvent::Pressure` carries the same values. Order is xx, yy, zz, xy, xz, yz. In real units, kinetic, model and constraint components (including the interval averages) are in kcal/mol, total and scalar pressure are in atm, volume is in Å³, and time and interval duration are in ps. A missing scalar or absent interval serializes as JSON null. A missing tensor serializes its reason list. Reading a serialized observation confers diagnostic status.

Use a distinct pressure output path from the log, trajectory, restore and restart files. The runner checks file identities before truncation. A failed check can leave an empty file at a previously absent pressure path and any newly created output directories.

The six-field `Energies` record and existing CSV columns retain their established meanings. Legacy `pressure` is `(2 KE + model virial) / (3 V)` and excludes constraint reaction and impulse stress.

## Updating Rust callers

- Add `pressure_jsonl: None` to exhaustive `RunSpec` literals, or use `..RunSpec::default()`. Omitted TOML and JSON fields default to `None`.
- Handle `RunEvent::Pressure { step, time, endpoint, interval }` in exhaustive event matches. A consumer interested only in existing events can use `RunEvent::Pressure { .. } => {}`.
- Implement `Barostat::apply` with return type `modena_core::Result<()>` and return `Ok(())` on success. Callers propagate or handle the result of `barostat.apply(system, &virial, dt)`. The compatibility `apply_with_model` method propagates errors. The full-cache hook `evaluate_move` returns the accepted cell's force-model evaluation, or `None` when the original cell and its caches apply.

## Controller and geometry limits

Scalar calls to anisotropic Berendsen and legacy piston controllers return a configuration error before changing any state. For direct experiments in an unconstrained orthorhombic cell, `anisotropic().apply_with_pressure(system, observation, dt)` uses the three diagonal components of a complete finite endpoint in reduced units at the same positive volume. A caller performing a direct move must refresh its own forces and model evaluation on the updated cell. Production runner configurations continue to reject anisotropic coupling, constrained box rescaling and the legacy piston. Tensor availability alone grants no sampling validation.

A periodic constraint with equally near bonded images, including an exact half-cell tie, now fails explicitly. Choose initial positions and a cell that define a unique bonded image. The v2 mechanical checkpoint keeps energy, scalar virial, cached forces and controller state. Its tensor is explicitly missing on continuation until the first ordinary endpoint evaluation, and the first interval carries that missing initial contribution.

## Continuing from checkpoints

Runner checkpoints now use `MODENA_RESTART 3`. Use `modena_core::io::restart::load_with_stress` to read them. It returns a `StressCheckpoint` with the v2 mechanical payload under `mechanical` and optional pressure state under `stress`. For a v2 export readable by an older binary, call `save(path, &checkpoint.mechanical)`. The existing `save` and `load` APIs retain the mechanical format. The new reader also accepts v1 and v2 files, returning `stress: None`.

Full v3 continuation preserves the endpoint after controllers, the interval before controllers, each projection ledger, OU bath exchange, and the cached model evaluation. Estimator version 1 binds these timings and the constraint solver options to the mechanical physical signature. Unavailable history at initialization, after a state-only import, or at an immediately re-saved legacy seam has no interval, stages or bath. State-only imports start a new protocol. Reading a checkpoint gives its pressure observations diagnostic status.

On full v2 continuation the endpoint reports `legacy_restart` and the seam has no recorded interval. After one ordinary measured step the interval still lacks its initial model contribution. The next measured step can supply a complete interval when the model supplies both endpoints. This preserves the conservative force cache and RNG sequence, including for a live external provider. Supported isotropic scalar coupling and energy-based Monte Carlo use their existing mechanical inputs. Direct tensor controllers require a complete endpoint.

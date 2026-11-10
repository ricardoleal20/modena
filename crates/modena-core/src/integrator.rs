//! Integrators.
//!
//! An `Integrator` owns the time-stepping scheme and drives force evaluation
//! through a `ForceModel`. Phase 1 ships:
//! - `VelocityVerlet` — symplectic, NVE.
//! - `LangevinBAOAB`  — stochastic, NVT (Leimkuhler–Matthews BAOAB).
//!
//! Convention: `system.forces` holds F(t) on entry and F(t+dt) on exit, so the
//! runner computes forces once before the loop and the first half-kick uses the
//! cached F(t). This avoids a redundant force evaluation per step.
//!
//! M7: the kick/drift loops operate on the flat `3N` SoA buffers directly.

use chacha20::ChaCha12Rng;
use rand::SeedableRng;
use rand_distr::{Distribution, StandardNormal};

use crate::constraint_stress::{constraint_geometry, ProjectionReport, ProjectionStage};
use crate::constraints::Constraints;
use crate::error::Result;
use crate::force_model::{ForceModel, Virial};
use crate::neighbor::NeighborList;
use crate::stress::{
    endpoint_pressure, kinetic_tensor, pressure_total, ConstraintSolverOptions, IntervalPressure,
    ModelEvaluation, StepEvaluation, StepStress, SymTensor6, TensorAvailability, ValidationStatus,
};
use crate::system::System;
use crate::thermo::Energies;

/// SHAKE/RATTLE convergence: residual tolerance (length² for SHAKE, velocity·
/// length for RATTLE) and the iteration cap. Tight defaults — the constrained
/// molecule test verifies energy conservation on the constraint manifold.
const SHAKE_TOL: f64 = ConstraintSolverOptions::DEFAULT.tolerance;
const SHAKE_MAX_ITER: usize = ConstraintSolverOptions::DEFAULT.max_iterations;

/// Time-stepping scheme. `step` advances the system by `dt` and returns the
/// post-step thermodynamic observables.
pub trait Integrator: Send {
    /// Stateless integrators encode null; stateful implementations override.
    fn checkpoint_state(&self) -> Result<serde_json::Value> {
        Ok(serde_json::Value::Null)
    }
    fn restore_state(&mut self, state: &serde_json::Value) -> Result<()> {
        crate::thermostat::restore_stateless(state)
    }
    fn step(
        &mut self,
        system: &mut System,
        model: &dyn ForceModel,
        neighbors: &mut NeighborList,
        dt: f64,
    ) -> Result<Energies>;

    /// Advance while reporting mechanical stress when the scheme supports it.
    /// The caller supplies an evaluation for the cached input forces, positions and cell.
    fn step_with_stress(
        &mut self,
        system: &mut System,
        model: &dyn ForceModel,
        neighbors: &mut NeighborList,
        dt: f64,
        _initial: &ModelEvaluation,
    ) -> Result<StepEvaluation> {
        Ok(StepEvaluation {
            energies: self.step(system, model, neighbors, dt)?,
            stress: None,
        })
    }

    /// External integrators can retain the default solver and legacy trajectory.
    fn step_with_stress_options(
        &mut self,
        system: &mut System,
        model: &dyn ForceModel,
        neighbors: &mut NeighborList,
        dt: f64,
        initial: &ModelEvaluation,
        solver: ConstraintSolverOptions,
    ) -> Result<StepEvaluation> {
        if solver != ConstraintSolverOptions::DEFAULT {
            return Err(crate::ModenaError::Config(
                "this integrator does not support custom constraint solver options".into(),
            ));
        }
        self.step_with_stress(system, model, neighbors, dt, initial)
    }
}

/// Velocity-Verlet (NVE). Expects F(t) already in `system.forces`.
///
/// With `constraints` set, applies SHAKE to the post-drift positions and RATTLE
/// to the post-kick velocities — the symplectic RATTLE scheme on the constraint
/// manifold.
pub struct VelocityVerlet {
    pub constraints: Option<Constraints>,
}

impl VelocityVerlet {
    pub fn new() -> Self {
        Self { constraints: None }
    }

    pub fn with_constraints(constraints: Constraints) -> Self {
        Self {
            constraints: Some(constraints),
        }
    }
}

impl Default for VelocityVerlet {
    fn default() -> Self {
        Self::new()
    }
}

impl Integrator for VelocityVerlet {
    fn step(
        &mut self,
        system: &mut System,
        model: &dyn ForceModel,
        neighbors: &mut NeighborList,
        dt: f64,
    ) -> Result<Energies> {
        Ok(self
            .advance(
                system,
                model,
                neighbors,
                dt,
                None,
                ConstraintSolverOptions::DEFAULT,
            )?
            .energies)
    }

    fn step_with_stress(
        &mut self,
        system: &mut System,
        model: &dyn ForceModel,
        neighbors: &mut NeighborList,
        dt: f64,
        initial: &ModelEvaluation,
    ) -> Result<StepEvaluation> {
        self.step_with_stress_options(
            system,
            model,
            neighbors,
            dt,
            initial,
            ConstraintSolverOptions::DEFAULT,
        )
    }

    fn step_with_stress_options(
        &mut self,
        system: &mut System,
        model: &dyn ForceModel,
        neighbors: &mut NeighborList,
        dt: f64,
        initial: &ModelEvaluation,
        solver: ConstraintSolverOptions,
    ) -> Result<StepEvaluation> {
        self.advance(system, model, neighbors, dt, Some(initial), solver)
    }
}

impl VelocityVerlet {
    fn advance(
        &self,
        system: &mut System,
        model: &dyn ForceModel,
        neighbors: &mut NeighborList,
        dt: f64,
        initial: Option<&ModelEvaluation>,
        solver: ConstraintSolverOptions,
    ) -> Result<StepEvaluation> {
        validate_step_inputs(system, dt, initial, solver)?;
        let constraints = self.constraints.as_ref();
        prepare_constraints(system, constraints)?;
        let n = system.n();
        let half = 0.5 * dt;
        kick(system, n, half);
        let position = position_drift(
            system,
            constraints,
            dt,
            solver,
            initial.map(|_| ProjectionStage::VvPosition),
        )?;
        system.wrap_positions();
        // This is the straight drift velocity, including SHAKE's displacement / dt.
        let measured = initial.map(|m| (m, kinetic_tensor(system)));
        system.forces.fill(0.);
        let mut virial = Virial::default();
        let (pe, measured) = if let Some((initial, kinetic_average)) = measured {
            let evaluation = model.evaluate(system, neighbors)?;
            virial.w = evaluation.scalar_virial;
            (
                evaluation.energy,
                Some((initial, kinetic_average, evaluation)),
            )
        } else {
            (model.compute(system, &mut virial, neighbors)?, None)
        };
        kick(system, n, half);
        let final_projection = project_tangent(
            system,
            constraints,
            solver,
            initial.map(|_| ProjectionStage::VvFinal),
        )?;
        let stress = if let Some((initial, kinetic_average, final_model)) = measured {
            let stages: Vec<_> = position.into_iter().chain(final_projection).collect();
            let impulse = stages
                .iter()
                .fold(SymTensor6::ZERO, |sum, r| sum.plus(r.tensor_impulse));
            let constraint_average = SymTensor6(impulse.0.map(|x| x / dt));
            let model_average = average_model_tensors(&initial.tensor, &final_model.tensor);
            let volume = system.box_size.volume();
            let (total, scalar) =
                pressure_total(kinetic_average, &model_average, constraint_average, volume)?;
            let interval = IntervalPressure {
                kinetic_average,
                model_average,
                constraint_average,
                total,
                scalar,
                volume,
                dt,
                validation: ValidationStatus::Diagnostic,
            };
            let endpoint = endpoint_pressure(system, &final_model, constraints)?;
            Some(StepStress {
                solver,
                final_model,
                interval,
                endpoint,
                stages,
                bath: None,
            })
        } else {
            None
        };
        Ok(StepEvaluation {
            energies: Energies::from(system, pe, &virial),
            stress,
        })
    }
}

fn validate_step_inputs(
    system: &System,
    dt: f64,
    initial: Option<&ModelEvaluation>,
    solver: ConstraintSolverOptions,
) -> Result<()> {
    solver.validate()?;
    if !dt.is_finite() || dt <= 0. {
        return Err(crate::ModenaError::Config(
            "timestep must be finite and positive".into(),
        ));
    }
    let volume = system.box_size.volume();
    if !volume.is_finite() || volume <= 0. {
        return Err(crate::ModenaError::Physics(
            "volume must be finite and positive".into(),
        ));
    }
    if let Some(m) = initial {
        if !m.energy.is_finite()
            || !m.scalar_virial.is_finite()
            || matches!(m.tensor,TensorAvailability::Complete(t) if !t.finite())
        {
            return Err(crate::ModenaError::Physics(
                "nonfinite initial model evaluation".into(),
            ));
        }
    }
    Ok(())
}

fn average_model_tensors(
    initial: &TensorAvailability,
    final_value: &TensorAvailability,
) -> TensorAvailability {
    use TensorAvailability::{Complete, Missing};
    match (initial, final_value) {
        (Complete(a), Complete(b)) => Complete(a.scaled(0.5).plus(b.scaled(0.5))),
        (Missing(a), Missing(b)) => Missing(a.iter().chain(b).cloned().collect()),
        (Missing(a), _) | (_, Missing(a)) => Missing(a.clone()),
    }
}

/// BAOAB Langevin velocity-Verlet (Leimkuhler & Matthews, 2013). NVT.
///
/// Ordering per step (with F(t) cached): B½ A½ O A½ B½(recompute F).
/// The O step applies exact Ornstein–Uhlenbeck relaxation to velocities:
///   v ← c1·v + sqrt((1−c1²)·kT/m)·N(0,1),  c1 = exp(−γ·dt).
/// Seeded for reproducibility — important for research-grade runs.
pub struct LangevinBAOAB {
    pub gamma: f64,
    pub target_temperature: f64,
    rng: ChaCha12Rng,
    seed: u64,
    pub constraints: Option<Constraints>,
}

impl LangevinBAOAB {
    pub fn new(gamma: f64, target_temperature: f64, seed: u64) -> Self {
        Self {
            gamma,
            target_temperature,
            rng: ChaCha12Rng::seed_from_u64(seed),
            seed,
            constraints: None,
        }
    }

    pub fn with_constraints(
        gamma: f64,
        target_temperature: f64,
        seed: u64,
        constraints: Constraints,
    ) -> Self {
        Self {
            gamma,
            target_temperature,
            rng: ChaCha12Rng::seed_from_u64(seed),
            seed,
            constraints: Some(constraints),
        }
    }
}

impl Integrator for LangevinBAOAB {
    fn checkpoint_state(&self) -> Result<serde_json::Value> {
        if let Some(cons) = &self.constraints {
            cons.validate_parameters()?;
        }
        if !self.gamma.is_finite() || !self.target_temperature.is_finite() {
            return Err(crate::thermostat::state_error(
                "nonfinite Langevin parameters",
            ));
        }
        Ok(
            serde_json::json!({"kind":"langevin_baoab","version":1,"gamma":self.gamma,
            "target_temperature":self.target_temperature,"seed":self.seed,"constraints":self.constraints,
            "rng":crate::thermostat::checkpoint_rng(&self.rng)?}),
        )
    }
    fn restore_state(&mut self, state: &serde_json::Value) -> Result<()> {
        crate::thermostat::validate_state_parameters(
            state,
            &self.checkpoint_state()?,
            &[
                "kind",
                "version",
                "gamma",
                "target_temperature",
                "seed",
                "constraints",
            ],
        )?;
        self.rng = crate::thermostat::restore_rng(&state["rng"], self.rng.get_seed())?;
        Ok(())
    }
    fn step(
        &mut self,
        system: &mut System,
        model: &dyn ForceModel,
        neighbors: &mut NeighborList,
        dt: f64,
    ) -> Result<Energies> {
        system.set_conserved_com_dof(0)?;
        prepare_constraints(system, self.constraints.as_ref())?;
        let n = system.n();
        let half = 0.5 * dt;
        let c1 = (-self.gamma * dt).exp();
        let kt = self.target_temperature;
        let c1sq = c1 * c1;

        // B½: v += (dt/2) F/m   (using cached F(t))
        kick(system, n, half);

        if let Some(cons) = &self.constraints {
            cons.rattle(system, SHAKE_TOL, SHAKE_MAX_ITER)?;
        }
        // A half-step on the constraint manifold, with only its projection
        // displacement applied as a velocity correction.
        project_drift(system, self.constraints.as_ref(), half, true)?;
        // O: exact OU velocity update.
        {
            let v = &mut system.velocities;
            let m = &system.masses;
            for (i, &mi) in m.iter().enumerate() {
                let s = ((1.0 - c1sq) * kt / mi).sqrt();
                for k in 0..3 {
                    let gauss: f64 = StandardNormal.sample(&mut self.rng);
                    let idx = 3 * i + k;
                    v[idx] = c1 * v[idx] + s * gauss;
                }
            }
        }
        // OU noise must be tangent, preserving the mass-weighted Gaussian
        // law of every unconstrained mode (including molecular COM).
        if let Some(cons) = &self.constraints {
            cons.rattle(system, SHAKE_TOL, SHAKE_MAX_ITER)?;
        }
        project_drift(system, self.constraints.as_ref(), half, true)?;

        // Recompute F(t+dt). Zero forces first (M3 accumulate contract).
        system.forces.fill(0.0);
        let mut virial = Virial::default();
        let pe = model.compute(system, &mut virial, neighbors)?;

        // B½: v += (dt/2) F/m
        kick(system, n, half);

        // RATTLE: project velocities onto the constraint tangent.
        if let Some(cons) = &self.constraints {
            cons.rattle(system, SHAKE_TOL, SHAKE_MAX_ITER)?;
        }

        Ok(Energies::from(system, pe, &virial))
    }
}

/// Record the constraint manifold before kinetic observables are used.
pub(crate) fn prepare_constraints(system: &mut System, cons: Option<&Constraints>) -> Result<()> {
    if let Some(cons) = cons {
        cons.validate(system)?;
    }
    system.set_constraint_count(cons.map_or(0, Constraints::count))
}

/// Drift and add only the position projection's velocity impulse. Measuring
/// the unwrapped correction avoids replacing the preceding OU update, and
/// leaves atoms absent from the constraint set exactly unchanged.
pub(crate) fn project_drift(
    system: &mut System,
    cons: Option<&Constraints>,
    interval: f64,
    tangent: bool,
) -> Result<()> {
    position_drift(
        system,
        cons,
        interval,
        ConstraintSolverOptions::DEFAULT,
        None,
    )?;
    if tangent {
        project_tangent(system, cons, ConstraintSolverOptions::DEFAULT, None)?;
    }
    system.wrap_positions();
    Ok(())
}

/// Freeze the image before prediction and leave the corrected drift unwrapped.
fn position_drift(
    system: &mut System,
    cons: Option<&Constraints>,
    interval: f64,
    solver: ConstraintSolverOptions,
    stage: Option<ProjectionStage>,
) -> Result<Option<ProjectionReport>> {
    let Some(cons) = cons else {
        drift(system, system.n(), interval);
        return Ok(None);
    };
    let geometry = constraint_geometry(system, cons)?;
    drift(system, system.n(), interval);
    let predicted = system.positions.clone();
    let report = if let Some(stage) = stage {
        Some(cons.shake_report(
            system,
            &geometry,
            interval,
            stage,
            solver.tolerance,
            solver.max_iterations,
        )?)
    } else {
        cons.shake_at_geometry(system, &geometry, solver.tolerance, solver.max_iterations)?;
        None
    };
    for ((v, x), predicted) in system
        .velocities
        .iter_mut()
        .zip(&system.positions)
        .zip(predicted)
    {
        *v += (*x - predicted) / interval;
    }
    Ok(report)
}

fn project_tangent(
    system: &mut System,
    cons: Option<&Constraints>,
    solver: ConstraintSolverOptions,
    stage: Option<ProjectionStage>,
) -> Result<Option<ProjectionReport>> {
    let Some(cons) = cons else {
        return Ok(None);
    };
    if let Some(stage) = stage {
        let geometry = constraint_geometry(system, cons)?;
        Ok(Some(cons.rattle_report(
            system,
            &geometry,
            stage,
            solver.tolerance,
            solver.max_iterations,
        )?))
    } else {
        cons.rattle(system, solver.tolerance, solver.max_iterations)?;
        Ok(None)
    }
}

/// Half-kick: `v += factor * F / m` over all atoms (flat SoA).
#[inline]
fn kick(system: &mut System, n: usize, factor: f64) {
    let v = &mut system.velocities;
    let f = &system.forces;
    let m = &system.masses;
    for i in 0..n {
        let inv_m = factor / m[i];
        for k in 0..3 {
            v[3 * i + k] += f[3 * i + k] * inv_m;
        }
    }
}

/// Drift: `x += factor * v` over all atoms (flat SoA).
#[inline]
fn drift(system: &mut System, n: usize, factor: f64) {
    let x = &mut system.positions;
    let v = &system.velocities;
    for k in 0..3 * n {
        x[k] += v[k] * factor;
    }
}

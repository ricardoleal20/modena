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

use crate::constraint_stress::{
    constraint_geometry, BathImpulse, ProjectionReport, ProjectionStage,
};
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
    /// Once built-in VV or BAOAB starts advancing, an error is terminal for that step.
    /// The system, force cache or integrator may have advanced partway. Restore a known
    /// state before retrying. Only a successful step returns a completed record.
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
    /// The matching input cache and partial-state error contract of
    /// [`Self::step_with_stress`] also apply when selecting solver options.
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
            let interval = interval_pressure(
                system,
                dt,
                kinetic_average,
                initial,
                &final_model,
                stages.iter(),
            )?;
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

fn interval_pressure<'a>(
    system: &System,
    dt: f64,
    kinetic_average: SymTensor6,
    initial: &ModelEvaluation,
    final_model: &ModelEvaluation,
    mechanical: impl Iterator<Item = &'a ProjectionReport>,
) -> Result<IntervalPressure> {
    let impulse = mechanical.fold(SymTensor6::ZERO, |sum, r| sum.plus(r.tensor_impulse));
    let constraint_average = SymTensor6(impulse.0.map(|x| x / dt));
    let model_average = average_model_tensors(&initial.tensor, &final_model.tensor);
    let volume = system.box_size.volume();
    let (total, scalar) =
        pressure_total(kinetic_average, &model_average, constraint_average, volume)?;
    Ok(IntervalPressure {
        kinetic_average,
        model_average,
        constraint_average,
        total,
        scalar,
        volume,
        dt,
        validation: ValidationStatus::Diagnostic,
    })
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

impl LangevinBAOAB {
    fn advance(
        &mut self,
        system: &mut System,
        model: &dyn ForceModel,
        neighbors: &mut NeighborList,
        dt: f64,
        initial: Option<&ModelEvaluation>,
        solver: ConstraintSolverOptions,
    ) -> Result<StepEvaluation> {
        use ProjectionStage::*;
        validate_step_inputs(system, dt, initial, solver)?;
        system.set_conserved_com_dof(0)?;
        let cons = self.constraints.as_ref();
        prepare_constraints(system, cons)?;
        let half = 0.5 * dt;
        let n = system.n();
        let mut stages = Vec::new();
        let measured = initial.is_some();
        kick(system, n, half);
        let b0 = project_tangent(system, cons, solver, initial.map(|_| BaoabB0))?;
        if measured {
            stages.push(b0.unwrap_or_else(|| ProjectionReport::new(BaoabB0)));
        }
        let k1 = baoab_drift(
            system,
            cons,
            half,
            solver,
            initial.map(|_| (BaoabA1Pre, BaoabA1Post)),
            &mut stages,
        )?;
        let bath = projected_ou(
            system,
            cons,
            &mut self.rng,
            (-self.gamma * dt).exp(),
            self.target_temperature,
            solver,
            measured,
        )?;
        let bath = bath.map(|(bath, projection)| {
            stages.push(projection);
            bath
        });
        let k2 = baoab_drift(
            system,
            cons,
            half,
            solver,
            initial.map(|_| (BaoabA2Pre, BaoabA2Post)),
            &mut stages,
        )?;

        system.forces.fill(0.);
        let mut virial = Virial::default();
        let (pe, models) = if let Some(initial) = initial {
            let final_model = model.evaluate(system, neighbors)?;
            virial.w = final_model.scalar_virial;
            (final_model.energy, Some((initial, final_model)))
        } else {
            (model.compute(system, &mut virial, neighbors)?, None)
        };
        kick(system, n, half);
        let b1 = project_tangent(system, cons, solver, initial.map(|_| BaoabB1))?;
        let stress = if let Some((initial, final_model)) = models {
            stages.push(b1.unwrap_or_else(|| ProjectionReport::new(BaoabB1)));
            let kinetic = k1.scaled(0.5).plus(k2.scaled(0.5));
            let mechanical = stages.iter().filter(|r| r.stage != BaoabO);
            let interval =
                interval_pressure(system, dt, kinetic, initial, &final_model, mechanical)?;
            let endpoint = endpoint_pressure(system, &final_model, cons)?;
            Some(StepStress {
                solver,
                final_model,
                interval,
                endpoint,
                stages,
                bath,
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

/// Sample straight drift momentum after SHAKE but before the tangent impulse.
fn baoab_drift(
    system: &mut System,
    cons: Option<&Constraints>,
    half: f64,
    solver: ConstraintSolverOptions,
    tags: Option<(ProjectionStage, ProjectionStage)>,
    stages: &mut Vec<ProjectionReport>,
) -> Result<SymTensor6> {
    let pre = position_drift(system, cons, half, solver, tags.map(|(pre, _)| pre))?;
    let kinetic = tags.map_or(SymTensor6::ZERO, |_| kinetic_tensor(system));
    let post = project_tangent(system, cons, solver, tags.map(|(_, post)| post))?;
    if let Some((a, b)) = tags {
        stages.push(pre.unwrap_or_else(|| ProjectionReport::new(a)));
        stages.push(post.unwrap_or_else(|| ProjectionReport::new(b)));
    }
    system.wrap_positions();
    Ok(kinetic)
}

/// Exact scalar OU transition followed by the mass-weighted tangent projection at fixed q.
fn projected_ou(
    system: &mut System,
    cons: Option<&Constraints>,
    rng: &mut ChaCha12Rng,
    c: f64,
    kt: f64,
    solver: ConstraintSolverOptions,
    measured: bool,
) -> Result<Option<(BathImpulse, ProjectionReport)>> {
    use crate::Vec3;
    let csq = c * c;
    let mut unprojected = Vec::new();
    let mut unprojected_moment = SymTensor6::ZERO;
    for i in 0..system.n() {
        let mass = system.masses[i];
        let sigma = ((1.0 - csq) * kt / mass).sqrt();
        let v = &mut system.velocities[3 * i..3 * i + 3];
        let mut impulse = Vec3::zero();
        for (k, value) in v.iter_mut().enumerate() {
            // Consume every normal even at zero amplitude for RNG continuation.
            let normal: f64 = StandardNormal.sample(rng);
            let before = *value;
            *value = c * before + sigma * normal;
            if measured {
                impulse.0[k] = mass * (*value - before);
            }
        }
        if measured {
            unprojected.push(impulse);
            unprojected_moment = unprojected_moment.plus(SymTensor6::dyad(system.pos(i), impulse));
        }
    }
    let before_projection = measured.then(|| system.velocities.clone());
    let stage = ProjectionStage::BaoabO;
    let projection = project_tangent(system, cons, solver, measured.then_some(stage))?;
    let Some(before) = before_projection else {
        return Ok(None);
    };
    let projection = projection.unwrap_or_else(|| ProjectionReport::new(stage));
    let reaction: Vec<_> = before
        .chunks_exact(3)
        .enumerate()
        .map(|(i, v)| (system.vel(i) - Vec3::new(v[0], v[1], v[2])) * system.masses[i])
        .collect();
    if !unprojected_moment.finite()
        || reaction
            .iter()
            .chain(&unprojected)
            .any(|j| j.0.iter().any(|v| !v.is_finite()))
    {
        return Err(crate::ModenaError::Physics(
            "nonfinite OU bath impulse".into(),
        ));
    }
    let bath = BathImpulse {
        unprojected,
        reaction,
        unprojected_moment,
        reaction_moment: projection.tensor_impulse,
    };
    Ok(Some((bath, projection)))
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

#[cfg(test)]
#[path = "../tests/support/conditional_ou.rs"]
mod conditional_ou;

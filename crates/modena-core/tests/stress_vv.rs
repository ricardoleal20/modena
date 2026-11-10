//! Finite-step VV pressure against independently evaluated rotor and moment oracles.
#[allow(dead_code)]
#[path = "support/stress.rs"]
mod support;
use modena_core::stress::{ConstraintSolverOptions, ModelEvaluation};
use modena_core::{ForceModel, Integrator, NeighborList, Result, System, VelocityVerlet, Virial};
use support::{near, rotor};

struct ZeroModel;
impl ForceModel for ZeroModel {
    fn compute(&self, _: &mut System, _: &mut Virial, _: &mut NeighborList) -> Result<f64> {
        Ok(0.)
    }
    fn evaluate(&self, _: &mut System, _: &mut NeighborList) -> Result<ModelEvaluation> {
        Ok(ModelEvaluation::zero())
    }
}

#[test]
fn vv_interval_uses_drift_kinetic_not_endpoint_kinetic() {
    let (mut s, cons) = rotor(1., 3.);
    let h = 0.08 / 0.7;
    let mut vv = VelocityVerlet::with_constraints(cons);
    let out = vv
        .step_with_stress_options(
            &mut s,
            &ZeroModel,
            &mut NeighborList::new(1.),
            h,
            &ModelEvaluation::zero(),
            ConstraintSolverOptions::TIGHT,
        )
        .unwrap();
    let stress = out.stress.unwrap();
    near(stress.interval.scalar.unwrap(), 0.56 / 24000., 5e-10);
    let wrong = stress.endpoint.kinetic.trace() + stress.interval.constraint_average.trace();
    let c = (1. - (0.7_f64 * h).powi(2)).sqrt();
    near(wrong - 0.56, -0.621075 * (1. - c) / (1. + c), 5e-10);
    assert!((wrong - 0.56).abs() > 1e-4);
}

use modena_core::constraint_stress::ProjectionStage;
use modena_core::stress::{
    IntervalPressure, MissingTensor, MissingTensorReason, StepStress, SymTensor6,
    TensorAvailability, ValidationStatus,
};
use modena_core::{Constraints, Energies, ModenaError, SimBox, Vec3};
use std::sync::atomic::{AtomicUsize, Ordering};
use support::complete;

const CENTER: Vec3 = Vec3::new(5., 6., 7.);
const COM_V: Vec3 = Vec3::new(0.2, -0.1, 0.3);
const H: f64 = 0.08 / 0.7;
const IMPULSE_TOL: f64 = 5e-10;
const ENDPOINT_TOL: f64 = 2e-9;

fn dyad(a: Vec3, b: Vec3) -> [f64; 6] {
    [(0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2)]
        .map(|(i, j)| (a.0[i] * b.0[j] + a.0[j] * b.0[i]) * 0.5)
}
fn check_tensor(actual: SymTensor6, expected: [f64; 6], tol: f64) {
    for (got, want) in actual.0.into_iter().zip(expected) {
        near(got, want, tol);
    }
}
fn check_vector(actual: Vec3, expected: Vec3, tol: f64) {
    for (got, want) in actual.0.into_iter().zip(expected.0) {
        near(got, want, tol);
    }
}
fn moment(s: &System) -> [f64; 6] {
    let mut g = [0.; 6];
    for i in 0..s.n() {
        let term = dyad(s.pos(i), s.vel(i) * s.masses[i]);
        for k in 0..6 {
            g[k] += term[k];
        }
    }
    g
}
fn check_delta_g(before: &System, after: &System, interval: &IntervalPressure) {
    let g0 = moment(before);
    let g1 = moment(after);
    let total = complete(&interval.total);
    let kinetic = interval.kinetic_average.0;
    let model = complete(&interval.model_average).0;
    let constraint = interval.constraint_average.0;
    for k in 0..6 {
        let ledger = interval.dt * interval.volume * total.0[k];
        let scale = 1. + interval.dt * (kinetic[k].abs() + model[k].abs() + constraint[k].abs());
        assert!(
            (g1[k] - g0[k] - ledger).abs() <= IMPULSE_TOL * scale,
            "component {k}: delta G = {}, ledger = {ledger}",
            g1[k] - g0[k]
        );
    }
}

#[test]
fn finite_rotor_map_stages_and_six_component_identity_across_images() {
    let e = Vec3::new(1., 2., 2.) * (1. / 3.);
    let t = Vec3::new(2., -1., 0.) * (1. / 5_f64.sqrt());
    let a = 0.7 * H;
    let c = (1. - a * a).sqrt();
    let e1 = e * c + t * a;
    let r1 = e1 * 1.3;
    let ua = (e * (c - 1.) + t * a) * (1.3 / H);
    let u1 = (-e * a + t * c) * (1.3 * 0.7);
    for (m1, m2) in [(1., 3.), (2., 5.)] {
        let (base, cons) = rotor(m1, m2);
        let mass = m1 + m2;
        let mu = m1 * m2 / mass;
        let ipos = dyad(e, e).map(|x| x * mu * 1.3_f64.powi(2) * (c - 1.) / H);
        let ifinal = dyad(e1, e1).map(|x| -x * mu * 1.3_f64.powi(2) * (1. - c) / H);
        let com = dyad(COM_V, COM_V).map(|x| x * mass);
        let internal = dyad(ua, ua).map(|x| x * mu);
        let kin = std::array::from_fn(|k| com[k] + internal[k]);
        let kend = dyad(u1, u1).map(|x| x * mu);
        let reaction = dyad(r1, r1).map(|x| -mu * 0.7_f64.powi(2) * x);
        for x in [-1, 0, 1] {
            for y in [-1, 0, 1] {
                for z in [-1, 0, 1] {
                    let axes = [x, y, z];
                    let center = Vec3(std::array::from_fn(|i| match axes[i] {
                        -1 => 0.,
                        1 => 20.,
                        _ => CENTER.0[i],
                    }));
                    for shift in [Vec3::zero(), Vec3::new(0.31, -0.27, 0.19)] {
                        let mut s = base.clone();
                        for i in 0..2 {
                            s.set_pos(i, s.pos(i) - CENTER + center + shift);
                        }
                        let unwrapped = s.clone();
                        s.wrap_positions();
                        let wrapped = s.clone();
                        let start = std::time::Instant::now();
                        let out = VelocityVerlet::with_constraints(cons.clone())
                            .step_with_stress_options(
                                &mut s,
                                &ZeroModel,
                                &mut NeighborList::new(1.),
                                H,
                                &ModelEvaluation::zero(),
                                ConstraintSolverOptions::TIGHT,
                            )
                            .unwrap();
                        let ns = start.elapsed().as_nanos();
                        let report = out.stress.unwrap();
                        assert_eq!(report.solver, ConstraintSolverOptions::TIGHT);
                        assert_eq!(report.endpoint.validation, ValidationStatus::Diagnostic);
                        assert_eq!(report.interval.validation, ValidationStatus::Diagnostic);
                        assert!(report.bath.is_none());
                        let [position, final_proj] = report.stages.as_slice() else {
                            panic!("expected exactly two VV projections")
                        };
                        assert_eq!(position.stage, ProjectionStage::VvPosition);
                        assert_eq!(final_proj.stage, ProjectionStage::VvFinal);
                        assert!(
                            position.max_position_residual <= 1e-13
                                && final_proj.max_velocity_residual <= 1e-13
                        );
                        for (stage, lever, j) in [
                            (position, e * 1.3, e * (mu * 1.3 * (c - 1.) / H)),
                            (final_proj, r1, e1 * (-mu * 1.3 * (1. - c) / H)),
                        ] {
                            let mut sum = Vec3::zero();
                            let mut torque = Vec3::zero();
                            for update in &stage.impulses {
                                assert_eq!(update.edge, 0);
                                check_vector(update.lever, lever, IMPULSE_TOL);
                                sum += update.impulse_on_j;
                                torque += update.lever.cross(update.impulse_on_j);
                            }
                            check_vector(sum, j, IMPULSE_TOL);
                            check_vector(torque, Vec3::zero(), IMPULSE_TOL);
                        }
                        check_tensor(position.tensor_impulse, ipos, IMPULSE_TOL);
                        check_tensor(final_proj.tensor_impulse, ifinal, IMPULSE_TOL);
                        check_tensor(
                            report.interval.constraint_average,
                            std::array::from_fn(|k| (ipos[k] + ifinal[k]) / H),
                            IMPULSE_TOL,
                        );
                        check_tensor(report.interval.kinetic_average, kin, IMPULSE_TOL);
                        check_tensor(
                            report.endpoint.kinetic,
                            std::array::from_fn(|k| com[k] + kend[k]),
                            ENDPOINT_TOL,
                        );
                        check_tensor(report.endpoint.constraint, reaction, ENDPOINT_TOL);
                        check_vector(s.total_momentum(), COM_V * mass, IMPULSE_TOL);
                        // Lift the tiny step into the continuous initial molecular image for G.
                        for i in 0..2 {
                            let delta = s.box_size.min_image(s.pos(i) - wrapped.pos(i));
                            s.set_pos(i, unwrapped.pos(i) + delta);
                            let fraction = if i == 0 { -m2 / mass } else { m1 / mass };
                            check_vector(delta, (COM_V + ua * fraction) * H, IMPULSE_TOL);
                            check_vector(s.vel(i), COM_V + u1 * fraction, IMPULSE_TOL);
                        }
                        check_vector(s.pos(1) - s.pos(0), r1, IMPULSE_TOL);
                        check_delta_g(&unwrapped, &s, &report.interval);
                        if std::env::var_os("B07_VV_EVIDENCE").is_some()
                            && axes == [0, 0, 0]
                            && shift == Vec3::zero()
                        {
                            println!(
                                "{}",
                                serde_json::json!({"masses":[m1,m2],"step_wall_ns":ns,"report":report,"before_G":moment(&unwrapped),"after_G":moment(&s),"expected_position_impulse":ipos,"expected_final_impulse":ifinal,"expected_drift_kinetic":kin,"serialized_bytes":serde_json::to_vec(&report).unwrap().len(),"impulse_item_bytes":std::mem::size_of::<modena_core::constraint_stress::ConstraintImpulse>()})
                            );
                        }
                    }
                }
            }
        }
    }
}

// Frozen external harmonic trap. Its configurational moment uses absolute q,
// the same origin as G. The independent force/tension oracle below uses F=-0.1(q-center).
#[derive(Default)]
struct Trap {
    calls: AtomicUsize,
}
impl Trap {
    fn apply(&self, s: &mut System) -> ModelEvaluation {
        self.calls.fetch_add(1, Ordering::SeqCst);
        let mut model = ModelEvaluation::zero();
        let mut tensor = [0.; 6];
        for i in 0..s.n() {
            let q = s.pos(i);
            let d = q - CENTER;
            let f = d * (-0.1);
            for k in 0..3 {
                s.forces[3 * i + k] += f.0[k];
            }
            model.energy += 0.05 * d.norm2();
            model.scalar_virial += q.dot(f);
            let term = dyad(q, f);
            for k in 0..6 {
                tensor[k] += term[k];
            }
        }
        model.tensor = TensorAvailability::Complete(SymTensor6(tensor));
        model
    }
}
impl ForceModel for Trap {
    fn compute(&self, s: &mut System, w: &mut Virial, _: &mut NeighborList) -> Result<f64> {
        let m = self.apply(s);
        w.w += m.scalar_virial;
        Ok(m.energy)
    }
    fn evaluate(&self, s: &mut System, _: &mut NeighborList) -> Result<ModelEvaluation> {
        Ok(self.apply(s))
    }
}

#[test]
fn forced_rotor_and_unconstrained_moments_use_both_conservative_evaluations() {
    for constrained in [true, false] {
        let (mut s, cons) = rotor(2., 5.);
        let mut vv = if constrained {
            VelocityVerlet::with_constraints(cons)
        } else {
            VelocityVerlet::new()
        };
        let trap = Trap::default();
        let mut nl = NeighborList::new(1.);
        let mut initial = trap.evaluate(&mut s, &mut nl).unwrap();
        for step in 1..=8 {
            let before = s.clone();
            let out = vv
                .step_with_stress_options(
                    &mut s,
                    &trap,
                    &mut nl,
                    H,
                    &initial,
                    ConstraintSolverOptions::TIGHT,
                )
                .unwrap();
            let stress = out.stress.unwrap();
            assert_eq!(trap.calls.load(Ordering::SeqCst), step + 1);
            for i in 0..2 {
                check_vector(s.force(i), (s.pos(i) - CENTER) * (-0.1), 1e-14);
            }
            let r = s.pos(1) - s.pos(0);
            let u = s.vel(1) - s.vel(0);
            let b_force = r.dot(s.force(1) * (1. / 5.) - s.force(0) * (1. / 2.));
            assert!(b_force.abs() > 1e-3);
            let lambda = if constrained {
                (u.norm2() + b_force) / (r.norm2() * (1. / 2. + 1. / 5.))
            } else {
                0.
            };
            check_tensor(
                stress.endpoint.constraint,
                dyad(r, r).map(|x| -lambda * x),
                ENDPOINT_TOL,
            );
            if !constrained {
                assert!(stress.stages.is_empty());
                assert_eq!(stress.interval.constraint_average, SymTensor6::ZERO);
            }
            let s0 = complete(&initial.tensor).0;
            let s1 = complete(&stress.final_model.tensor).0;
            assert!((s0[0] - s1[0]).abs() > 1e-6);
            check_tensor(
                complete(&stress.interval.model_average),
                std::array::from_fn(|k| (s0[k] + s1[k]) / 2.),
                IMPULSE_TOL,
            );
            check_delta_g(&before, &s, &stress.interval);
            initial = stress.final_model;
        }
    }
}

// Legacy contract at the pre-extension source epoch, using the public projection wrappers.
fn reference_vv(
    s: &mut System,
    cons: &Constraints,
    model: &dyn ForceModel,
    nl: &mut NeighborList,
) -> Energies {
    let kick = |s: &mut System| {
        for i in 0..s.n() {
            let inv = 0.5 * H / s.masses[i];
            for k in 0..3 {
                s.velocities[3 * i + k] += s.forces[3 * i + k] * inv;
            }
        }
    };
    s.set_constraint_count(cons.count()).unwrap();
    kick(s);
    let reference = s.positions.clone();
    for (q, v) in s.positions.iter_mut().zip(&s.velocities) {
        *q += v * H;
    }
    let predicted = s.positions.clone();
    cons.shake(s, &reference, 1e-10, 100).unwrap();
    for ((v, q), p) in s.velocities.iter_mut().zip(&s.positions).zip(predicted) {
        *v += (q - p) / H;
    }
    s.wrap_positions();
    s.forces.fill(0.);
    let mut virial = Virial::default();
    let pe = model.compute(s, &mut virial, nl).unwrap();
    kick(s);
    cons.rattle(s, 1e-10, 100).unwrap();
    Energies::from(s, pe, &virial)
}

#[test]
fn default_steps_match_the_legacy_trajectory_energies_and_force_call_count() {
    let (mut s, cons) = rotor(1., 3.);
    let trap = Trap::default();
    let mut nl = NeighborList::new(1.);
    let mut initial = trap.evaluate(&mut s, &mut nl).unwrap();
    let mut legacy = s.clone();
    let mut reference = s.clone();
    let mut vv = VelocityVerlet::with_constraints(cons.clone());
    for n in 1..=8 {
        let expected = reference_vv(&mut reference, &cons, &trap, &mut nl);
        let e = vv.step(&mut legacy, &trap, &mut nl, H).unwrap();
        let out = vv
            .step_with_stress(&mut s, &trap, &mut nl, H, &initial)
            .unwrap();
        let stress = out.stress.unwrap();
        assert_eq!(stress.solver, ConstraintSolverOptions::DEFAULT);
        assert_eq!(trap.calls.load(Ordering::SeqCst), 1 + 3 * n);
        for state in [&legacy, &s] {
            assert_eq!(state.positions, reference.positions);
            assert_eq!(state.velocities, reference.velocities);
            assert_eq!(state.forces, reference.forces);
        }
        for energies in [e, out.energies] {
            assert_eq!(
                serde_json::to_value(energies).unwrap(),
                serde_json::to_value(expected).unwrap()
            );
        }
        initial = stress.final_model;
    }
}

#[test]
fn rigid_water_interval_matches_unwrapped_delta_g_for_both_mass_sets() {
    for masses in [[15.9994, 1.008, 1.008], [16., 1., 2.]] {
        let (mut s, cons) = support::water(masses);
        let before = s.clone();
        let out = VelocityVerlet::with_constraints(cons)
            .step_with_stress_options(
                &mut s,
                &ZeroModel,
                &mut NeighborList::new(1.),
                H,
                &ModelEvaluation::zero(),
                ConstraintSolverOptions::TIGHT,
            )
            .unwrap();
        let stress = out.stress.unwrap();
        check_delta_g(&before, &s, &stress.interval);
        near(
            stress.interval.constraint_average.trace() + stress.interval.kinetic_average.trace(),
            s.masses.iter().sum::<f64>() * COM_V.norm2(),
            IMPULSE_TOL,
        );
    }
}

struct ScalarOnly;
impl ForceModel for ScalarOnly {
    fn compute(&self, _: &mut System, _: &mut Virial, _: &mut NeighborList) -> Result<f64> {
        Ok(0.)
    }
}
#[test]
fn unavailable_initial_or_final_model_retains_reasons_and_endpoint_independence() {
    let mut vv = VelocityVerlet::new();
    for mi in [true, false] {
        for mf in [true, false] {
            let (mut s, _) = rotor(1., 3.);
            let mut nl = NeighborList::new(1.);
            let mut initial = ModelEvaluation::zero();
            let first = MissingTensor {
                component: "initial".into(),
                reason: MissingTensorReason::MissingInitialEvaluation,
            };
            let last = MissingTensor {
                component: std::any::type_name::<ScalarOnly>().into(),
                reason: MissingTensorReason::LegacyForceModel,
            };
            if mi {
                initial.tensor = TensorAvailability::Missing(vec![first.clone()]);
            }
            let model: &dyn ForceModel = if mf { &ScalarOnly } else { &ZeroModel };
            let out = vv
                .step_with_stress(&mut s, model, &mut nl, H, &initial)
                .unwrap()
                .stress
                .unwrap();
            if mi || mf {
                let why = TensorAvailability::Missing(
                    mi.then_some(first)
                        .into_iter()
                        .chain(mf.then_some(last))
                        .collect(),
                );
                assert_eq!(out.interval.model_average, why);
                assert_eq!(out.interval.total, why);
                assert!(out.interval.scalar.is_none());
            } else {
                assert!(out.interval.scalar.is_some());
            }
            assert_eq!(out.endpoint.scalar.is_none(), mf);
        }
    }
}

#[test]
fn options_dt_and_volume_fail_before_mutation_and_stage_failure_returns_error() {
    let trap = Trap::default();
    let mut nl = NeighborList::new(1.);
    let initial = ModelEvaluation::zero();
    let default = ConstraintSolverOptions::DEFAULT;
    for (dt, volume_edge, solver) in [
        (0., 20., default),
        (-1., 20., default),
        (f64::NAN, 20., default),
        (f64::INFINITY, 20., default),
        (H, 0., default),
        (H, f64::INFINITY, default),
        (H, f64::NAN, default),
        (
            H,
            20.,
            ConstraintSolverOptions {
                tolerance: 0.,
                ..default
            },
        ),
        (
            H,
            20.,
            ConstraintSolverOptions {
                tolerance: -1.,
                ..default
            },
        ),
        (
            H,
            20.,
            ConstraintSolverOptions {
                tolerance: f64::NAN,
                ..default
            },
        ),
        (
            H,
            20.,
            ConstraintSolverOptions {
                tolerance: f64::INFINITY,
                ..default
            },
        ),
        (
            H,
            20.,
            ConstraintSolverOptions {
                max_iterations: 0,
                ..default
            },
        ),
    ] {
        let (mut s, cons) = rotor(1., 3.);
        s.box_size = SimBox::cubic(volume_edge);
        let before = s.clone();
        let mut vv = VelocityVerlet::with_constraints(cons);
        assert!(vv
            .step_with_stress_options(&mut s, &trap, &mut nl, dt, &initial, solver)
            .is_err());
        assert_eq!(before.positions, s.positions);
        assert_eq!(before.velocities, s.velocities);
        assert_eq!(before.forces, s.forces);
        assert_eq!(s.constraint_count(), before.constraint_count());
    }
    let (mut s, cons) = rotor(1., 3.);
    let err = VelocityVerlet::with_constraints(cons)
        .step_with_stress_options(
            &mut s,
            &trap,
            &mut nl,
            H,
            &initial,
            ConstraintSolverOptions {
                max_iterations: 1,
                ..ConstraintSolverOptions::TIGHT
            },
        )
        .unwrap_err();
    assert!(matches!(err, ModenaError::Physics(_)));
    assert_eq!(trap.calls.load(Ordering::SeqCst), 0);
}

#[test]
fn legacy_custom_integrator_defaults_to_no_stress_and_rejects_custom_options() {
    struct Legacy(usize);
    impl Integrator for Legacy {
        fn step(
            &mut self,
            _: &mut System,
            _: &dyn ForceModel,
            _: &mut NeighborList,
            _: f64,
        ) -> Result<Energies> {
            self.0 += 1;
            Ok(Energies {
                potential: 17.,
                ..Default::default()
            })
        }
    }
    let mut driver = Legacy(0);
    let (mut s, _) = rotor(1., 3.);
    let mut nl = NeighborList::new(1.);
    let m = ModelEvaluation::zero();
    let out = driver
        .step_with_stress_options(
            &mut s,
            &ZeroModel,
            &mut nl,
            H,
            &m,
            ConstraintSolverOptions::DEFAULT,
        )
        .unwrap();
    assert_eq!(out.energies.potential, 17.);
    assert!(out.stress.is_none());
    assert!(matches!(
        driver.step_with_stress_options(
            &mut s,
            &ZeroModel,
            &mut nl,
            H,
            &m,
            ConstraintSolverOptions::TIGHT
        ),
        Err(ModenaError::Config(_))
    ));
    assert_eq!(driver.0, 1);
}

#[test]
fn step_record_preserves_actual_solver_and_deserializes_only_diagnostic_observations() {
    let (mut s, cons) = rotor(1., 3.);
    let solver = ConstraintSolverOptions {
        tolerance: 1e-12,
        max_iterations: 120,
    };
    let report = VelocityVerlet::with_constraints(cons)
        .step_with_stress_options(
            &mut s,
            &ZeroModel,
            &mut NeighborList::new(1.),
            H,
            &ModelEvaluation::zero(),
            solver,
        )
        .unwrap()
        .stress
        .unwrap();
    assert_eq!(report.solver, solver);
    let mut raw = serde_json::to_value(&report).unwrap();
    for field in ["endpoint", "interval"] {
        raw[field]["validation"] =
            serde_json::json!({"kind":"validated","value":{"pub_record":"unverified"}});
    }
    let restored: StepStress = serde_json::from_value(raw.clone()).unwrap();
    assert_eq!(restored.endpoint.validation, ValidationStatus::Diagnostic);
    assert_eq!(restored.interval.validation, ValidationStatus::Diagnostic);
    assert_eq!(restored.solver, solver);
    for (field, value) in [
        ("dt", serde_json::json!(0)),
        ("volume", serde_json::json!(-1)),
        ("scalar", serde_json::Value::Null),
        ("unknown", serde_json::json!(0)),
    ] {
        let mut bad = raw.clone();
        bad["interval"][field] = value;
        assert!(serde_json::from_value::<StepStress>(bad).is_err());
    }
    for (field, value) in [
        ("tolerance", serde_json::json!(0)),
        ("max_iterations", serde_json::json!(0)),
        ("unknown", serde_json::json!(0)),
    ] {
        let mut bad = raw.clone();
        bad["solver"][field] = value;
        assert!(serde_json::from_value::<StepStress>(bad).is_err());
    }
}

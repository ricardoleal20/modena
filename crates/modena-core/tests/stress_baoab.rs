//! Projected Langevin stage ledger, deterministic rotor and bath exchange.
#[allow(dead_code)]
#[path = "support/stress.rs"]
mod support;
use modena_core::stress::{ConstraintSolverOptions, ModelEvaluation};
use modena_core::{ForceModel, Integrator, LangevinBAOAB, NeighborList, Result, System, Virial};
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

use modena_core::constraint_stress::ProjectionStage::*;
use modena_core::stress::{SymTensor6, TensorAvailability, ValidationStatus};
use modena_core::Vec3;
use std::sync::atomic::{AtomicUsize, Ordering};
use support::complete;

const H: f64 = 0.08 / 0.7;
const TOL: f64 = 5e-10;
const COM_V: Vec3 = Vec3::new(0.2, -0.1, 0.3);

fn outer(a: Vec3, b: Vec3) -> SymTensor6 {
    SymTensor6(
        [(0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2)]
            .map(|(i, j)| (a.0[i] * b.0[j] + a.0[j] * b.0[i]) * 0.5),
    )
}
fn tensor_near(a: SymTensor6, b: SymTensor6) {
    for (x, y) in a.0.into_iter().zip(b.0) {
        near(x, y, TOL);
    }
}
fn vec_near(a: Vec3, b: Vec3) {
    for (x, y) in a.0.into_iter().zip(b.0) {
        near(x, y, TOL);
    }
}
fn moment(s: &System) -> SymTensor6 {
    (0..s.n()).fold(SymTensor6::ZERO, |sum, i| {
        sum.plus(outer(s.pos(i), s.vel(i) * s.masses[i]))
    })
}

#[test]
fn two_half_drifts_cancel_rotational_interval_trace() {
    let tau = H / 2.;
    let a = 0.7 * tau;
    let c = (1. - a * a).sqrt();
    for (m1, m2) in [(1., 3.), (2., 5.)] {
        let (mut s, cons) = rotor(m1, m2);
        let mut e = (s.pos(1) - s.pos(0)) * (1. / 1.3);
        let mut t = (s.vel(1) - s.vel(0)) * (1. / (1.3 * 0.7));
        let m = m1 + m2;
        let mu = m1 * m2 / m;
        let mut kinetic = outer(COM_V, COM_V).scaled(m);
        let out = LangevinBAOAB::with_constraints(0., 1., 1001, cons)
            .step_with_stress_options(
                &mut s,
                &ZeroModel,
                &mut NeighborList::new(1.),
                H,
                &ModelEvaluation::zero(),
                ConstraintSolverOptions::TIGHT,
            )
            .unwrap();
        let st = out.stress.unwrap();
        assert_eq!(
            st.stages.iter().map(|p| p.stage).collect::<Vec<_>>(),
            [
                BaoabB0,
                BaoabA1Pre,
                BaoabA1Post,
                BaoabO,
                BaoabA2Pre,
                BaoabA2Post,
                BaoabB1
            ]
        );
        for (pre, post) in [(1, 2), (4, 5)] {
            let e1 = e * c + t * a;
            let ua = (e * (c - 1.) + t * a) * (1.3 / tau);
            kinetic = kinetic.plus(outer(ua, ua).scaled(0.5 * mu));
            let factor = mu * 1.3_f64.powi(2) * (c - 1.) / tau;
            tensor_near(st.stages[pre].tensor_impulse, outer(e, e).scaled(factor));
            tensor_near(st.stages[post].tensor_impulse, outer(e1, e1).scaled(factor));
            for entry in &st.stages[pre].impulses {
                vec_near(entry.lever, e * 1.3);
            }
            for entry in &st.stages[post].impulses {
                vec_near(entry.lever, e1 * 1.3);
            }
            t = t * c - e * a;
            e = e1;
        }
        vec_near(s.pos(1) - s.pos(0), e * 1.3);
        vec_near(s.vel(1) - s.vel(0), t * (1.3 * 0.7));
        tensor_near(st.interval.kinetic_average, kinetic);
        near(st.interval.scalar.unwrap() * 24000., m * 0.14, TOL);
        assert_eq!(st.solver, ConstraintSolverOptions::TIGHT);
        assert_eq!(st.endpoint.validation, ValidationStatus::Diagnostic);
        assert_eq!(st.interval.validation, ValidationStatus::Diagnostic);
        let bath = st.bath.unwrap();
        assert_eq!(bath.unprojected, vec![Vec3::zero(); 2]);
        assert_eq!(bath.reaction, vec![Vec3::zero(); 2]);
        assert_eq!(bath.unprojected_moment, SymTensor6::ZERO);
        assert_eq!(bath.reaction_moment, SymTensor6::ZERO);
        println!(
            "rotor masses={m1},{m2} numerator={} kinetic={:?} reaction={:?}",
            st.interval.scalar.unwrap() * 24000.,
            kinetic.0,
            st.interval.constraint_average.0
        );
    }
}

#[test]
fn baoab_interval_excludes_ou_bath_exchange() {
    use chacha20::ChaCha12Rng;
    use rand::SeedableRng;
    use rand_distr::{Distribution, StandardNormal};
    let (mut s, cons) = rotor(1., 3.);
    let before = moment(&s);
    let h = 0.02;
    let tau = h / 2.;
    let a = 0.7_f64 * tau;
    let cos = (1. - a * a).sqrt();
    let e = (s.pos(1) - s.pos(0)) * (1. / 1.3);
    let t = (s.vel(1) - s.vel(0)) * (1. / (1.3 * 0.7));
    let r = (e * cos + t * a) * 1.3;
    let u = (t * cos - e * a) * (1.3 * 0.7);
    let center = Vec3::new(5., 6., 7.) + COM_V * tau;
    let q = [center - r * 0.75, center + r * 0.25];
    let v = [COM_V - u * 0.75, COM_V + u * 0.25];
    let mut rng = ChaCha12Rng::seed_from_u64(1001);
    let decay = (-h).exp();
    let raw: [Vec3; 2] = std::array::from_fn(|i| {
        let noise = Vec3(std::array::from_fn(|_| StandardNormal.sample(&mut rng)));
        let m = s.masses[i];
        (v[i] * (decay - 1.) + noise * ((1. - decay * decay) / m).sqrt()) * m
    });
    let rel = raw[1] * (1. / 3.) - raw[0];
    let on_j = r * (-0.75 * rel.dot(r) / r.norm2());
    let st = LangevinBAOAB::with_constraints(1., 1., 1001, cons)
        .step_with_stress_options(
            &mut s,
            &ZeroModel,
            &mut NeighborList::new(1.),
            h,
            &ModelEvaluation::zero(),
            ConstraintSolverOptions::TIGHT,
        )
        .unwrap()
        .stress
        .unwrap();
    let bath = st.bath.unwrap();
    assert_eq!((bath.unprojected.len(), bath.reaction.len()), (2, 2));
    for i in 0..2 {
        vec_near(bath.unprojected[i], raw[i]);
        vec_near(bath.reaction[i], [-on_j, on_j][i]);
    }
    let raw_moment = outer(q[0], raw[0]).plus(outer(q[1], raw[1]));
    let reaction_moment = outer(r, on_j);
    tensor_near(bath.unprojected_moment, raw_moment);
    tensor_near(bath.reaction_moment, reaction_moment);
    assert_eq!(st.stages[3].stage, BaoabO);
    tensor_near(st.stages[3].tensor_impulse, reaction_moment);
    assert!(reaction_moment.trace().abs() > 1e-4 && raw_moment.trace().abs() > 1e-4);
    let mechanical = complete(&st.interval.total).scaled(h * 8000.);
    let delta = moment(&s).plus(before.scaled(-1.));
    tensor_near(delta, mechanical.plus(raw_moment).plus(reaction_moment));
    println!(
        "delta={:?} mechanical={:?} raw={:?} reaction={:?}",
        delta.0, mechanical.0, raw_moment.0, reaction_moment.0
    );
}

#[test]
fn legacy_checkpoints() {
    for gamma in [0., 1.] {
        let (mut s, c) = rotor(1., 3.);
        let mut measured = s.clone();
        let mut b = LangevinBAOAB::with_constraints(gamma, 1., 1001, c.clone());
        let mut reporting = LangevinBAOAB::with_constraints(gamma, 1., 1001, c);
        let energies = b
            .step(&mut s, &ZeroModel, &mut NeighborList::new(1.), H)
            .unwrap();
        let out = reporting
            .step_with_stress(
                &mut measured,
                &ZeroModel,
                &mut NeighborList::new(1.),
                H,
                &ModelEvaluation::zero(),
            )
            .unwrap();
        assert_eq!(s.positions, measured.positions);
        assert_eq!(s.velocities, measured.velocities);
        assert_eq!(s.forces, measured.forces);
        assert_eq!(
            serde_json::to_value(energies).unwrap(),
            serde_json::to_value(out.energies).unwrap()
        );
        assert_eq!(out.stress.unwrap().solver, ConstraintSolverOptions::DEFAULT);
        let state = b.checkpoint_state().unwrap();
        assert_eq!(state["rng"]["word_pos"], "12");
        assert_eq!(state, reporting.checkpoint_state().unwrap());
        println!(
            "legacy {}",
            serde_json::json!([gamma, s.positions, s.velocities, state])
        );
    }
}

/// Conservative, anisotropic pair potential U = (r.x)^2 / 2.
struct Pair(AtomicUsize);
impl ForceModel for Pair {
    fn compute(&self, _: &mut System, _: &mut Virial, _: &mut NeighborList) -> Result<f64> {
        panic!("report must call evaluate once")
    }
    fn evaluate(&self, s: &mut System, _: &mut NeighborList) -> Result<ModelEvaluation> {
        self.0.fetch_add(1, Ordering::SeqCst);
        let r = s.pos(1) - s.pos(0);
        let f = Vec3::new(-r.0[0], 0., 0.);
        s.set_force(0, -f);
        s.set_force(1, f);
        let tensor = outer(r, f);
        Ok(ModelEvaluation {
            energy: r.0[0].powi(2) * 0.5,
            scalar_virial: tensor.trace(),
            tensor: TensorAvailability::Complete(tensor),
        })
    }
}
#[test]
fn endpoint_pair_quadrature_uses_cached_initial_evaluation() {
    let (mut s, cons) = rotor(1., 3.);
    let model = Pair(AtomicUsize::new(0));
    let mut nl = NeighborList::new(1.);
    let initial = model.evaluate(&mut s, &mut nl).unwrap();
    let before = moment(&s);
    let st = LangevinBAOAB::with_constraints(0., 1., 1001, cons)
        .step_with_stress_options(
            &mut s,
            &model,
            &mut nl,
            H,
            &initial,
            ConstraintSolverOptions::TIGHT,
        )
        .unwrap()
        .stress
        .unwrap();
    assert_eq!(model.0.load(Ordering::SeqCst), 2);
    let r = s.pos(1) - s.pos(0);
    let f = Vec3::new(-r.0[0], 0., 0.);
    assert_eq!(s.force(1), f);
    tensor_near(
        complete(&st.interval.model_average),
        complete(&initial.tensor).plus(outer(r, f)).scaled(0.5),
    );
    assert!(
        st.stages[0].tensor_impulse.trace().abs() > 1e-5
            && st.stages[6].tensor_impulse.trace().abs() > 1e-5
    );
    tensor_near(
        moment(&s).plus(before.scaled(-1.)),
        complete(&st.interval.total).scaled(H * 8000.),
    );
}

#[test]
fn invalid_timestep_or_options_precedes_state_and_rng_mutation() {
    let (mut s, c) = rotor(1., 3.);
    let mut b = LangevinBAOAB::with_constraints(1., 1., 1001, c);
    let state = b.checkpoint_state().unwrap();
    let before = s.clone();
    for (dt, solver) in [
        (0., ConstraintSolverOptions::TIGHT),
        (
            H,
            ConstraintSolverOptions {
                tolerance: 0.,
                max_iterations: 1,
            },
        ),
    ] {
        assert!(b
            .step_with_stress_options(
                &mut s,
                &ZeroModel,
                &mut NeighborList::new(1.),
                dt,
                &ModelEvaluation::zero(),
                solver
            )
            .is_err());
        assert_eq!(
            (s.positions.as_slice(), s.velocities.as_slice()),
            (before.positions.as_slice(), before.velocities.as_slice())
        );
        assert_eq!(b.checkpoint_state().unwrap(), state);
    }
}

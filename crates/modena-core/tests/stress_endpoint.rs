#[allow(dead_code)]
#[path = "support/stress.rs"]
mod support;

use modena_core::constraint_stress::{constraint_geometry, endpoint_reaction};
use modena_core::stress::{endpoint_pressure, ModelEvaluation, ValidationStatus};
use support::{near, rotor};

#[test]
fn endpoint_rotor_cancels_internal_scalar_pressure() {
    let (s, c) = rotor(1., 3.);
    let geometry = constraint_geometry(&s, &c).unwrap();
    let cached = s.forces.clone();
    let reaction = endpoint_reaction(&s, &geometry).unwrap();
    near(reaction.multipliers[0], 0.3675, 2e-9);
    near(reaction.tensor.trace(), -0.621075, 2e-9);
    let p = endpoint_pressure(&s, &ModelEvaluation::zero(), Some(&c)).unwrap();
    near(p.scalar.unwrap() * 3. * 8000., 0.56, 2e-9);
    assert_eq!(p.validation, ValidationStatus::Diagnostic);
    assert_eq!(s.forces, cached);
}

use modena_core::constraint_stress::EndpointReaction;
use modena_core::stress::{
    kinetic_tensor, MissingTensor, MissingTensorReason, PressureObservation, SymTensor6,
    TensorAvailability,
};
use modena_core::{Constraint, Constraints, SimBox, System, Vec3};
use support::{complete, water, water_rotation};

const CENTER: Vec3 = Vec3::new(5., 6., 7.);
const COM_V: Vec3 = Vec3::new(0.2, -0.1, 0.3);
const COMPONENTS: [(usize, usize); 6] = [(0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2)];

/// Independent Cartesian force moments and kinetic numerators, using the
/// analytic unwrapped coordinates (including for face-split fixtures).
fn moment(positions: &[Vec3], forces: &[Vec3]) -> [f64; 6] {
    COMPONENTS.map(|(a, b)| {
        positions
            .iter()
            .zip(forces)
            .map(|(q, f)| q.0[a] * f.0[b])
            .sum()
    })
}
fn check_case(
    s: &System,
    c: &Constraints,
    model: &ModelEvaluation,
    lambda: &[f64],
    forces: &[Vec3],
    virial: [f64; 6],
) {
    let g = constraint_geometry(s, c).unwrap();
    let mut position = 0_f64;
    let mut tangent = 0_f64;
    for e in &g.edges {
        let r = e.reference;
        position = position.max((r.norm2() - e.distance * e.distance).abs());
        tangent = tangent.max(r.dot(s.vel(e.j) - s.vel(e.i)).abs());
    }
    assert!(position <= 1e-13 && tangent <= 1e-13);
    let cached = s.forces.clone();
    let start = std::time::Instant::now();
    let reaction = endpoint_reaction(s, &g).unwrap();
    let solve_wall_ns = start.elapsed().as_nanos();
    let p = endpoint_pressure(s, model, Some(c)).unwrap();
    assert_eq!(s.forces, cached);
    assert_eq!(p.validation, ValidationStatus::Diagnostic);
    assert_eq!(reaction.multipliers.len(), lambda.len());
    assert_eq!(reaction.forces.len(), forces.len());
    assert!(reaction.relative_linear_residual <= 1e-11);
    for (actual, expected) in reaction.multipliers.iter().zip(lambda) {
        near(*actual, *expected, 2e-9);
    }
    for (actual, expected) in reaction.forces.iter().zip(forces) {
        for k in 0..3 {
            near(actual.0[k], expected.0[k], 2e-9);
        }
    }
    let kinetic = COMPONENTS.map(|(a, b)| {
        (0..s.n())
            .map(|i| s.masses[i] * s.vel(i).0[a] * s.vel(i).0[b])
            .sum::<f64>()
    });
    let expected_pressure: [f64; 6] =
        std::array::from_fn(|k| kinetic[k] + complete(&model.tensor).0[k] + virial[k]);
    for k in 0..6 {
        near(reaction.tensor.0[k], virial[k], 2e-9);
        near(p.constraint.0[k], virial[k], 2e-9);
        near(p.kinetic.0[k], kinetic[k], 2e-9);
        near(complete(&p.total).0[k] * 8000., expected_pressure[k], 2e-9);
    }
    near(
        p.scalar.unwrap() * 3. * 8000.,
        expected_pressure[..3].iter().sum(),
        2e-9,
    );
    if std::env::var_os("B07_ENDPOINT_EVIDENCE").is_some() {
        println!(
            "endpoint={}",
            serde_json::json!({"masses":s.masses,"positions":s.positions,"velocities":s.velocities,
            "conservative":s.forces,"edge_count":g.edges.len(),"expected_lambda":lambda,"expected_forces":forces,
            "expected_constraint":virial,"position_residual":position,"tangent_residual":tangent,"reaction":reaction,
            "observation":p,"endpoint_reaction_wall_ns":solve_wall_ns})
        );
    }
}

fn split_at_periodic_boundary(original: &System, mut check: impl FnMut(&System)) {
    check(original);
    let translated = Vec3::new(0.31, -0.27, 0.19);
    for x in [-1, 0, 1] {
        for y in [-1, 0, 1] {
            for z in [-1, 0, 1] {
                let indices = [x, y, z];
                let center = Vec3(std::array::from_fn(|a| match indices[a] {
                    -1 => 0.,
                    1 => 20.,
                    _ => CENTER.0[a],
                }));
                let mut s = original.clone();
                for i in 0..s.n() {
                    let q = original.pos(i) - CENTER + center;
                    s.set_pos(i, Vec3(q.0.map(|v| v.rem_euclid(20.))));
                }
                for (a, &index) in indices.iter().enumerate() {
                    if index != 0 {
                        assert!((0..s.n()).any(|i| s.pos(i).0[a] < 1.));
                        assert!((0..s.n()).any(|i| s.pos(i).0[a] > 19.));
                    }
                }
                check(&s);
                for i in 0..s.n() {
                    s.set_pos(
                        i,
                        Vec3((s.pos(i) + translated).0.map(|v| v.rem_euclid(20.))),
                    );
                }
                check(&s);
            }
        }
    }
}

#[test]
fn free_rotor_across_faces_and_corners() {
    let r = Vec3::new(1., 2., 2.) * (1.3 / 3.);
    for (m1, m2) in [(1., 3.), (2., 5.)] {
        let (s, c) = rotor(m1, m2);
        let lambda = m1 * m2 / (m1 + m2) * 0.7 * 0.7;
        let f = r * lambda;
        let virial = moment(&[r], &[-f]);
        split_at_periodic_boundary(&s, |state| {
            check_case(
                state,
                &c,
                &ModelEvaluation::zero(),
                &[lambda],
                &[f, -f],
                virial,
            )
        });
    }
}

#[test]
fn harmonic_trap_changes_pair_tension_at_two_orientations() {
    for (m1, m2) in [(1., 3.), (2., 5.)] {
        for rotated in [false, true] {
            let (mut s, c) = rotor(m1, m2);
            let mut rho = [Vec3::zero(); 2];
            let mut external = [Vec3::zero(); 2];
            for i in 0..2 {
                rho[i] = s.pos(i) - CENTER;
                if rotated {
                    rho[i] = water_rotation(rho[i]);
                    s.set_pos(i, CENTER + rho[i]);
                    s.set_vel(i, COM_V + water_rotation(s.vel(i) - COM_V));
                }
                external[i] = rho[i] * -0.1;
                s.set_force(i, external[i]);
            }
            let reduced = m1 * m2 / (m1 + m2);
            // The radial external relative acceleration is inward, magnitude
            // k*d*(m2/m1 + m1/m2)/M. It supplies part of the centripetal force.
            let lambda = reduced * (0.7 * 0.7 - 0.1 * (m2 / m1 + m1 / m2) / (m1 + m2));
            assert!((lambda - reduced * 0.7 * 0.7).abs() > 0.01);
            let f = (rho[1] - rho[0]) * lambda;
            let external_moment = SymTensor6(moment(&rho, &external));
            let model = ModelEvaluation {
                energy: 0.05 * (rho[0].norm2() + rho[1].norm2()),
                scalar_virial: external_moment.trace(),
                tensor: TensorAvailability::Complete(external_moment),
            };
            check_case(&s, &c, &model, &[lambda], &[f, -f], moment(&rho, &[f, -f]));
        }
    }
}

/// Invert the three Cartesian inertia columns by scalar triple products.
/// This is Euler's rigid-body balance, independent of edge gradients/masses.
fn euler_alpha(rho: &[Vec3; 3], mass: &[f64; 3], omega: Vec3) -> Vec3 {
    let inertia = |v: Vec3| {
        rho.iter()
            .zip(mass)
            .fold(Vec3::zero(), |sum, (r, &m)| sum + r.cross(v.cross(*r)) * m)
    };
    let torque = -omega.cross(inertia(omega));
    let [a, b, c] = [
        Vec3::new(1., 0., 0.),
        Vec3::new(0., 1., 0.),
        Vec3::new(0., 0., 1.),
    ]
    .map(inertia);
    let determinant = a.dot(b.cross(c));
    assert!(determinant > 0.);
    Vec3::new(
        torque.dot(b.cross(c)),
        torque.dot(c.cross(a)),
        torque.dot(a.cross(b)),
    ) * (1. / determinant)
}

#[test]
fn euler_water_reactions_principal_and_nonprincipal() {
    for masses in [[15.9994, 1.008, 1.008], [16., 1., 2.]] {
        for principal in [true, false] {
            let (mut s, c) = water(masses);
            let rho: [Vec3; 3] = std::array::from_fn(|i| s.pos(i) - CENTER);
            let omega = if principal {
                water_rotation(Vec3::new(0., 0., 0.7))
            } else {
                Vec3::new(0.2, 0.3, 0.4)
            };
            let alpha = euler_alpha(&rho, &masses, omega);
            if principal {
                assert!(alpha.norm() < 1e-13);
            } else {
                assert!(alpha.norm() > 0.01);
            }
            let forces: [Vec3; 3] = std::array::from_fn(|i| {
                s.set_vel(i, COM_V + omega.cross(rho[i]));
                if principal {
                    rho[i] * (-0.7 * 0.7 * masses[i])
                } else {
                    (alpha.cross(rho[i]) + omega.cross(omega.cross(rho[i]))) * masses[i]
                }
            });
            // At O, F_O = λ01*r01 + λ02*r02. Take oriented planar
            // areas to isolate each coefficient. H1 then fixes λ12.
            let a = rho[1] - rho[0];
            let b = rho[2] - rho[0];
            let hh = b - a;
            let normal = a.cross(b);
            let l01 = forces[0].cross(b).dot(normal) / normal.norm2();
            let l02 = a.cross(forces[0]).dot(normal) / normal.norm2();
            let l12 = (forces[1] + a * l01).dot(hh) / hh.norm2();
            let virial = moment(&rho, &forces);
            let internal_ke: f64 = (0..3)
                .map(|i| masses[i] * omega.cross(rho[i]).norm2())
                .sum();
            near(internal_ke + virial[..3].iter().sum::<f64>(), 0., 2e-9);
            split_at_periodic_boundary(&s, |state| {
                check_case(
                    state,
                    &c,
                    &ModelEvaluation::zero(),
                    &[l01, l02, l12],
                    &forces,
                    virial,
                )
            });
        }
    }
}

#[test]
fn empty_constraints_are_exact_zero_and_missing_model_stays_missing() {
    let (s, c) = rotor(2., 5.);
    let empty = Constraints::new();
    let reaction = endpoint_reaction(&s, &constraint_geometry(&s, &empty).unwrap()).unwrap();
    assert!(reaction.multipliers.is_empty());
    assert_eq!(reaction.forces, vec![Vec3::zero(); 2]);
    assert_eq!(reaction.tensor, SymTensor6::ZERO);
    assert_eq!(reaction.relative_linear_residual, 0.);
    let p = endpoint_pressure(&s, &ModelEvaluation::zero(), None).unwrap();
    assert_eq!(p.constraint, SymTensor6::ZERO);
    assert_eq!(p.kinetic, kinetic_tensor(&s));
    assert_eq!(
        p,
        endpoint_pressure(&s, &ModelEvaluation::zero(), Some(&empty)).unwrap()
    );
    let missing = TensorAvailability::Missing(vec![MissingTensor {
        component: "legacy".into(),
        reason: MissingTensorReason::LegacyForceModel,
    }]);
    let model = ModelEvaluation {
        energy: 4.,
        scalar_virial: 123.,
        tensor: missing.clone(),
    };
    let p = endpoint_pressure(&s, &model, Some(&c)).unwrap();
    assert_eq!(
        (p.model, p.total, p.scalar),
        (missing.clone(), missing, None)
    );
    assert_ne!(p.constraint, SymTensor6::ZERO);
}

#[test]
fn independent_connected_blocks_preserve_global_edge_order() {
    let (a, _) = rotor(1., 3.);
    let (b, _) = rotor(2., 5.);
    let mut s = System::new(SimBox::cubic(20.));
    s.resize(4);
    for (offset, other) in [(0, a), (2, b)] {
        for i in 0..2 {
            s.set_pos(offset + i, other.pos(i) + Vec3::splat(offset as f64));
            s.set_vel(offset + i, other.vel(i));
            s.masses[offset + i] = other.masses[i];
        }
    }
    let c = Constraints {
        bonds: vec![
            Constraint { i: 3, j: 2, d: 1.3 },
            Constraint { i: 0, j: 1, d: 1.3 },
        ],
        waters: vec![],
    };
    let g = constraint_geometry(&s, &c).unwrap();
    let r = endpoint_reaction(&s, &g).unwrap();
    near(r.multipliers[0], 10. / 7. * 0.49, 2e-9);
    near(r.multipliers[1], 0.3675, 2e-9);
    for i in 0..4 {
        s.set_vel(i, COM_V);
    }
    let static_result = endpoint_reaction(&s, &g).unwrap();
    assert_eq!(static_result.relative_linear_residual, 0.);
    assert_eq!(static_result.forces, vec![Vec3::zero(); 4]);
}

#[test]
fn invalid_endpoint_states_and_overflows_return_physics_errors() {
    let (s, c) = rotor(1., 3.);
    let g = constraint_geometry(&s, &c).unwrap();
    let zero = ModelEvaluation::zero();
    let physics = |err| assert!(matches!(err, modena_core::ModenaError::Physics(_)), "{err}");
    let mut bad = s.clone();
    bad.add_vel(1, Vec3::splat(1e-6));
    let err = endpoint_reaction(&bad, &g).unwrap_err();
    assert!(err.to_string().contains("tangent residual"));
    physics(err);
    bad = s.clone();
    bad.add_pos(1, Vec3::splat(1e-6));
    physics(endpoint_reaction(&bad, &g).unwrap_err());
    let err = endpoint_pressure(&bad, &zero, Some(&c)).unwrap_err();
    assert!(err.to_string().contains("position residual"));
    physics(err);
    let mut forged = g.clone();
    forged.edges[0].image_shift[1] = 1;
    physics(endpoint_reaction(&s, &forged).unwrap_err());
    forged = g.clone();
    forged.edges.push(g.edges[0].clone());
    physics(endpoint_reaction(&s, &forged).unwrap_err());
    forged = g.clone();
    forged.edges[0].j = usize::MAX;
    physics(endpoint_reaction(&s, &forged).unwrap_err());
    bad = s.clone();
    bad.forces.pop();
    physics(endpoint_reaction(&bad, &g).unwrap_err());
    bad = s.clone();
    bad.forces[0] = f64::NAN;
    physics(endpoint_pressure(&bad, &zero, Some(&c)).unwrap_err());
    bad = s.clone();
    bad.masses[0] = 0.;
    physics(endpoint_pressure(&bad, &zero, None).unwrap_err());
    bad = s.clone();
    bad.box_size = SimBox::cubic(0.);
    physics(endpoint_pressure(&bad, &zero, None).unwrap_err());
    bad = s.clone();
    bad.set_vel(0, Vec3::splat(f64::MAX));
    physics(endpoint_pressure(&bad, &zero, None).unwrap_err());
    let mut model = zero;
    model.tensor = TensorAvailability::Complete(SymTensor6([f64::NAN; 6]));
    physics(endpoint_pressure(&s, &model, None).unwrap_err());

    let line = support::atoms(&[[1., 1., 1.], [2., 1., 1.], [3., 1., 1.]], &[1.; 3]);
    let redundant = Constraints {
        bonds: vec![
            Constraint { i: 0, j: 1, d: 1. },
            Constraint { i: 1, j: 2, d: 1. },
            Constraint { i: 0, j: 2, d: 2. },
        ],
        waters: vec![],
    };
    let err = endpoint_pressure(&line, &ModelEvaluation::zero(), Some(&redundant)).unwrap_err();
    assert!(err.to_string().contains("scaled pivot"));
    physics(err);
}

#[test]
fn serde_rejects_invalid_observations_and_deserialized_claims_are_diagnostic() {
    let (s, c) = rotor(1., 3.);
    let p = endpoint_pressure(&s, &ModelEvaluation::zero(), Some(&c)).unwrap();
    let claim = ValidationStatus::Validated {
        pub_record: "claimed".into(),
    };
    let mut wire = serde_json::to_value(&p).unwrap();
    assert_eq!(
        serde_json::to_value(&claim).unwrap(),
        serde_json::json!({"kind":"validated","value":{"pub_record":"claimed"}})
    );
    wire["validation"] = serde_json::to_value(claim).unwrap();
    assert_eq!(
        serde_json::from_value::<PressureObservation>(wire.clone()).unwrap(),
        p
    );
    let mut bad = wire.clone();
    bad["extra"] = 0.into();
    assert!(serde_json::from_value::<PressureObservation>(bad).is_err());
    let mut bad = wire.clone();
    bad["volume"] = 0.into();
    assert!(serde_json::from_value::<PressureObservation>(bad).is_err());
    wire["scalar"] = serde_json::Value::Null;
    assert!(serde_json::from_value::<PressureObservation>(wire).is_err());
    let reaction = endpoint_reaction(&s, &constraint_geometry(&s, &c).unwrap()).unwrap();
    let mut value = toml::Value::try_from(&reaction).unwrap();
    for (key, bad) in [
        (
            "multipliers",
            toml::Value::Array(vec![toml::Value::Float(f64::INFINITY)]),
        ),
        ("relative_linear_residual", toml::Value::Float(f64::NAN)),
    ] {
        let mut invalid = value.clone();
        invalid[key] = bad;
        assert!(invalid.try_into::<EndpointReaction>().is_err());
    }
    assert_eq!(
        value.clone().try_into::<EndpointReaction>().unwrap(),
        reaction
    );
    value
        .as_table_mut()
        .unwrap()
        .insert("extra".into(), 0.into());
    assert!(value.try_into::<EndpointReaction>().is_err());
}

#[allow(dead_code)]
#[path = "support/stress.rs"]
mod support;
use modena_core::constraint_stress::{constraint_geometry, ProjectionStage};
use modena_core::Vec3;
use support::{near, rotor};

#[test]
fn rattle_report_reconstructs_actual_momentum_change() {
    let (mut s, c) = rotor(1., 3.);
    let g = constraint_geometry(&s, &c).unwrap();
    s.add_vel(1, g.edges[0].reference * 0.2);
    let before = s.velocities.clone();
    let rep = c
        .rattle_report(&mut s, &g, ProjectionStage::VvFinal, 1e-13, 100)
        .unwrap();
    let mut dp = vec![Vec3::zero(); s.n()];
    for j in &rep.impulses {
        let e = &g.edges[j.edge];
        dp[e.j] += j.impulse_on_j;
        dp[e.i] -= j.impulse_on_j;
    }
    for i in 0..s.n() {
        for a in 0..3 {
            near(
                dp[i].0[a],
                s.masses[i] * (s.velocities[3 * i + a] - before[3 * i + a]),
                5e-10,
            );
        }
    }
    near(rep.tensor_impulse.trace(), -0.2535, 5e-10);
    check_ledger(&s, &g, &rep, &before, None);
}

#[test]
fn position_report_reconstructs_mass_correction_divided_by_interval() {
    let (mut s, c) = rotor(2., 5.);
    let g = constraint_geometry(&s, &c).unwrap();
    let interval = 0.08 / 0.7;
    for i in 0..s.n() {
        s.add_pos(i, s.vel(i) * interval);
    }
    let predicted = s.positions.clone();
    let velocities = s.velocities.clone();
    let rep = c
        .shake_report(
            &mut s,
            &g,
            interval,
            ProjectionStage::VvPosition,
            1e-13,
            100,
        )
        .unwrap();
    let mut dp = vec![Vec3::zero(); s.n()];
    for j in &rep.impulses {
        let e = &g.edges[j.edge];
        dp[e.j] += j.impulse_on_j;
        dp[e.i] -= j.impulse_on_j;
    }
    for i in 0..s.n() {
        for a in 0..3 {
            near(
                dp[i].0[a],
                s.masses[i] * (s.positions[3 * i + a] - predicted[3 * i + a]) / interval,
                5e-10,
            );
        }
    }
    assert_eq!(s.velocities, velocities);
    let mu = 10. / 7.;
    let cosine = (1_f64 - 0.08 * 0.08).sqrt();
    near(
        rep.tensor_impulse.trace(),
        mu * 1.3 * 1.3 * (cosine - 1.) / interval,
        5e-10,
    );
}

use modena_core::constraint_stress::{ConstraintGeometry, ProjectionReport};
use modena_core::constraints::{Constraint, Constraints};
use modena_core::{SimBox, System};
use support::water;

// Independent full (unsymmetrized) moments catch errors hidden by compression.
fn check_ledger(
    s: &System,
    g: &ConstraintGeometry,
    rep: &ProjectionReport,
    before: &[f64],
    interval: Option<f64>,
) {
    let mut impulses = vec![Vec3::zero(); s.n()];
    let mut moment = [[0.; 3]; 3];
    let mut torque = Vec3::zero();
    for step in &rep.impulses {
        let e = &g.edges[step.edge];
        impulses[e.j] += step.impulse_on_j;
        impulses[e.i] -= step.impulse_on_j;
        torque += step.lever.cross(step.impulse_on_j);
        for (a, row) in moment.iter_mut().enumerate() {
            for (b, value) in row.iter_mut().enumerate() {
                *value += step.lever.0[a] * step.impulse_on_j.0[b];
            }
        }
    }
    let mut total = Vec3::zero();
    for (i, impulse) in impulses.iter().enumerate() {
        total += *impulse;
        for a in 0..3 {
            let after = if interval.is_some() {
                s.positions[3 * i + a]
            } else {
                s.velocities[3 * i + a]
            };
            near(
                impulse.0[a],
                s.masses[i] * (after - before[3 * i + a]) / interval.unwrap_or(1.),
                5e-10,
            );
        }
    }
    for x in total.0.into_iter().chain(torque.0) {
        near(x, 0., 5e-10);
    }
    for (k, (a, b)) in [(0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2)]
        .into_iter()
        .enumerate()
    {
        near(
            rep.tensor_impulse.0[k],
            0.5 * (moment[a][b] + moment[b][a]),
            5e-10,
        );
        near(moment[a][b], moment[b][a], 5e-10);
    }
    // Unwrap the small fixture independently from its anchored reference edges.
    let mut q = vec![Vec3::zero(); s.n()];
    for e in &g.edges {
        if e.i == 0 {
            q[e.j] = e.reference;
        }
    }
    let mut atom_moment = [[0.; 3]; 3];
    for i in 0..s.n() {
        for (a, row) in atom_moment.iter_mut().enumerate() {
            for (b, x) in row.iter_mut().enumerate() {
                *x += q[i].0[a] * impulses[i].0[b];
            }
        }
    }
    for a in 0..3 {
        for b in 0..3 {
            near(atom_moment[a][b], moment[a][b], 5e-10);
        }
    }
    if std::env::var_os("B07_IMPULSE_EVIDENCE").is_some() {
        println!(
            "ledger={}",
            serde_json::json!({"masses":s.masses,"edges":g.edges,"report":rep,"per_atom_impulses":impulses,"full_moment":moment,"torque":torque,"momentum":total,"interval":interval,"entry_bytes":std::mem::size_of::<modena_core::constraint_stress::ConstraintImpulse>(),"capacity":rep.impulses.capacity(),"rollback_bytes":std::mem::size_of_val(before)})
        );
    }
    println!("stage={:?} iterations={} entries={} tensor={:?} torque={:?} position_residual={:.17e} velocity_residual={:.17e}",rep.stage,rep.iterations,rep.impulses.len(),rep.tensor_impulse.0,torque.0,rep.max_position_residual,rep.max_velocity_residual);
}

#[test]
fn rotor_records_every_actual_newton_increment() {
    let (mut s, c) = rotor(1., 3.);
    let g = constraint_geometry(&s, &c).unwrap();
    let h = 0.08 / 0.7;
    for i in 0..s.n() {
        s.add_pos(i, s.vel(i) * h);
    }
    let before = s.positions.clone();
    let mut replay = s.clone();
    let rep = c
        .shake_report(&mut s, &g, h, ProjectionStage::VvPosition, 1e-13, 100)
        .unwrap();
    assert!(rep.impulses.len() > 1);
    assert_eq!(rep.iterations, rep.impulses.len());
    for impulse in &rep.impulses {
        let e = &g.edges[impulse.edge];
        let r = replay.pos(e.j) - replay.pos(e.i);
        let expected = -(r.norm2() - 1.3 * 1.3) / (2. * (1. + 1. / 3.) * r.dot(e.reference));
        for a in 0..3 {
            near(
                impulse.impulse_on_j.0[a],
                e.reference.0[a] * expected / h,
                5e-10,
            );
        }
        replay.add_pos(e.i, impulse.impulse_on_j * (-h / replay.masses[e.i]));
        replay.add_pos(e.j, impulse.impulse_on_j * (h / replay.masses[e.j]));
    }
    for (a, b) in s.positions.iter().zip(&replay.positions) {
        near(*a, *b, 5e-10);
    }
    check_ledger(&s, &g, &rep, &before, Some(h));
}

#[test]
fn water_newton_and_three_bonds_agree_for_both_mass_sets() {
    for masses in [[15.9994, 1.008, 1.008], [16., 1., 2.]] {
        let (initial, water_c) = water(masses);
        let water_g = constraint_geometry(&initial, &water_c).unwrap();
        let bonds = Constraints {
            bonds: water_g
                .edges
                .iter()
                .map(|e| Constraint {
                    i: e.i,
                    j: e.j,
                    d: e.distance,
                })
                .collect(),
            waters: vec![],
        };
        let mut water_s = initial.clone();
        let mut bond_s = initial.clone();
        let h = 0.08 / 0.7;
        for s in [&mut water_s, &mut bond_s] {
            for i in 0..s.n() {
                s.add_pos(i, s.vel(i) * h);
            }
        }
        let before = water_s.positions.clone();
        let a = water_c
            .shake_report(
                &mut water_s,
                &water_g,
                h,
                ProjectionStage::VvPosition,
                1e-13,
                100,
            )
            .unwrap();
        let b = bonds
            .shake_report(
                &mut bond_s,
                &water_g,
                h,
                ProjectionStage::VvPosition,
                1e-13,
                100,
            )
            .unwrap();
        assert_eq!(a.impulses.len(), 3 * a.iterations);
        assert!(a.iterations > 1);
        check_ledger(&water_s, &water_g, &a, &before, Some(h));
        check_ledger(&bond_s, &water_g, &b, &before, Some(h));
        assert!(a.max_position_residual <= 1e-13 && b.max_position_residual <= 1e-13);
        // Replay every simultaneous Newton triple, not only its net multiplier.
        let mut replay = initial.clone();
        replay.positions = before.clone();
        let mut prior = f64::INFINITY;
        for pass in a.impulses.as_chunks::<3>().0 {
            let residual = water_c.max_residual(&replay);
            assert!(residual < prior);
            prior = residual;
            for impulse in pass {
                let e = &water_g.edges[impulse.edge];
                replay.add_pos(e.i, impulse.impulse_on_j * (-h / replay.masses[e.i]));
                replay.add_pos(e.j, impulse.impulse_on_j * (h / replay.masses[e.j]));
            }
        }
        for ((x, y), z) in water_s
            .positions
            .iter()
            .zip(&bond_s.positions)
            .zip(&replay.positions)
        {
            near(*x, *y, 5e-10);
            near(*x, *z, 5e-10);
        }
        for (x, y) in a.tensor_impulse.0.into_iter().zip(b.tensor_impulse.0) {
            near(x, y, 5e-10);
        }
        for s in [&mut water_s, &mut bond_s] {
            for (i, x) in before.iter().enumerate() {
                s.velocities[i] += (s.positions[i] - x) / h;
            }
        }
        let vg = constraint_geometry(&water_s, &water_c).unwrap();
        let vb = water_s.velocities.clone();
        let a = water_c
            .rattle_report(&mut water_s, &vg, ProjectionStage::VvFinal, 1e-13, 100)
            .unwrap();
        let bg = constraint_geometry(&bond_s, &bonds).unwrap();
        let bb = bond_s.velocities.clone();
        let b = bonds
            .rattle_report(&mut bond_s, &bg, ProjectionStage::VvFinal, 1e-13, 100)
            .unwrap();
        assert_eq!(a.iterations, 1);
        assert_eq!(a.impulses.len(), 3);
        check_ledger(&water_s, &vg, &a, &vb, None);
        check_ledger(&bond_s, &bg, &b, &bb, None);
        assert!(a.max_velocity_residual <= 1e-13 && b.max_velocity_residual <= 1e-13);
        for (x, y) in water_s.velocities.iter().zip(&bond_s.velocities) {
            near(*x, *y, 5e-10);
        }
        for (x, y) in a.tensor_impulse.0.into_iter().zip(b.tensor_impulse.0) {
            near(x, y, 5e-10);
        }
    }
}

#[test]
fn reports_are_invariant_across_faces_corners_and_rigid_translation() {
    for (initial, c) in [rotor(1., 3.), rotor(2., 5.), water([16., 1., 2.])] {
        let mut expected = None;
        let mut expected_velocity = None;
        for shift in [
            Vec3::zero(),
            Vec3::new(-5., 0., 0.),
            Vec3::new(0., -6., 0.),
            Vec3::new(0., 0., -7.),
            Vec3::new(-5., -6., -7.),
            Vec3::new(0.31, -0.27, 0.19),
        ] {
            let mut s = initial.clone();
            for i in 0..s.n() {
                s.add_pos(i, shift);
            }
            s.wrap_positions();
            let g = constraint_geometry(&s, &c).unwrap();
            if g.edges.len() == 3 {
                assert!(
                    (g.edges[1].reference - g.edges[0].reference - g.edges[2].reference).norm()
                        < 1e-12
                );
            }
            let h = 0.08 / 0.7;
            for i in 0..s.n() {
                s.add_pos(i, s.vel(i) * h);
            }
            let before = s.positions.clone();
            let velocities = s.velocities.clone();
            let report = c
                .shake_report(&mut s, &g, h, ProjectionStage::VvPosition, 1e-13, 100)
                .unwrap();
            check_ledger(&s, &g, &report, &before, Some(h));
            assert_eq!(velocities, s.velocities);
            if let Some(t) = expected {
                for (x, y) in report.tensor_impulse.0.into_iter().zip(t) {
                    near(x, y, 5e-10);
                }
            } else {
                expected = Some(report.tensor_impulse.0);
            }
            for (k, x) in before.iter().enumerate() {
                s.velocities[k] += (s.positions[k] - x) / h;
            }
            s.wrap_positions();
            let stage = constraint_geometry(&s, &c).unwrap();
            let vb = s.velocities.clone();
            let post = c
                .rattle_report(&mut s, &stage, ProjectionStage::VvFinal, 1e-13, 100)
                .unwrap();
            check_ledger(&s, &stage, &post, &vb, None);
            if let Some(t) = expected_velocity {
                for (x, y) in post.tensor_impulse.0.into_iter().zip(t) {
                    near(x, y, 5e-10);
                }
            } else {
                expected_velocity = Some(post.tensor_impulse.0);
            }
            let mut translated = initial.clone();
            let tg = constraint_geometry(&translated, &c).unwrap();
            for i in 0..translated.n() {
                translated.add_pos(i, Vec3::new(45., -60., 80.) + translated.vel(i) * h);
            }
            let rep = c
                .shake_report(
                    &mut translated,
                    &tg,
                    h,
                    ProjectionStage::VvPosition,
                    1e-13,
                    100,
                )
                .unwrap();
            for (x, y) in rep
                .tensor_impulse
                .0
                .into_iter()
                .zip(report.tensor_impulse.0)
            {
                near(x, y, 5e-10);
            }
        }
    }
}

#[test]
fn rejects_ambiguous_duplicate_winding_and_singular_graphs() {
    let (mut s, mut c) = rotor(1., 3.);
    s.set_pos(0, Vec3::zero());
    s.set_pos(1, Vec3::new(10., 0., 0.));
    assert!(constraint_geometry(&s, &c).is_err());
    let (s, base) = rotor(1., 3.);
    c = base;
    c.bonds.push(Constraint { i: 1, j: 0, d: 1.3 });
    assert!(constraint_geometry(&s, &c).is_err());
    let mut s = support::atoms(&[[0., 0., 0.], [4., 0., 0.], [8., 0., 0.]], &[1., 1., 1.]);
    s.box_size = SimBox::cubic(10.);
    let c = Constraints {
        bonds: vec![
            Constraint { i: 0, j: 1, d: 4. },
            Constraint { i: 1, j: 2, d: 4. },
            Constraint { i: 0, j: 2, d: 2. },
        ],
        waters: vec![],
    };
    assert!(constraint_geometry(&s, &c)
        .unwrap_err()
        .to_string()
        .contains("winding"));
    s.set_pos(1, Vec3::new(1., 0., 0.));
    s.set_pos(2, Vec3::new(2., 0., 0.));
    assert!(constraint_geometry(&s, &c)
        .unwrap_err()
        .to_string()
        .contains("singular"));
    s.masses[1] = 0.;
    assert!(constraint_geometry(&s, &c).is_err());
}

#[test]
fn skew_cell_search_finds_far_shift_and_true_shortest_vector() {
    let mut s = support::atoms(&[[0., 0., 0.], [19.6, 1., 0.2]], &[1., 3.]);
    s.box_size = SimBox::triclinic([
        Vec3::new(10., 0., 0.),
        Vec3::new(39., 2., 0.),
        Vec3::new(0., 0., 10.),
    ]);
    let c = Constraints {
        bonds: vec![Constraint {
            i: 0,
            j: 1,
            d: 1.2_f64.sqrt(),
        }],
        waters: vec![],
    };
    let g = constraint_geometry(&s, &c).unwrap();
    assert_eq!(g.edges[0].image_shift, [-2, 0, 0]);
    for (x, y) in g.edges[0].reference.0.into_iter().zip([-0.4, 1., 0.2]) {
        near(x, y, 1e-13);
    }
    // Shortest lattice vector b-4a=(-1,2,0), not a basis vector.
    s.add_pos(1, Vec3::new(0., 1.2, 0.));
    let before = s.positions.clone();
    assert!(c
        .shake_report(&mut s, &g, 0.1, ProjectionStage::VvPosition, 1e-13, 100)
        .is_err());
    assert_eq!(s.positions, before);
    s.set_pos(1, Vec3::new(-0.5, 1., 0.));
    assert!(constraint_geometry(&s, &c).is_err());
}

#[test]
fn report_failures_restore_source_arrays_and_reject_stale_geometry() {
    let (initial, c) = rotor(1., 3.);
    let g = constraint_geometry(&initial, &c).unwrap();
    for (interval, iterations) in [
        (0.1, 1),
        (1e-320, 100),
        (0., 100),
        (f64::INFINITY, 100),
        (0.1, 0),
    ] {
        let mut s = initial.clone();
        for i in 0..s.n() {
            s.add_pos(i, s.vel(i) * 0.1);
        }
        let (positions, velocities) = (s.positions.clone(), s.velocities.clone());
        assert!(c
            .shake_report(
                &mut s,
                &g,
                interval,
                ProjectionStage::VvPosition,
                1e-13,
                iterations
            )
            .is_err());
        assert_eq!(s.positions, positions);
        assert_eq!(s.velocities, velocities);
    }
    let mut s = initial.clone();
    s.add_pos(1, Vec3::new(0.001, 0., 0.));
    let velocities = s.velocities.clone();
    assert!(c
        .rattle_report(&mut s, &g, ProjectionStage::VvFinal, 1e-13, 100)
        .is_err());
    assert_eq!(s.velocities, velocities);
    let (mut s, w) = water([16., 1., 2.]);
    let wg = constraint_geometry(&s, &w).unwrap();
    let c = Constraints {
        bonds: wg
            .edges
            .iter()
            .map(|e| Constraint {
                i: e.i,
                j: e.j,
                d: e.distance,
            })
            .collect(),
        waters: vec![],
    };
    s.add_vel(1, Vec3::new(0.4, 0.2, -0.3));
    let (positions, velocities) = (s.positions.clone(), s.velocities.clone());
    assert!(c
        .rattle_report(&mut s, &wg, ProjectionStage::VvFinal, 1e-13, 1)
        .is_err());
    assert_eq!(s.positions, positions);
    assert_eq!(s.velocities, velocities);
}

#[test]
fn wrappers_share_updates_without_a_ledger_and_settle_remains_repair() {
    for (mut s, c) in [rotor(1., 3.), water([16., 1., 2.])] {
        let reference = s.positions.clone();
        let g = constraint_geometry(&s, &c).unwrap();
        for i in 0..s.n() {
            s.add_pos(i, s.vel(i) * 0.1);
        }
        let mut legacy = s.clone();
        c.shake(&mut legacy, &reference, 1e-13, 100).unwrap();
        c.shake_report(&mut s, &g, 0.1, ProjectionStage::VvPosition, 1e-13, 100)
            .unwrap();
        assert_eq!(s.positions, legacy.positions);
        let g = constraint_geometry(&s, &c).unwrap();
        c.rattle(&mut legacy, 1e-13, 100).unwrap();
        c.rattle_report(&mut s, &g, ProjectionStage::VvFinal, 1e-13, 100)
            .unwrap();
        assert_eq!(s.velocities, legacy.velocities);
    }
    let (mut s, c) = water([16., 1., 2.]);
    s.add_pos(1, Vec3::new(0.15, -0.1, 0.07));
    let velocities = s.velocities.clone();
    c.settle(&mut s, 1e-13).unwrap();
    assert!(c.max_residual(&s) <= 1e-13);
    assert_eq!(s.velocities, velocities);
}

#[test]
fn report_serde_is_strict_and_vectors_are_transparent_finite_arrays() {
    use serde::Deserialize;
    assert_eq!(
        serde_json::to_string(&Vec3::new(1., 2., 3.)).unwrap(),
        "[1.0,2.0,3.0]"
    );
    let invalid = serde::de::value::SeqDeserializer::<_, serde::de::value::Error>::new(
        [f64::NAN, 0., 0.].into_iter(),
    );
    assert!(Vec3::deserialize(invalid).is_err());
    let (mut s, c) = rotor(1., 3.);
    let g = constraint_geometry(&s, &c).unwrap();
    s.add_vel(1, g.edges[0].reference * 0.2);
    let rep = c
        .rattle_report(&mut s, &g, ProjectionStage::VvFinal, 1e-13, 100)
        .unwrap();
    let mut value = serde_json::to_value(&rep).unwrap();
    assert_eq!(value["stage"]["kind"], "vv_final");
    assert_eq!(
        serde_json::from_value::<ProjectionReport>(value.clone()).unwrap(),
        rep
    );
    value["extra"] = serde_json::json!(1);
    assert!(serde_json::from_value::<ProjectionReport>(value).is_err());
    assert!(serde_json::from_str::<ProjectionStage>(r#"{"kind":"vv_final","extra":0}"#).is_err());
    let mut edge = serde_json::to_value(&g.edges[0]).unwrap();
    edge["extra"] = serde_json::json!(1);
    assert!(
        serde_json::from_value::<modena_core::constraint_stress::ConstraintEdge>(edge).is_err()
    );
}

#[test]
fn public_geometry_cannot_sneak_a_nonminimal_or_tied_reference_into_a_report() {
    let (mut s, c) = rotor(1., 3.);
    let mut g = constraint_geometry(&s, &c).unwrap();
    g.edges[0].image_shift[0] += 1;
    g.edges[0].reference += Vec3::new(20., 0., 0.);
    let positions = s.positions.clone();
    let velocities = s.velocities.clone();
    assert!(c
        .rattle_report(&mut s, &g, ProjectionStage::VvFinal, 1e-13, 100)
        .is_err());
    assert_eq!(s.positions, positions);
    assert_eq!(s.velocities, velocities);
    s.set_pos(0, Vec3::zero());
    s.set_pos(1, Vec3::new(10., 0., 0.));
    g.edges[0].reference = s.pos(1);
    g.edges[0].image_shift = [0; 3];
    let positions = s.positions.clone();
    assert!(c
        .shake_report(&mut s, &g, 0.1, ProjectionStage::VvPosition, 1e-13, 100)
        .is_err());
    assert_eq!(s.positions, positions);
    assert_eq!(s.velocities, velocities);
}

#[test]
fn velocity_overflow_after_an_earlier_correction_rolls_back_both_arrays() {
    let mut s = support::atoms(
        &[[0., 0., 0.], [1., 0., 0.], [4., 2., 1.], [7., 2., 1.]],
        &[1.; 4],
    );
    let c = Constraints {
        bonds: vec![
            Constraint { i: 0, j: 1, d: 1. },
            Constraint { i: 2, j: 3, d: 3. },
        ],
        waters: vec![],
    };
    let g = constraint_geometry(&s, &c).unwrap();
    s.set_vel(1, Vec3::new(0.2, 0., 0.));
    s.set_vel(3, Vec3::new(1e308, 0., 0.));
    let positions = s.positions.clone();
    let velocities = s.velocities.clone();
    assert!(matches!(
        c.rattle_report(&mut s, &g, ProjectionStage::VvFinal, 1e-13, 100),
        Err(modena_core::error::ModenaError::Physics(_))
    ));
    assert_eq!(s.positions, positions);
    assert_eq!(s.velocities, velocities);
}

#[test]
fn report_contract_shake_only_accepts_position_stages() {
    use ProjectionStage::*;
    for (stage, valid) in [
        (VvPosition, true),
        (VvFinal, false),
        (BaoabB0, false),
        (BaoabA1Pre, true),
        (BaoabA1Post, false),
        (BaoabO, false),
        (BaoabA2Pre, true),
        (BaoabA2Post, false),
        (BaoabB1, false),
    ] {
        let (mut s, c) = rotor(1., 3.);
        let g = constraint_geometry(&s, &c).unwrap();
        s.add_pos(1, g.edges[0].reference * 0.01);
        let before = s.clone();
        let result = c.shake_report(&mut s, &g, 0.1, stage, 1e-13, 100);
        if valid {
            assert_eq!(result.unwrap().stage, stage);
        } else {
            assert!(
                matches!(result, Err(modena_core::error::ModenaError::Physics(_))),
                "{stage:?}"
            );
            assert_eq!(s.positions, before.positions);
            assert_eq!(s.velocities, before.velocities);
        }
    }
}
#[test]
fn report_contract_rattle_only_accepts_velocity_stages() {
    use ProjectionStage::*;
    for (stage, valid) in [
        (VvPosition, false),
        (VvFinal, true),
        (BaoabB0, true),
        (BaoabA1Pre, false),
        (BaoabA1Post, true),
        (BaoabO, true),
        (BaoabA2Pre, false),
        (BaoabA2Post, true),
        (BaoabB1, true),
    ] {
        let (mut s, c) = rotor(1., 3.);
        let g = constraint_geometry(&s, &c).unwrap();
        s.add_vel(1, g.edges[0].reference * 0.2);
        let before = s.clone();
        let result = c.rattle_report(&mut s, &g, stage, 1e-13, 100);
        if valid {
            assert_eq!(result.unwrap().stage, stage);
        } else {
            assert!(
                matches!(result, Err(modena_core::error::ModenaError::Physics(_))),
                "{stage:?}"
            );
            assert_eq!(s.positions, before.positions);
            assert_eq!(s.velocities, before.velocities);
        }
    }
}
#[test]
fn report_contract_rejects_unrepresented_velocity_impulses() {
    let mut s = support::atoms(&[[0., 0., 0.], [1., 0., 0.]], &[1., 1.]);
    let c = Constraints {
        bonds: vec![Constraint { i: 0, j: 1, d: 1. }],
        waters: vec![],
    };
    let g = constraint_geometry(&s, &c).unwrap();
    let huge = 2_f64.powi(53);
    s.set_vel(0, Vec3::new(huge, 0., 0.));
    s.set_vel(1, Vec3::new(huge + 2., 0., 0.));
    let before = s.clone();
    // Both final velocities round to 2^53, so changes (0,-2) cannot equal (+1,-1).
    let result = c.rattle_report(&mut s, &g, ProjectionStage::VvFinal, 1e-13, 100);
    assert!(matches!(
        result,
        Err(modena_core::error::ModenaError::Physics(_))
    ));
    assert_eq!(s.positions, before.positions);
    assert_eq!(s.velocities, before.velocities);
}
#[test]
fn report_contract_rejects_unrepresented_position_impulses() {
    let mut s = support::atoms(&[[0., 0., 0.], [5., 0., 0.]], &[1., 3.]);
    let c = Constraints {
        bonds: vec![Constraint { i: 0, j: 1, d: 5. }],
        waters: vec![],
    };
    let g = constraint_geometry(&s, &c).unwrap();
    let huge = 2_f64.powi(53);
    s.set_pos(0, Vec3::new(huge - 2., 0., 0.));
    s.set_pos(1, Vec3::new(huge + 4., 0., 0.));
    let before = s.clone();
    // The unit change in represented length exactly satisfies SHAKE. Over a
    // half-unit interval the realized atom impulses (+2,0) differ from ±1.375.
    let result = c.shake_report(&mut s, &g, 0.5, ProjectionStage::VvPosition, 1e-13, 100);
    assert!(matches!(
        result,
        Err(modena_core::error::ModenaError::Physics(_))
    ));
    assert_eq!(s.positions, before.positions);
    assert_eq!(s.velocities, before.velocities);
}

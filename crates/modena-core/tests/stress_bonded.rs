#[path = "support/stress.rs"]
mod support;

use modena_core::bonded::{Angle, Bond, BondedForceModel, BondedTopology, Dihedral, DihedralTerm};
use modena_core::composite::{CompositeForceModel, LongRangeForceModel};
use modena_core::ewald::EwaldForceModel;
use modena_core::force_model::{ForceModel, PairForceModel, Virial};
use modena_core::neighbor::NeighborList;
use modena_core::spme::SpmeForceModel;
use modena_core::stress::{
    MissingTensor, MissingTensorReason, ModelEvaluation, SymTensor6, TensorAvailability,
};
use modena_core::system::System;
use modena_core::vec3::Vec3;
use support::{affine, atoms, complete, near, Harmonic};

fn evaluate(model: &dyn ForceModel, source: &System) -> (ModelEvaluation, Vec<f64>) {
    let mut s = source.clone();
    s.forces.fill(0.0);
    let e = model.evaluate(&mut s, &mut NeighborList::new(8.0)).unwrap();
    (e, s.forces)
}

#[test]
fn right_angle_has_zero_trace_and_nonzero_shear() {
    let mut s = atoms(&[[1., 0., 0.], [0., 0., 0.], [0., 1., 0.]], &[1., 1., 1.]);
    let mut top = BondedTopology::new();
    top.angles.push(Angle {
        i: 0,
        j: 1,
        k: 2,
        theta0: std::f64::consts::PI / 3.,
        ktheta: 3.,
    });
    let m = BondedForceModel::new(top);
    let e = ForceModel::evaluate(&m, &mut s, &mut NeighborList::new(8.)).unwrap();
    let expected = [0., 0., 0., std::f64::consts::PI / 2., 0., 0.];
    for (a, b) in complete(&e.tensor).0.into_iter().zip(expected) {
        near(a, b, 2e-12);
    }
    near(e.scalar_virial, 0., 2e-12);
}

#[test]
fn harmonic_bond_has_all_six_exact_components() {
    let s = atoms(&[[1., 2., 3.], [0., 0., 0.]], &[1., 1.]);
    let mut top = BondedTopology::new();
    top.bonds.push(Bond {
        i: 0,
        j: 1,
        r0: 0.,
        k: 1.,
    });
    let (e, _) = evaluate(&BondedForceModel::new(top), &s);
    near(e.energy, 7., 2e-12);
    near(e.scalar_virial, -14., 2e-12);
    for (a, b) in complete(&e.tensor)
        .0
        .into_iter()
        .zip([-1., -4., -9., -2., -3., -6.])
    {
        near(a, b, 2e-12);
    }
}

#[test]
fn rotated_right_angle_transforms_its_shear() {
    let e = [1. / 3., 2. / 3., 2. / 3.];
    let t = [2. / 5_f64.sqrt(), -1. / 5_f64.sqrt(), 0.];
    let s = atoms(&[e, [0.; 3], t], &[1.; 3]);
    let mut top = BondedTopology::new();
    top.angles.push(Angle {
        i: 0,
        j: 1,
        k: 2,
        theta0: std::f64::consts::PI / 3.,
        ktheta: 3.,
    });
    let (result, _) = evaluate(&BondedForceModel::new(top), &s);
    let p = std::f64::consts::PI;
    let expected = [
        p * e[0] * t[0],
        p * e[1] * t[1],
        p * e[2] * t[2],
        p / 2. * (e[0] * t[1] + e[1] * t[0]),
        p / 2. * (e[0] * t[2] + e[2] * t[0]),
        p / 2. * (e[1] * t[2] + e[2] * t[1]),
    ];
    for (a, b) in complete(&result.tensor).0.into_iter().zip(expected) {
        near(a, b, 2e-12);
    }
    near(result.scalar_virial, 0., 2e-12);
}

fn fixture() -> System {
    atoms(
        &[
            [1.2, 0.2, -0.1],
            [0.1, 0.3, 0.2],
            [-0.2, 1.5, 0.4],
            [0.5, 1.8, 1.6],
        ],
        &[1.; 4],
    )
}

fn extended_chain() -> System {
    // Each connected bond is shorter than half the box, while j-to-l is not.
    // This catches independently wrapping the last atom around the anchor.
    atoms(
        &[[-1., 1., 1.], [0., 0., 0.], [8., 1., 0.5], [15., 3.2, 1.7]],
        &[1.; 4],
    )
}

fn topologies() -> Vec<(&'static str, BondedTopology)> {
    let mut bond = BondedTopology::new();
    bond.bonds.push(Bond {
        i: 0,
        j: 1,
        r0: 0.9,
        k: 1.7,
    });
    let mut angle = BondedTopology::new();
    angle.angles.push(Angle {
        i: 0,
        j: 1,
        k: 2,
        theta0: 1.1,
        ktheta: 2.3,
    });
    let mut proper = BondedTopology::new();
    proper.dihedrals.push(Dihedral {
        i: 0,
        j: 1,
        k: 2,
        l: 3,
        terms: vec![
            DihedralTerm {
                n: 1,
                k: 0.7,
                delta: 0.2,
            },
            DihedralTerm {
                n: 3,
                k: 0.2,
                delta: -0.3,
            },
        ],
    });
    let mut improper = BondedTopology::new();
    improper.dihedrals.push(Dihedral {
        i: 0,
        j: 2,
        k: 1,
        l: 3,
        terms: vec![DihedralTerm {
            n: 2,
            k: 0.9,
            delta: std::f64::consts::PI,
        }],
    });
    vec![
        ("bond", bond),
        ("angle", angle),
        ("proper", proper),
        ("improper", improper),
    ]
}

#[test]
fn bonded_six_strain_derivatives_pass_adjacent_frozen_epsilons() {
    const EPS: [f64; 4] = [1e-4, 3e-5, 1e-5, 3e-6];
    const TOL: f64 = 2e-7;
    for (geometry, source) in [("compact", fixture()), ("extended", extended_chain())] {
        for (name, top) in topologies() {
            let model = BondedForceModel::new(top);
            let (base, _) = evaluate(&model, &source);
            let tensor = complete(&base.tensor);
            near(tensor.trace(), base.scalar_virial, 5e-12);
            let mut rows = [[0.; 6]; 4];
            for (row, eps) in EPS.into_iter().enumerate() {
                for (component, value) in rows[row].iter_mut().enumerate() {
                    let plus = evaluate(&model, &affine(&source, component, eps)).0.energy;
                    let minus = evaluate(&model, &affine(&source, component, -eps)).0.energy;
                    *value = -(plus - minus) / (2. * eps);
                }
                println!(
                    "six_strain fixture={geometry}/{name} eps={eps:.8e} analytic={:?} fd={:?}",
                    tensor.0, rows[row]
                );
            }
            for component in 0..6 {
                let expected = tensor.0[component];
                let adjacent = (0..3).any(|r| {
                    (rows[r][component] - expected).abs() <= TOL * (1. + expected.abs())
                        && (rows[r + 1][component] - expected).abs() <= TOL * (1. + expected.abs())
                        && (rows[r + 1][component] - rows[r][component]).abs()
                            <= 2. * TOL * (1. + expected.abs())
                });
                assert!(adjacent, "{geometry}/{name} component {component}: no adjacent eps pass: {rows:?}, expected {expected}");
            }
        }
    }
}

#[test]
fn connected_images_and_translation_preserve_each_bonded_tensor() {
    for source in [fixture(), extended_chain()] {
        let mut moved = source.clone();
        let shifts = [[2., -1., 3.], [-2., 1., 0.], [1., 2., -2.], [-1., -2., 1.]];
        for (i, shift) in shifts.into_iter().enumerate() {
            moved.set_pos(
                i,
                source.pos(i) + Vec3::new(0.31, -0.27, 0.19) + Vec3(shift) * 20.,
            );
        }
        for (name, top) in topologies() {
            let model = BondedForceModel::new(top);
            let (base, base_f) = evaluate(&model, &source);
            let (translated, translated_f) = evaluate(&model, &moved);
            near(base.energy, translated.energy, 2e-12);
            near(base.scalar_virial, translated.scalar_virial, 5e-12);
            for (a, b) in complete(&base.tensor)
                .0
                .into_iter()
                .zip(complete(&translated.tensor).0)
            {
                near(a, b, 2e-12);
            }
            for (a, b) in base_f.into_iter().zip(translated_f) {
                near(a, b, 2e-12);
            }
            println!("image_and_translation fixture={name} passed");
        }
    }
}

#[test]
fn extended_dihedral_rejects_independently_wrapped_anchor_negative_control() {
    let source = extended_chain();
    let top = topologies().remove(2).1;
    let (evaluation, forces) = evaluate(&BondedForceModel::new(top), &source);
    let anchor = source.pos(1);
    let wrong = (0..source.n()).fold(SymTensor6::ZERO, |sum, i| {
        let f = Vec3::new(forces[3 * i], forces[3 * i + 1], forces[3 * i + 2]);
        sum.plus(SymTensor6::dyad(
            source.box_size.min_image(source.pos(i) - anchor),
            f,
        ))
    });
    let correct = complete(&evaluation.tensor);
    assert!(correct
        .0
        .into_iter()
        .zip(wrong.0)
        .any(|(a, b)| (a - b).abs() > 1e-3));
}

#[test]
fn evaluation_adds_only_term_forces_and_compute_preserves_legacy_accumulators() {
    let source = fixture();
    let mut all = BondedTopology::new();
    let mut expected = ModelEvaluation::zero();
    let mut forces = vec![0.; source.forces.len()];
    for (_, top) in topologies() {
        let (e, f) = evaluate(&BondedForceModel::new(top.clone()), &source);
        expected = expected.combine(e);
        for (a, b) in forces.iter_mut().zip(f) {
            *a += b;
        }
        all.bonds.extend(top.bonds);
        all.angles.extend(top.angles);
        all.dihedrals.extend(top.dihedrals);
    }
    let model = BondedForceModel::new(all);
    let mut s = source.clone();
    for (i, f) in s.forces.iter_mut().enumerate() {
        *f = 0.125 * (i + 1) as f64;
    }
    let initial = s.forces.clone();
    let result = model.evaluate(&mut s, &mut NeighborList::new(8.)).unwrap();
    near(result.energy, expected.energy, 2e-12);
    near(result.scalar_virial, expected.scalar_virial, 5e-12);
    for (a, b) in complete(&result.tensor)
        .0
        .into_iter()
        .zip(complete(&expected.tensor).0)
    {
        near(a, b, 2e-12);
    }
    for ((a, b), initial) in s.forces.iter().zip(&forces).zip(&initial) {
        near(*a, b + initial, 2e-12);
    }
    let mut legacy = source.clone();
    legacy.forces = initial;
    let mut w = Virial { w: 7. };
    let energy = model
        .compute(&mut legacy, &mut w, &mut NeighborList::new(8.))
        .unwrap();
    near(energy, result.energy, 2e-12);
    near(w.w, 7. + result.scalar_virial, 5e-12);
    assert_eq!(legacy.forces, s.forces);
}

#[test]
fn combine_keeps_missing_provenance_and_known_empty_identity() {
    let e = ModelEvaluation {
        energy: 7.,
        scalar_virial: -14.,
        tensor: TensorAvailability::Complete(SymTensor6([-1., -4., -9., -2., -3., -6.])),
    };
    assert_eq!(e.clone().combine(ModelEvaluation::zero()), e);
    let missing = MissingTensor {
        component: "legacy".into(),
        reason: MissingTensorReason::LegacyForceModel,
    };
    let legacy = ModelEvaluation {
        energy: 2.,
        scalar_virial: 3.,
        tensor: TensorAvailability::Missing(vec![missing.clone()]),
    };
    let sum = e.combine(legacy);
    assert_eq!(sum.energy, 9.);
    assert_eq!(sum.scalar_virial, -11.);
    assert_eq!(
        sum.tensor,
        TensorAvailability::Missing(vec![missing.clone()])
    );
    let other = MissingTensor {
        component: "other".into(),
        reason: MissingTensorReason::UnsupportedBackend,
    };
    let both = sum.combine(ModelEvaluation {
        energy: 0.,
        scalar_virial: 0.,
        tensor: TensorAvailability::Missing(vec![other.clone()]),
    });
    assert_eq!(
        both.tensor,
        TensorAvailability::Missing(vec![missing, other])
    );
}

#[test]
fn known_empty_bonded_and_composite_are_complete_zero() {
    let source = atoms(&[], &[]);
    let bonded = BondedForceModel::new(BondedTopology::new());
    assert_eq!(evaluate(&bonded, &source).0, ModelEvaluation::zero());
    assert_eq!(evaluate(&bonded, &fixture()).0, ModelEvaluation::zero());
    let composite = CompositeForceModel::new(PairForceModel::new(Harmonic)).with_bonded(bonded);
    assert_eq!(evaluate(&composite, &source).0, ModelEvaluation::zero());
}

#[test]
fn concrete_composite_combines_pair_bonded_and_long_range_once() {
    let mut source = fixture();
    source.charges = vec![1., -1., 0.2, -0.2];
    let top = topologies().remove(1).1;
    let pair = PairForceModel::new(Harmonic);
    let bonded = BondedForceModel::new(top.clone());
    let (p, pf) = evaluate(&pair, &source);
    let (b, bf) = evaluate(&bonded, &source);
    let expected = p.clone().combine(b.clone());
    let composite = CompositeForceModel::new(pair).with_bonded(bonded);
    let (actual, forces) = evaluate(&composite, &source);
    assert_eq!(actual, expected);
    for ((actual, p), b) in forces.into_iter().zip(&pf).zip(&bf) {
        near(actual, p + b, 2e-12);
    }
    for long_range in [
        LongRangeForceModel::Ewald(EwaldForceModel::new(0.8, 3.5, 3.).unwrap()),
        LongRangeForceModel::Spme(SpmeForceModel::new(0.8, 3.5, 3.).unwrap()),
    ] {
        let (underlying, uf) = match &long_range {
            LongRangeForceModel::Ewald(m) => evaluate(m, &source),
            LongRangeForceModel::Spme(m) => evaluate(m, &source),
        };
        let (forwarded, forwarded_f) = evaluate(&long_range, &source);
        assert_eq!(forwarded, underlying);
        assert_eq!(forwarded_f, uf);
        match &long_range {
            LongRangeForceModel::Ewald(_) => {
                let tensor = complete(&forwarded.tensor);
                near(tensor.trace(), forwarded.scalar_virial, 5e-12);
            }
            LongRangeForceModel::Spme(_) => {
                assert!(matches!(forwarded.tensor, TensorAvailability::Missing(_)))
            }
        }
        let expected = p.clone().combine(b.clone()).combine(underlying);
        let combined = CompositeForceModel {
            pair: PairForceModel::new(Harmonic),
            bonded: Some(BondedForceModel::new(top.clone())),
            long_range: Some(long_range),
        };
        let (actual, forces) = evaluate(&combined, &source);
        assert_eq!(actual, expected);
        for (((actual, p), b), u) in forces.into_iter().zip(&pf).zip(&bf).zip(uf) {
            near(actual, p + b + u, 2e-12);
        }
    }
}

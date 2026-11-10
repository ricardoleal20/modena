#[path = "support/stress.rs"]
mod support;

use modena_core::force_model::{ForceModel, PairForceModel, Virial};
use modena_core::neighbor::NeighborList;
use modena_core::stress::{
    MissingTensor, MissingTensorReason, ModelEvaluation, SymTensor6, TensorAvailability,
};
use modena_core::{LennardJones, PairScale, System, Vec3};
use support::{affine, atoms, complete, near, Harmonic};

#[test]
fn pair_tensor_is_additive_and_not_isotropic() {
    let mut s = atoms(&[[0., 0., 0.], [1., 2., 3.]], &[1., 3.]);
    s.forces.fill(7.);
    let model = PairForceModel::new(Harmonic);
    let mut nl = NeighborList::new(8.);
    let a = ForceModel::evaluate(&model, &mut s, &mut nl).unwrap();
    assert_eq!(a.energy, 7.);
    assert_eq!(a.scalar_virial, -14.);
    assert_eq!(complete(&a.tensor).0, [-1., -4., -9., -2., -3., -6.]);
    assert_eq!(s.forces, vec![8., 9., 10., 6., 5., 4.]);
    let b = ForceModel::evaluate(&model, &mut s, &mut nl).unwrap();
    assert_eq!(b, a);
    assert_eq!(s.forces, vec![9., 11., 13., 5., 3., 1.]);
}

struct Legacy;
impl ForceModel for Legacy {
    fn compute(
        &self,
        s: &mut System,
        w: &mut Virial,
        _: &mut NeighborList,
    ) -> modena_core::Result<f64> {
        for f in &mut s.forces {
            *f += 1.0;
        }
        w.w += 3.0;
        Ok(2.0)
    }
}

#[test]
fn legacy_evaluation_is_missing_even_for_empty_systems() {
    for mut s in [atoms(&[[0., 0., 0.]], &[1.]), atoms(&[], &[])] {
        s.forces.fill(7.0);
        let a = (&Legacy as &dyn ForceModel)
            .evaluate(&mut s, &mut NeighborList::new(8.))
            .unwrap();
        assert_eq!(a.energy, 2.0);
        assert_eq!(a.scalar_virial, 3.0);
        assert_eq!(s.forces, vec![8.0; 3 * s.n()]);
        assert_eq!(
            a.tensor,
            TensorAvailability::Missing(vec![MissingTensor {
                component: std::any::type_name::<Legacy>().into(),
                reason: MissingTensorReason::LegacyForceModel,
            }])
        );
    }
}

#[test]
fn known_empty_pair_is_complete_zero() {
    let model = PairForceModel::new(Harmonic);
    for mut s in [atoms(&[], &[]), atoms(&[[1., 2., 3.]], &[1.])] {
        s.forces.fill(7.0);
        assert_eq!(
            model.evaluate(&mut s, &mut NeighborList::new(8.)).unwrap(),
            ModelEvaluation::zero()
        );
        assert_eq!(s.forces, vec![7.0; 3 * s.n()]);
    }
}

#[test]
fn symmetric_dyad_and_combination_preserve_all_components_and_missing_reasons() {
    let t = SymTensor6::dyad(Vec3::new(1., 2., 3.), Vec3::new(4., 5., 6.));
    assert_eq!(t.0, [4., 10., 18., 6.5, 9., 13.5]);
    assert_eq!(t.trace(), 32.);
    assert_eq!(t.scaled(2.).plus(t).0, [12., 30., 54., 19.5, 27., 40.5]);
    assert!(t.finite());
    assert!(!SymTensor6([f64::INFINITY; 6]).finite());
    let a = ModelEvaluation {
        energy: 2.,
        scalar_virial: 32.,
        tensor: TensorAvailability::Complete(t),
    };
    assert_eq!(a.clone().combine(ModelEvaluation::zero()), a);
    assert_eq!(
        a.clone().combine(a.clone()),
        ModelEvaluation {
            energy: 4.,
            scalar_virial: 64.,
            tensor: TensorAvailability::Complete(t.scaled(2.)),
        }
    );
    let missing = |component: &str| ModelEvaluation {
        energy: 1.,
        scalar_virial: 3.,
        tensor: TensorAvailability::Missing(vec![MissingTensor {
            component: component.into(),
            reason: MissingTensorReason::UnsupportedBackend,
        }]),
    };
    for combined in [
        a.clone().combine(missing("first")),
        missing("first").combine(a),
    ] {
        assert_eq!(combined.energy, 3.);
        assert_eq!(combined.scalar_virial, 35.);
        assert_eq!(combined.tensor, missing("first").tensor);
    }
    assert_eq!(
        missing("first").combine(missing("second")).tensor,
        TensorAvailability::Missing(vec![
            MissingTensor {
                component: "first".into(),
                reason: MissingTensorReason::UnsupportedBackend
            },
            MissingTensor {
                component: "second".into(),
                reason: MissingTensorReason::UnsupportedBackend
            },
        ])
    );
}

#[test]
fn tensor_records_round_trip_and_reject_nonfinite_or_unknown_fields() {
    let a = ModelEvaluation {
        energy: 7.,
        scalar_virial: -14.,
        tensor: TensorAvailability::Complete(SymTensor6([-1., -4., -9., -2., -3., -6.])),
    };
    let json = serde_json::to_string(&a).unwrap();
    assert!(json.contains("\"kind\":\"complete\""));
    assert_eq!(serde_json::from_str::<ModelEvaluation>(&json).unwrap(), a);
    assert!(serde_json::from_str::<ModelEvaluation>(
        &json.replace("\"energy\":7.0", "\"extra\":0,\"energy\":7.0")
    )
    .is_err());
    for field in ["energy", "scalar_virial", "tensor"] {
        let energy = if field == "energy" { "nan" } else { "7.0" };
        let w = if field == "scalar_virial" {
            "inf"
        } else {
            "-14.0"
        };
        let t = if field == "tensor" {
            "[0.0, 0.0, 0.0, 0.0, 0.0, -inf]"
        } else {
            "[0.0, 0.0, 0.0, 0.0, 0.0, 0.0]"
        };
        let text = format!(
            "energy = {energy}\nscalar_virial = {w}\n[tensor]\nkind = 'complete'\nvalue = {t}\n"
        );
        assert!(
            toml::from_str::<ModelEvaluation>(&text).is_err(),
            "accepted nonfinite {field}"
        );
    }
    let missing = TensorAvailability::Missing(vec![MissingTensor {
        component: "legacy".into(),
        reason: MissingTensorReason::LegacyForceModel,
    }]);
    let json = serde_json::to_string(&missing).unwrap();
    assert!(json.contains("legacy_force_model"));
    assert_eq!(
        serde_json::from_str::<TensorAvailability>(&json).unwrap(),
        missing
    );
    assert!(serde_json::from_str::<MissingTensor>(
        r#"{"component":"legacy","reason":{"kind":"legacy_force_model"},"extra":0}"#
    )
    .is_err());
}

fn six_strain_gate(
    model: &PairForceModel,
    source: &System,
    cutoff: f64,
    scales: &[PairScale],
    label: &str,
    threads: usize,
    parallel: bool,
) {
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .unwrap();
    pool.install(|| {
        let mut s = source.clone();
        s.forces.fill(7.0);
        let mut nl = NeighborList::with_pair_scales(cutoff, 0.3, scales).unwrap();
        let cold = model.evaluate(&mut s, &mut nl).unwrap();
        if parallel { assert!(nl.n_pairs() >= 100_000); }
        let once = s.forces.clone();
        let warm = model.evaluate(&mut s, &mut nl).unwrap();
        assert_eq!(warm, cold, "reused scratch changed evaluation");
        for (f, once) in s.forces.iter().zip(once) { near(*f, 2.*once-7., 5e-12); }
        let t = complete(&cold.tensor);
        assert!(t.finite() && cold.energy.is_finite() && cold.scalar_virial.is_finite());
        let scale = t.0.iter().map(|x| x.abs()).sum::<f64>();
        assert!((t.trace()-cold.scalar_virial).abs() <= 5e-12*(1.+scale));
        let mut deterministic = source.clone();
        deterministic.forces.fill(7.);
        let mut w = Virial { w: 11. };
        let e = model.compute_deterministic(&mut deterministic, &mut w, &mut nl).unwrap();
        near(cold.energy, e, 5e-12);
        near(cold.scalar_virial, w.w-11., 5e-12);
        for (a,b) in deterministic.forces.iter().zip(&s.forces) { near(*b, 2.*a-7., 5e-12); }
        for component in 0..6 {
            let mut derivatives = Vec::new();
            for eps in [1e-4,3e-5,1e-5,3e-6] {
                let energy = |eta| {
                    let mut deformed = affine(source, component, eta);
                    let mut list = NeighborList::with_pair_scales(cutoff, 0.3, scales).unwrap();
                    model.compute(&mut deformed, &mut Virial::default(), &mut list).unwrap()
                };
                let fd = -(energy(eps)-energy(-eps))/(2.*eps);
                assert!(fd.is_finite());
                eprintln!("six_strain label={label} threads={threads} component={component} eps={eps:e} fd={fd:.17e} tensor={:.17e}", t.0[component]);
                derivatives.push(fd);
            }
            let expected = t.0[component];
            assert!(derivatives.windows(2).any(|pair| {
                pair.iter().all(|x| (x-expected).abs() <= 2e-7*(1.+expected.abs()))
                    && (pair[0]-pair[1]).abs() <= 4e-7*(1.+expected.abs())
            }), "no adjacent eps passed: {label} component={component} {derivatives:?} expected={expected}");
        }
    });
}

#[test]
fn pair_six_strains_periodic_scaled_and_parallel_reused_scratch() {
    let periodic = atoms(
        &[[0.2, 0.3, 0.5], [19.1, 19.4, 19.6], [1.5, 1.9, 1.7]],
        &[1., 2., 3.],
    );
    let scales = [PairScale {
        i: 0,
        j: 1,
        lj: 0.25,
        coulomb: 1.,
    }];
    let points: Vec<_> = (0..450)
        .map(|i| {
            let x = i as f64 / 450.;
            [1. + x, 2. + x * x, 3. + x * x * x]
        })
        .collect();
    let large = atoms(&points, &vec![1.; 450]);
    for threads in [1, 4] {
        six_strain_gate(
            &PairForceModel::new(Harmonic),
            &periodic,
            8.,
            &scales,
            "periodic_scaled_harmonic",
            threads,
            false,
        );
        six_strain_gate(
            &PairForceModel::new(LennardJones::new(1., 1., 4.).unwrap()),
            &periodic,
            4.,
            &scales,
            "periodic_scaled_lj",
            threads,
            false,
        );
        six_strain_gate(
            &PairForceModel::new(Harmonic),
            &large,
            8.,
            &scales,
            "parallel_harmonic",
            threads,
            true,
        );
    }
}

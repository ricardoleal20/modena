use modena_core::decomp::DomainForceModel;
use modena_core::potential::{InterpKind, Tabulated};
use modena_core::{
    CompositeForceModel, EwaldForceModel, ForceModel, LennardJones, NeighborList, PairForceModel,
    PairPotential, SimBox, SpmeForceModel, System, Vec3, Virial,
};

fn charged_system() -> System {
    let mut s = System::new(SimBox::cubic(9.0));
    s.resize(3);
    for (i, p) in [
        Vec3::new(1.1, 2.3, 3.4),
        Vec3::new(2.7, 3.2, 3.8),
        Vec3::new(6.2, 5.8, 7.0),
    ]
    .into_iter()
    .enumerate()
    {
        s.set_pos(i, p);
    }
    s.charges = vec![1.0, -0.7, -0.3];
    s
}

fn eval(model: &dyn ForceModel, s: &System) -> (f64, f64, Vec<f64>) {
    let mut s = s.clone();
    s.zero_forces();
    let mut v = Virial::default();
    let mut nl = NeighborList::new(4.2);
    let u = model.compute(&mut s, &mut v, &mut nl).unwrap();
    (u, v.w, s.forces)
}

fn strain(s: &System, lam: f64) -> System {
    let mut s = s.clone();
    s.box_size.scale_uniform(lam);
    for x in &mut s.positions {
        *x *= lam;
    }
    s
}

#[test]
fn tabulated_joint_gradient_including_lower_extrapolation() {
    // A cubic and its analytical derivative must be reproduced by endpoint
    // Hermite data, including shifted-force cutoff. Below the table, the
    // endpoint tangent must be extended consistently.
    let t = Tabulated::new(
        vec![1.0, 2.0, 3.0],
        vec![1.0, 8.0, 27.0],
        vec![-3.0, -12.0, -27.0],
        InterpKind::Linear,
    )
    .unwrap();
    for r in [0.7_f64, 1.3, 1.8, 2.4, 2.999] {
        let h = 1e-6;
        let (u, fr) = t.energy_force(r * r, 0, 0);
        let fd = -(t.energy_force((r + h).powi(2), 0, 0).0
            - t.energy_force((r - h).powi(2), 0, 0).0)
            / (2.0 * h);
        assert!(
            (fr * r - fd).abs() < 1e-7,
            "r={r}: force={} derivative={fd}",
            fr * r
        );
        if r >= 1.0 {
            assert!((u - (r.powi(3) - 27.0 - (r - 3.0) * 27.0)).abs() < 1e-12);
        }
    }
}

#[test]
fn direct_ewald_virial_is_isotropic_strain_derivative() {
    let m = EwaldForceModel::new(0.8, 4.2, 6.4).unwrap();
    let s = charged_system();
    let h = 1e-6;
    let (_, w, _) = eval(&m, &s);
    let fd = -(eval(&m, &strain(&s, 1.0 + h)).0 - eval(&m, &strain(&s, 1.0 - h)).0) / (2.0 * h);
    assert!((fd - w).abs() < 1e-7, "W={w} strain derivative={fd}");
}

#[test]
fn charged_background_virial_is_three_times_energy() {
    let mut s = System::new(SimBox::cubic(9.0));
    s.resize(1);
    s.charges[0] = 1.3;
    s.set_pos(0, Vec3::new(2.0, 3.0, 4.0));
    let m = EwaldForceModel::new(0.8, 4.2, 0.1).unwrap();
    let (_, w, _) = eval(&m, &s);
    let ubg = -std::f64::consts::PI * 1.3_f64.powi(2) / (2.0 * 729.0 * 0.8_f64.powi(2));
    assert!(
        (w - 3.0 * ubg).abs() < 1e-14,
        "W={w}, expected {}",
        3.0 * ubg
    );
}

#[test]
fn spme_full_net_charge_virial_is_strain_derivative() {
    let mut s = charged_system();
    s.charges[2] = 0.4;
    let m = SpmeForceModel::new(0.8, 4.2, 6.4).unwrap();
    let h = 1e-6;
    let (_, w, _) = eval(&m, &s);
    let fd = -(eval(&m, &strain(&s, 1.0 + h)).0 - eval(&m, &strain(&s, 1.0 - h)).0) / (2.0 * h);
    assert!((w - fd).abs() < 1e-7, "W={w} strain derivative={fd}");
}

#[test]
fn excluded_ewald_energy_only_never_writes_a_force_slice() {
    let s = charged_system();
    let m = EwaldForceModel::new(0.8, 4.2, 6.4)
        .unwrap()
        .with_excluded_pairs(vec![(0, 1)]);
    let mut nl = NeighborList::with_skin_exclusions(4.2, 0.0, &[(0, 1)]);
    let u = m.energy(&s, &mut nl).unwrap();
    let mut ss = s.clone();
    let mut v = Virial::default();
    assert!((u - m.compute(&mut ss, &mut v, &mut nl).unwrap()).abs() < 1e-12);
}

#[test]
fn triclinic_neighbor_search_finds_nonadjacent_cartesian_bins() {
    let mut s = System::new(SimBox::triclinic([
        Vec3::new(10.0, 0.0, 0.0),
        Vec3::new(49.0, 10.0, 0.0),
        Vec3::new(0.0, 0.0, 40.0),
    ]));
    s.resize(2);
    s.set_pos(0, Vec3::zero());
    s.set_pos(1, Vec3::new(19.7, 0.1, 0.0));
    let mut nl = NeighborList::new(4.8);
    nl.build(&s);
    assert_eq!(
        nl.pairs(),
        &[(0, 1)],
        "lattice image (-0.3,0.1,0) is inside rc"
    );
    s.box_size = SimBox::triclinic([
        Vec3::new(10.0, 0.0, 0.0),
        Vec3::new(4.9, 10.0, 0.0),
        Vec3::new(0.0, 0.0, 40.0),
    ]);
    s.set_pos(1, Vec3::new(6.37, 3.0, 0.0));
    nl.force_rebuild(&s);
    assert_eq!(
        nl.pairs(),
        &[(0, 1)],
        "audited image (-3.63,3,0) is inside rc"
    );
}

#[test]
fn ewald_reciprocal_cache_keys_full_cell_basis() {
    let m = EwaldForceModel::new(0.8, 4.2, 6.4).unwrap();
    let mut s = charged_system();
    s.box_size = SimBox::triclinic([
        Vec3::new(9.0, 0.0, 0.0),
        Vec3::new(3.0, 9.0, 0.0),
        Vec3::new(0.0, 0.0, 9.0),
    ]);
    eval(&m, &s);
    s.box_size = SimBox::triclinic([
        Vec3::new(9.0, 0.0, 0.0),
        Vec3::new(-3.0, 9.0, 0.0),
        Vec3::new(0.0, 0.0, 9.0),
    ]);
    let fresh = EwaldForceModel::new(0.8, 4.2, 6.4).unwrap();
    assert!(
        (eval(&m, &s).0 - eval(&fresh, &s).0).abs() < 1e-12,
        "same bounding box does not mean same reciprocal basis"
    );
}

#[test]
fn all_force_models_accumulate_existing_force_and_virial() {
    let mut s = System::new(SimBox::cubic(10.0));
    s.resize(2);
    s.set_pos(1, Vec3::new(1.3, 0.2, 0.1));
    let models: Vec<Box<dyn ForceModel>> = vec![
        Box::new(DomainForceModel::single_rank(
            LennardJones::new(1.0, 1.0, 4.2).unwrap(),
            s.box_size,
            4.2,
        )),
        Box::new(CompositeForceModel::new(PairForceModel::new(
            LennardJones::new(1.0, 1.0, 4.2).unwrap(),
        ))),
    ];
    for m in models {
        let (_, w, f) = eval(&*m, &s);
        let mut ss = s.clone();
        ss.forces.fill(7.0);
        let mut v = Virial { w: 13.0 };
        m.compute(&mut ss, &mut v, &mut NeighborList::new(4.2))
            .unwrap();
        assert!((v.w - w - 13.0).abs() < 1e-12, "existing virial erased");
        for (a, b) in ss.forces.iter().zip(f) {
            assert!((a - b - 7.0).abs() < 1e-12, "existing force erased");
        }
    }
}

#[test]
fn concrete_cutoffs_are_real_model_requirements() {
    let pair = PairForceModel::new(LennardJones::new(1.0, 1.0, 3.7).unwrap());
    assert_eq!(pair.neighbor_cutoff(), Some(3.7));
    let c = CompositeForceModel::new(pair).with_ewald(EwaldForceModel::new(0.8, 4.2, 6.4).unwrap());
    assert_eq!(c.neighbor_cutoff(), Some(4.2));
}

#[test]
fn electrostatic_constructors_reject_nonfinite_parameters() {
    for x in [f64::INFINITY, f64::NAN, -1.0, 0.0] {
        for p in [(x, 4.2, 6.4), (0.8, x, 6.4), (0.8, 4.2, x)] {
            assert!(EwaldForceModel::new(p.0, p.1, p.2).is_err());
            assert!(SpmeForceModel::new(p.0, p.1, p.2).is_err());
        }
    }
}

#[test]
fn empty_system_invalidates_neighbor_cache_before_repopulation() {
    let mut s = System::new(SimBox::cubic(12.0));
    s.resize(2);
    s.set_pos(1, Vec3::new(1.0, 0.0, 0.0));
    let mut nl = NeighborList::with_skin(4.2, 0.3);
    nl.build(&s);
    assert_eq!(nl.pairs(), &[(0, 1)]);
    s.resize(0);
    nl.build(&s);
    assert!(nl.pairs().is_empty());
    s.resize(2);
    s.set_pos(1, Vec3::new(1.0, 0.0, 0.0));
    nl.build(&s);
    assert_eq!(nl.pairs(), &[(0, 1)]);
}

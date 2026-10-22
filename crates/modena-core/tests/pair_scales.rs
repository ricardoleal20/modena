use modena_core::{
    EwaldForceModel, ForceModel, LennardJones, NeighborList, PairForceModel, PairPotential,
    PairScale, SimBox, SpmeForceModel, System, Vec3, Virial,
};

fn system(r: f64) -> System {
    let mut s = System::new(SimBox::cubic(12.0));
    s.resize(2);
    s.set_pos(0, Vec3::new(1.0, 2.0, 3.0));
    s.set_pos(1, Vec3::new(1.0 + r, 2.0, 3.0));
    s.charges = vec![1.0, -0.7];
    s
}
fn compute(m: &dyn ForceModel, s: &System, scales: &[PairScale]) -> (f64, f64, Vec<f64>) {
    let mut s = s.clone();
    s.zero_forces();
    let mut v = Virial::default();
    let mut n = NeighborList::with_pair_scales(4.2, 0.3, scales).unwrap();
    let u = m.compute(&mut s, &mut v, &mut n).unwrap();
    (u, v.w, s.forces)
}

#[test]
fn independent_scales_keep_a_pair_when_only_one_interaction_is_zero() {
    for (lj, coulomb) in [(0.0, 1.0), (1.0, 0.0), (0.3, 0.7), (0.0, 0.0)] {
        let scale = PairScale {
            i: 1,
            j: 0,
            lj,
            coulomb,
        };
        let s = system(1.3);
        let mut n = NeighborList::with_pair_scales(4.2, 0.3, &[scale]).unwrap();
        n.build(&s);
        assert_eq!(n.pair_scale(0, 1), (lj, coulomb));
        assert_eq!(n.pair_scale(1, 0), (lj, coulomb));
        assert_eq!(n.n_pairs(), usize::from(lj != 0.0 || coulomb != 0.0));
        let pot = LennardJones::new(1.0, 1.0, 4.2).unwrap();
        let (bare_u, bare_fr) = pot.energy_force(1.3_f64.powi(2), 0, 0);
        let (u, w, f) = compute(&PairForceModel::new(pot), &s, &[scale]);
        assert!((u - lj * bare_u).abs() < 1e-13);
        assert!((f[0] + lj * bare_fr * 1.3).abs() < 1e-13);
        assert!((w - lj * bare_fr * 1.3_f64.powi(2)).abs() < 1e-13);
    }
}

#[test]
fn constructors_reject_invalid_and_conflicting_pair_policies() {
    let p = PairScale {
        i: 0,
        j: 1,
        lj: 0.3,
        coulomb: 0.7,
    };
    for bad in [
        PairScale { i: 0, j: 0, ..p },
        PairScale { lj: f64::NAN, ..p },
        PairScale { coulomb: -0.1, ..p },
        PairScale {
            lj: f64::INFINITY,
            ..p
        },
    ] {
        assert!(NeighborList::with_pair_scales(4.2, 0.3, &[bad]).is_err());
        assert!(EwaldForceModel::new(0.8, 4.2, 6.4)
            .unwrap()
            .with_pair_scales(vec![bad])
            .is_err());
        assert!(SpmeForceModel::new(0.8, 4.2, 6.4)
            .unwrap()
            .with_pair_scales(vec![bad])
            .is_err());
    }
    assert!(NeighborList::with_pair_scales(
        4.2,
        0.3,
        &[
            p,
            PairScale {
                i: 1,
                j: 0,
                lj: 0.5,
                ..p
            }
        ]
    )
    .is_err());
    assert!(NeighborList::with_pair_scales(4.2, 0.3, &[p, PairScale { i: 1, j: 0, ..p }]).is_ok());
    assert!(NeighborList::with_pair_scales(f64::INFINITY, 0.3, &[]).is_err());
}

#[test]
fn weighted_ewald_correction_removes_bare_fraction_inside_and_outside_real_cutoff() {
    for r in [1.3, 5.0] {
        for c in [0.0, 0.4, 1.0] {
            let s = system(r);
            let p = PairScale {
                i: 0,
                j: 1,
                lj: 1.0,
                coulomb: c,
            };
            let plain = EwaldForceModel::new(0.8, 4.2, 6.4).unwrap();
            let scaled = EwaldForceModel::new(0.8, 4.2, 6.4)
                .unwrap()
                .with_pair_scales(vec![p])
                .unwrap();
            let (u0, w0, f0) = compute(&plain, &s, &[]);
            let (u, w, f) = compute(&scaled, &s, &[p]);
            let bare = -0.7 / r;
            let removed = (1.0 - c) * bare;
            assert!(
                (u - (u0 - removed)).abs() < 1e-12,
                "r={r} c={c}: U={u} expected {}",
                u0 - removed
            );
            assert!((w - (w0 - removed)).abs() < 1e-12);
            assert!((f[0] - (f0[0] - (1.0 - c) * 0.7 / (r * r))).abs() < 1e-12);
            // A neighbor-only policy must be equally valid; duplicate model/list
            // policy is applied once, and the energy-only path matches compute.
            assert!((compute(&plain, &s, &[p]).0 - u).abs() < 1e-12);
            let mut n = NeighborList::with_pair_scales(4.2, 0.3, &[p]).unwrap();
            assert!((scaled.energy(&s, &mut n).unwrap() - u).abs() < 1e-12);
        }
    }
}

#[test]
fn spme_weighted_correction_and_strain_match_bare_fraction() {
    for r in [1.3, 5.0] {
        let s = system(r);
        let p = PairScale {
            i: 0,
            j: 1,
            lj: 0.0,
            coulomb: 0.4,
        };
        let plain = SpmeForceModel::new(0.8, 4.2, 6.4).unwrap();
        let scaled = SpmeForceModel::new(0.8, 4.2, 6.4)
            .unwrap()
            .with_pair_scales(vec![p])
            .unwrap();
        let (u0, w0, f0) = compute(&plain, &s, &[]);
        let (u, w, f) = compute(&scaled, &s, &[p]);
        let removed = 0.6 * (-0.7) / r;
        assert!((u - u0 + removed).abs() < 1e-12);
        assert!((w - w0 + removed).abs() < 1e-12);
        assert!((f[0] - f0[0] + 0.6 * 0.7 / (r * r)).abs() < 1e-12);
        let h = 1e-6;
        let strain = |lambda: f64| {
            let mut t = s.clone();
            t.box_size.scale_uniform(lambda);
            for x in &mut t.positions {
                *x *= lambda;
            }
            compute(&scaled, &t, &[p]).0
        };
        assert!((w + (strain(1.0 + h) - strain(1.0 - h)) / (2.0 * h)).abs() < 1e-7);
    }
}

#[test]
fn conflicting_model_and_neighbor_scales_and_bad_indices_are_errors() {
    let p = PairScale {
        i: 0,
        j: 1,
        lj: 1.0,
        coulomb: 0.4,
    };
    let m = EwaldForceModel::new(0.8, 4.2, 6.4)
        .unwrap()
        .with_pair_scales(vec![p])
        .unwrap();
    let mut s = system(1.3);
    let mut n =
        NeighborList::with_pair_scales(4.2, 0.3, &[PairScale { coulomb: 0.7, ..p }]).unwrap();
    assert!(m.compute(&mut s, &mut Virial::default(), &mut n).is_err());
    let mut n = NeighborList::with_pair_scales(4.2, 0.3, &[PairScale { j: 10, ..p }]).unwrap();
    assert!(m.compute(&mut s, &mut Virial::default(), &mut n).is_err());
}

#[test]
fn domain_backend_rejects_special_pair_policy_and_changed_box() {
    let mut s = system(1.3);
    let m = modena_core::decomp::DomainForceModel::single_rank(
        LennardJones::new(1.0, 1.0, 4.2).unwrap(),
        s.box_size,
        4.2,
    );
    let p = PairScale {
        i: 0,
        j: 1,
        lj: 0.5,
        coulomb: 1.0,
    };
    let mut n = NeighborList::with_pair_scales(4.2, 0.3, &[p]).unwrap();
    assert!(m
        .compute(&mut s, &mut Virial::default(), &mut n)
        .unwrap_err()
        .to_string()
        .contains("special pair"));
    s.box_size.scale_uniform(1.1);
    assert!(m
        .compute(&mut s, &mut Virial::default(), &mut NeighborList::new(4.2))
        .unwrap_err()
        .to_string()
        .contains("original box"));
}

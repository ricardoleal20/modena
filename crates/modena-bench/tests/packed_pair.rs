use std::collections::BTreeSet;

use modena_bench::packed_pair::{ClusterView, Layout, PackedForceModel, Permutation};
use modena_core::{ForceModel, NeighborList, PairScale, SimBox, System, Vec3, Virial};

const MODES: [Layout; 4] = [
    Layout::Interleaved,
    Layout::Soa,
    Layout::Aosoa4,
    Layout::Cluster4,
];

fn fixture(n: usize) -> System {
    let mut s = System::new(SimBox::cubic(12.));
    s.resize(n);
    for i in 0..n {
        // Reverse spatial order deliberately; never confuse a lane with an ID.
        s.set_pos(
            i,
            Vec3::new(
                0.2 + 1.15 * ((n - i - 1) % 5) as f64,
                0.3 + 1.17 * (i / 5) as f64,
                0.4,
            ),
        );
    }
    s
}

// Independent scalar half-pair oracle: no production potential, minimum-image
// helper, packed coordinates, masks or kernel arithmetic is reused here.
fn oracle(s: &System, nl: &NeighborList) -> (f64, f64, Vec<f64>) {
    let mut f = vec![0.; 3 * s.n()];
    let (mut e, mut w) = (0., 0.);
    let rc = 2.5_f64;
    let uc = 4. * (rc.powi(-12) - rc.powi(-6));
    let fc = 24. * (2. * rc.powi(-13) - rc.powi(-7));
    for &(i, j) in nl.pairs() {
        let mut d = [0.; 3];
        for (k, x) in d.iter_mut().enumerate() {
            let raw = s.positions[3 * i + k] - s.positions[3 * j + k];
            let l = s.box_size.size.0[k];
            *x = raw - l * (raw / l).round();
        }
        let r2 = d.iter().map(|x| x * x).sum::<f64>();
        if !(1e-14..6.25).contains(&r2) {
            continue;
        }
        let r = r2.sqrt();
        let scale = nl.pair_scale(i, j).0;
        e += scale * (4. * (r.powi(-12) - r.powi(-6)) - uc + (r - rc) * fc);
        let fr = scale * (24. * (2. * r.powi(-13) - r.powi(-7)) - fc) / r;
        for k in 0..3 {
            let force = d[k] * fr;
            f[3 * i + k] += force;
            f[3 * j + k] -= force;
            w += d[k] * force;
        }
    }
    (e, w, f)
}

fn near(a: f64, b: f64, abs: f64) {
    assert!((a - b).abs() <= abs + 1e-12 * b.abs(), "{a} versus {b}");
}

fn gate(model: &PackedForceModel, s: &mut System, nl: &mut NeighborList) {
    nl.build(s);
    let (e, w, f) = oracle(s, nl);
    s.forces.fill(0.125);
    let mut v = Virial { w: 0.375 };
    let actual = model.compute(s, &mut v, nl).unwrap();
    if let Some(pairs) = model.cached_cluster_pairs() {
        let expected: BTreeSet<_> = nl
            .pairs()
            .iter()
            .map(|&(i, j)| (i, j, nl.pair_scale(i, j).0.to_bits()))
            .collect();
        assert_eq!(pairs.len(), expected.len());
        assert_eq!(pairs.into_iter().collect::<BTreeSet<_>>(), expected);
    }
    near(actual, e, 1e-8);
    near(v.w, w + 0.375, 1e-8);
    for (actual, expected) in s.forces.iter().zip(f) {
        near(*actual, expected + 0.125, 1e-10);
    }
}

#[test]
fn replacement_list_same_epoch_and_pair_count_cannot_reuse_stale_masks() {
    let model = PackedForceModel::new(Layout::Cluster4);
    let mut s = fixture(5);
    let mut nl = NeighborList::with_pair_scales(
        2.5,
        0.3,
        &[PairScale {
            i: 0,
            j: 1,
            lj: 0.25,
            coulomb: 1.,
        }],
    )
    .unwrap();
    gate(&model, &mut s, &mut nl);
    let epoch = nl.rebuild_count();
    let count = nl.n_pairs();
    let reconstructions = model.stats().cluster_reconstructions;
    nl = NeighborList::with_pair_scales(
        2.5,
        0.3,
        &[PairScale {
            i: 0,
            j: 1,
            lj: 0.75,
            coulomb: 1.,
        }],
    )
    .unwrap();
    nl.build(&s);
    assert_eq!(nl.rebuild_count(), epoch);
    assert_eq!(nl.n_pairs(), count);
    gate(&model, &mut s, &mut nl);
    assert_eq!(model.stats().cluster_reconstructions, reconstructions + 1);
    nl.force_rebuild(&s);
    gate(&model, &mut s, &mut nl);
    assert_eq!(model.stats().cluster_reconstructions, reconstructions + 2);
}

#[test]
fn exact_roundtrip_and_scaled_mask_set_with_padding() {
    for n in [0, 1, 2, 3, 4, 5, 17] {
        let s = fixture(n);
        let p = Permutation::spatial(&s).unwrap();
        let logical_ids: Vec<_> = (0..n).collect();
        assert_eq!(
            p.unpack(&p.pack(&logical_ids).unwrap()).unwrap(),
            logical_ids
        );
        if n == 5 {
            assert_ne!(p.logical_ids()[0], Some(0));
        }
        let scales = if n >= 3 {
            vec![
                PairScale {
                    i: 0,
                    j: 1,
                    lj: 0.,
                    coulomb: 0.,
                },
                PairScale {
                    i: 1,
                    j: 2,
                    lj: 0.25,
                    coulomb: 1.,
                },
                PairScale {
                    i: 0,
                    j: 2,
                    lj: 0.,
                    coulomb: 1.,
                },
            ]
        } else {
            vec![]
        };
        let mut nl = NeighborList::with_pair_scales(2.5, 0.3, &scales).unwrap();
        nl.build(&s);
        let view = ClusterView::new(&s, &nl).unwrap();
        let decoded = view.logical_pairs();
        let expected: BTreeSet<_> = nl
            .pairs()
            .iter()
            .map(|&(i, j)| (i, j, nl.pair_scale(i, j).0.to_bits()))
            .collect();
        assert_eq!(decoded.len(), expected.len(), "duplicate masks");
        assert_eq!(decoded.into_iter().collect::<BTreeSet<_>>(), expected);
        for (i, j, _) in view.logical_pairs() {
            assert!(i < j && j < n);
        }
    }
}

#[test]
fn oracle_gates_live_coordinates_rebuild_box_reorder_resize_and_additivity() {
    for mode in MODES {
        let model = PackedForceModel::new(mode);
        let mut s = fixture(17);
        let mut nl = NeighborList::with_pair_scales(
            2.5,
            0.3,
            &[
                PairScale {
                    i: 0,
                    j: 1,
                    lj: 0.,
                    coulomb: 0.,
                },
                PairScale {
                    i: 1,
                    j: 2,
                    lj: 0.3,
                    coulomb: 1.,
                },
            ],
        )
        .unwrap();
        gate(&model, &mut s, &mut nl);
        let epoch = nl.rebuild_count();
        s.positions[0] += 0.01;
        gate(&model, &mut s, &mut nl);
        assert_eq!(
            nl.rebuild_count(),
            epoch,
            "current coordinate packing without rebuild"
        );
        s.positions[0] += 3.;
        nl.force_rebuild(&s);
        gate(&model, &mut s, &mut nl);
        s.box_size = SimBox::from_array([10., 11., 12.]);
        gate(&model, &mut s, &mut nl);
        let p0 = s.pos(0);
        let p4 = s.pos(4);
        s.set_pos(0, p4);
        s.set_pos(4, p0);
        gate(&model, &mut s, &mut nl);
        s = fixture(5);
        gate(&model, &mut s, &mut nl);
        // Different list, same rebuild count: epoch alone is not a cache key.
        nl = NeighborList::with_skin(2.5, 0.3);
        gate(&model, &mut s, &mut nl);
        s.resize(0);
        gate(&model, &mut s, &mut nl);
        s = fixture(3);
        gate(&model, &mut s, &mut nl);
    }
}

#[test]
fn cutoff_boundary_periodic_crossing_and_empty_clusters() {
    for mode in MODES {
        let model = PackedForceModel::new(mode);
        for r in [0., 2.5 - 1e-9, 2.5, 2.5 + 1e-9, 11.] {
            let mut s = fixture(2);
            s.set_pos(0, Vec3::new(0., 0., 0.));
            s.set_pos(1, Vec3::new(r, 0., 0.));
            let mut nl = NeighborList::with_skin(2.5, 0.3);
            gate(&model, &mut s, &mut nl);
            s.positions[0] = 11.99;
            s.positions[3] = 1.05;
            gate(&model, &mut s, &mut nl);
        }
    }
}

#[test]
fn unsupported_inputs_reject_explicitly() {
    for mode in MODES {
        let model = PackedForceModel::new(mode);
        for bad in 0..6 {
            let mut s = fixture(3);
            let mut nl = NeighborList::with_skin(2.5, 0.3);
            match bad {
                0 => {
                    s.box_size = SimBox::triclinic([
                        Vec3::new(12., 0., 0.),
                        Vec3::new(1., 12., 0.),
                        Vec3::new(0., 0., 12.),
                    ])
                }
                1 => s.species[0] = 1,
                2 => s.masses[0] = 2.,
                3 => s.positions[0] = f64::NAN,
                4 => nl = NeighborList::with_skin(2., 0.3),
                _ => s.forces.pop().map(|_| ()).unwrap(),
            }
            assert!(model
                .compute(&mut s, &mut Virial::default(), &mut nl)
                .is_err());
        }
    }
}

#[test]
fn full_step_matches_serial_reference_and_retains_timestep_convergence() {
    use modena_core::{Integrator, LennardJones, PairForceModel, VelocityVerlet};
    let mut initial = System::new(SimBox::cubic(20.));
    initial.resize(2);
    initial.set_pos(0, Vec3::new(9.4, 10., 10.));
    initial.set_pos(1, Vec3::new(10.6, 10., 10.));
    initial.set_vel(0, Vec3::new(0., 0.1, 0.));
    initial.set_vel(1, Vec3::new(0., -0.1, 0.));
    let run = |model: &dyn ForceModel, cutoff: f64, dt: f64| {
        let mut s = initial.clone();
        let mut nl = NeighborList::with_skin(cutoff, 0.3);
        let e0 = model
            .compute(&mut s, &mut Virial::default(), &mut nl)
            .unwrap()
            + s.kinetic_energy();
        let mut band = 0_f64;
        let mut vv = VelocityVerlet::new();
        for _ in 0..(2. / dt).round() as usize {
            let e = vv.step(&mut s, model, &mut nl, dt).unwrap();
            band = band.max((e.total - e0).abs());
        }
        (band, s)
    };
    // Reproduce the existing B06 observation exactly at its original cutoff4.2.
    // Packed modes use cutoff2.5; they must meet the same fixed 4 +/-0.3 criterion.
    let b06 = PairForceModel::new(LennardJones::reduced(4.2).unwrap());
    let mut models: Vec<(String, f64, Box<dyn ForceModel>)> =
        vec![("b06-original".into(), 4.2, Box::new(b06))];
    models.push((
        "interleaved-serial".into(),
        2.5,
        Box::new(PackedForceModel::serial_baseline()),
    ));
    models.extend(MODES.into_iter().map(|m| {
        (
            m.name().into(),
            2.5,
            Box::new(PackedForceModel::new(m)) as Box<dyn ForceModel>,
        )
    }));
    for (name, cutoff, model) in models {
        let mut bands = Vec::new();
        for dt in [0.004, 0.002, 0.001] {
            let (band, s) = run(&*model, cutoff, dt);
            bands.push(band);
            if cutoff == 2.5 {
                let (reference_band, reference) =
                    run(&PackedForceModel::serial_baseline(), 2.5, dt);
                near(band, reference_band, 1e-8);
                for (a, b) in s
                    .positions
                    .iter()
                    .chain(&s.velocities)
                    .chain(&s.forces)
                    .zip(
                        reference
                            .positions
                            .iter()
                            .chain(&reference.velocities)
                            .chain(&reference.forces),
                    )
                {
                    near(*a, *b, 1e-10);
                }
            }
        }
        eprintln!(
            "{name} timestep bands={bands:?} ratios={:?}",
            [bands[0] / bands[1], bands[1] / bands[2]]
        );
        assert!(bands.iter().all(|&x| x < 0.001));
        for adjacent in bands.windows(2) {
            assert!((adjacent[0] / adjacent[1] - 4.).abs() <= 0.3);
        }
    }
}

#[test]
fn large_pair_set_passes_oracle_in_four_thread_pool() {
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(4)
        .build()
        .unwrap();
    pool.install(|| {
        for mode in MODES {
            let mut s = modena_core::lattice::fcc(1.6796, [10, 5, 25], "Ar");
            for i in 0..s.n() {
                s.positions[3 * i] += 0.003 * (i as f64).sin();
            }
            let mut nl = NeighborList::with_skin(2.5, 0.3);
            gate(&PackedForceModel::new(mode), &mut s, &mut nl);
            assert!(nl.n_pairs() >= 100_000, "exercise production parallel path");
        }
    });
}

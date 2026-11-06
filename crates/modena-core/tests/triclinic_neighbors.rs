#![cfg(feature = "neighbor-diagnostics")]

use modena_core::snapshot::{evaluate_lj_snapshot, LjSnapshot, SnapshotUnits, SNAPSHOT_PROTOCOL};
use modena_core::{
    ForceModel, LennardJones, NeighborBuildFallbackReason, NeighborBuildStrategy, NeighborList,
    PairForceModel, PairScale, SimBox, System, Vec3, Virial,
};

fn cell(tilt: f64) -> SimBox {
    SimBox::triclinic([
        Vec3::new(10., 0., 0.),
        Vec3::new(tilt, 10., 0.),
        Vec3::new(0., 0., 10.),
    ])
}

fn system(points: &[Vec3]) -> System {
    let mut s = System::new(cell(3.));
    s.resize(points.len());
    for (i, &p) in points.iter().enumerate() {
        s.set_pos(i, p);
    }
    s
}

#[test]
fn scalar_order_is_lexicographic() {
    let s = system(&[
        Vec3::new(0., 0., 0.),
        Vec3::new(1., 0., 0.),
        Vec3::new(9.5, 0., 0.),
        Vec3::new(4., 0., 0.),
    ]);
    for strategy in [
        NeighborBuildStrategy::Auto,
        NeighborBuildStrategy::Scalar,
        NeighborBuildStrategy::FractionalStencil,
    ] {
        let mut nl = NeighborList::new(2.).with_triclinic_strategy_for_diagnostics(strategy);
        assert!(nl.triclinic_build_diagnostics().is_none());
        nl.build(&s);
        assert_eq!(nl.pairs(), &[(0, 1), (0, 2), (1, 2)]);
        let report = nl.triclinic_build_diagnostics().unwrap();
        assert_eq!(report.selected, NeighborBuildStrategy::Scalar);
        if strategy != NeighborBuildStrategy::FractionalStencil {
            assert_eq!(report.fallback, None);
            assert_eq!(report.requested_workspace_bytes, 0);
        } else {
            assert_eq!(
                report.fallback,
                Some(NeighborBuildFallbackReason::FractionalBoundary)
            );
            assert!(report.requested_workspace_bytes > 0);
        }
        assert_eq!(report.candidate_checks, 6);
    }
}

#[path = "support/periodic_oracle.rs"]
#[allow(dead_code)] // Shared graph wire-format helpers are not used by classical tests.
mod periodic_oracle;

fn interior_fixture(tilt: f64) -> periodic_oracle::Fixture {
    let basis = [[32., 0., 0.], [tilt, 32., 0.], [0., 0., 32.]];
    let mut positions = Vec::new();
    // Dyadic fractions sit strictly inside the m=12/12/12 bins (and m=10
    // bins with skin). Four separated two-atom clusters make cost eligible.
    for [x, y, z] in [[0., 0., 0.], [0.5, 0., 0.], [0., 0.5, 0.5], [0.5, 0.5, 0.5]] {
        for dx in [0., 1. / 32.] {
            let f = [x + 1. / 32. + dx, y + 3. / 64., z + 5. / 64.];
            positions.push(std::array::from_fn(|k| {
                (basis[0][k] * f[0] + basis[1][k] * f[1]) + basis[2][k] * f[2]
            }));
        }
    }
    periodic_oracle::Fixture {
        id: format!("interior-{tilt}"),
        cell: basis,
        atomic_numbers: vec![18; positions.len()],
        positions,
        pbc: [true; 3],
        cutoff: 2.5,
        periodic_self: false,
        cube: 16,
        official: false,
    }
}

fn forced(s: &System, cutoff: f64) -> NeighborList {
    let mut nl = NeighborList::new(cutoff)
        .with_triclinic_strategy_for_diagnostics(NeighborBuildStrategy::FractionalStencil);
    nl.build(s);
    nl
}

#[test]
fn stencil_matches_scalar_on_all_fixed_cells() {
    let mut cases = periodic_oracle::fixtures();
    cases.extend([-8., 0., 8., 24.].map(interior_fixture));
    let mut selected = 0;
    for f in cases {
        let s = periodic_oracle::system(&f);
        let stencil = forced(&s, f.cutoff);
        let mut scalar = NeighborList::new(f.cutoff)
            .with_triclinic_strategy_for_diagnostics(NeighborBuildStrategy::Scalar);
        scalar.build(&s);
        assert_eq!(stencil.pairs(), scalar.pairs(), "{}", f.id);
        assert!(stencil.pairs().windows(2).all(|p| p[0] < p[1]));
        eprintln!("{}: {:?}", f.id, stencil.triclinic_build_diagnostics());
        if f.id.starts_with("interior-") {
            let report = stencil.triclinic_build_diagnostics().unwrap();
            assert_eq!(
                report.selected,
                NeighborBuildStrategy::FractionalStencil,
                "{}: {report:?}",
                f.id
            );
            assert_eq!(report.fallback, None);
            assert!(report.candidate_checks < (s.n() * (s.n() - 1) / 2) as u64);
            assert!(
                report.requested_workspace_bytes > 0
                    && report.requested_workspace_bytes <= 134_217_728
            );
            selected += 1;
            let mut auto = NeighborList::new(f.cutoff);
            auto.build(&s);
            assert_eq!(
                auto.triclinic_build_diagnostics().unwrap().selected,
                NeighborBuildStrategy::Scalar
            );
        }
    }
    assert_eq!(selected, 4);
}

#[test]
fn stencil_matches_fixed_cube_membership() {
    let mut cases: Vec<_> = periodic_oracle::fixtures()
        .into_iter()
        .filter(|f| f.pbc == [true; 3])
        .collect();
    cases.extend([-8., 0., 8., 24.].map(interior_fixture));
    for f in cases {
        let s = periodic_oracle::system(&f);
        let nl = forced(&s, f.cutoff);
        let mut cube: Vec<_> = periodic_oracle::fixed_cube(&f)
            .into_iter()
            .filter_map(|(i, j, _)| (i < j).then_some((i, j)))
            .collect();
        cube.sort_unstable();
        cube.dedup(); // Independent multi-image oracle -> classical unordered membership.
        assert_eq!(nl.pairs(), cube, "{}", f.id);
        if f.id.starts_with("interior-") {
            assert_eq!(
                nl.triclinic_build_diagnostics().unwrap().selected,
                NeighborBuildStrategy::FractionalStencil
            );
        }
    }
}

#[test]
fn uncertain_fractional_face_falls_back() {
    let f = interior_fixture(8.);
    for p in [Vec3::zero(), Vec3::new(8., 0., 0.), Vec3::new(-32., 0., 0.)] {
        let mut s = periodic_oracle::system(&f);
        s.set_pos(0, p);
        let nl = forced(&s, f.cutoff);
        let report = nl.triclinic_build_diagnostics().unwrap();
        assert_eq!(report.selected, NeighborBuildStrategy::Scalar);
        assert_eq!(
            report.fallback,
            Some(NeighborBuildFallbackReason::FractionalBoundary)
        );
        let mut scalar = NeighborList::new(f.cutoff);
        scalar.build(&s);
        assert_eq!(nl.pairs(), scalar.pairs());
    }
}

#[test]
fn caps_fall_back_without_missing_pairs() {
    use NeighborBuildFallbackReason as F;
    let base = interior_fixture(8.);
    let mut huge = base.clone();
    huge.cell = huge.cell.map(|v| v.map(|x| x * 1000.));
    huge.positions = huge
        .positions
        .iter()
        .map(|v| v.map(|x| x * 1000.))
        .collect();
    let mut dense = base.clone();
    dense.positions.fill(base.positions[0]);
    let mut all = base.clone();
    all.cutoff = 20.;
    let mut image = base.clone();
    image.positions[0] = [1e20; 3];
    let mut bad = base.clone();
    bad.cell = [[1., 1., 0.], [1., 1.0_f64.next_up(), 0.], [0., 0., 1.]];
    // Coincident inputs avoid pathological scalar index enumeration after the
    // certificate rejects this ill-conditioned cell; fallback still executes.
    bad.positions.fill([0.; 3]);
    for (f, reason) in [
        (huge, F::GridCap),
        (dense, F::CandidateCost),
        (all, F::AllBins),
        (image, F::ImageDomain),
        (bad, F::UncertifiedGeometry),
    ] {
        let s = periodic_oracle::system(&f);
        let nl = forced(&s, f.cutoff);
        let report = nl.triclinic_build_diagnostics().unwrap();
        assert_eq!(report.selected, NeighborBuildStrategy::Scalar);
        assert_eq!(report.fallback, Some(reason), "{report:?}");
        eprintln!("fallback {reason:?}: {report:?}");
        let mut scalar = NeighborList::new(f.cutoff);
        scalar.build(&s);
        assert_eq!(nl.pairs(), scalar.pairs());
    }
}

#[test]
fn stencil_preserves_classical_force_energy_virial() {
    for tilt in [-8., 8., 24.] {
        for skin in [0., 0.5] {
            let f = interior_fixture(tilt);
            let mut scalar = periodic_oracle::system(&f);
            let mut stencil = periodic_oracle::system(&f);
            scalar.forces.fill(0.125);
            stencil.forces.fill(0.125);
            let scales = [
                PairScale {
                    i: 1,
                    j: 0,
                    lj: 0.,
                    coulomb: 0.,
                },
                PairScale {
                    i: 3,
                    j: 2,
                    lj: 1.3,
                    coulomb: 0.,
                },
                PairScale {
                    i: 5,
                    j: 4,
                    lj: 0.,
                    coulomb: 0.7,
                },
            ];
            let mut a = NeighborList::with_pair_scales(f.cutoff, skin, &scales)
                .unwrap()
                .with_triclinic_strategy_for_diagnostics(NeighborBuildStrategy::Scalar);
            let mut b = NeighborList::with_pair_scales(f.cutoff, skin, &scales)
                .unwrap()
                .with_triclinic_strategy_for_diagnostics(NeighborBuildStrategy::FractionalStencil);
            let model = PairForceModel::new(LennardJones::new(1., 1., f.cutoff).unwrap());
            let mut wa = Virial { w: 0.25 };
            let mut wb = wa;
            let ea = model.compute(&mut scalar, &mut wa, &mut a).unwrap();
            let eb = model.compute(&mut stencil, &mut wb, &mut b).unwrap();
            assert_eq!(
                b.triclinic_build_diagnostics().unwrap().selected,
                NeighborBuildStrategy::FractionalStencil
            );
            assert_eq!(a.pairs(), b.pairs());
            let mut list_fixture = f.clone();
            list_fixture.cutoff += skin;
            let mut cube: Vec<_> = periodic_oracle::fixed_cube(&list_fixture)
                .into_iter()
                .filter_map(|(i, j, _)| (i < j && (i, j) != (0, 1)).then_some((i, j)))
                .collect();
            cube.sort_unstable();
            cube.dedup();
            assert_eq!(b.pairs(), cube);

            for (x, y) in [ea, wa.w]
                .into_iter()
                .chain(scalar.forces)
                .zip([eb, wb.w].into_iter().chain(stencil.forces))
            {
                assert!(
                    (x - y).abs() <= 1e-10 + 1e-10 * x.abs().max(y.abs()),
                    "{x} != {y}"
                );
            }
        }
    }
}

#[test]
fn skin_half_threshold_and_box_change() {
    let f = interior_fixture(8.);
    let mut s = periodic_oracle::system(&f);
    let mut nl = NeighborList::with_skin(f.cutoff, 0.5)
        .with_triclinic_strategy_for_diagnostics(NeighborBuildStrategy::FractionalStencil);
    nl.build(&s);
    let report = nl.triclinic_build_diagnostics().unwrap();
    assert_eq!(report.selected, NeighborBuildStrategy::FractionalStencil);
    let r = s.pos(0);
    s.set_pos(0, r + Vec3::new(0.25, 0., 0.));
    nl.build(&s);
    assert_eq!(nl.rebuild_count(), 1);
    assert_eq!(nl.triclinic_build_diagnostics(), Some(report));
    s.set_pos(0, r + Vec3::new(0.250001, 0., 0.));
    nl.build(&s);
    assert_eq!(nl.rebuild_count(), 2);
    s.box_size = SimBox::triclinic(interior_fixture(-8.).cell.map(Vec3));
    nl.build(&s);
    assert_eq!(nl.rebuild_count(), 3);
    let mut scalar = NeighborList::with_skin(f.cutoff, 0.5);
    scalar.build(&s);
    assert_eq!(nl.pairs(), scalar.pairs());
}

#[test]
fn stencil_wraps_small_dimensions_unreduced_cells_and_atom_gauge() {
    let mut cases = Vec::new();
    for tilt in [-128., 128.] {
        let mut f = interior_fixture(tilt);
        // Move individual atoms by integral cell vectors. Fractional-bin
        // certification must keep their wrapped assignment and cube parity.
        for k in 0..3 {
            f.positions[1][k] += 2. * f.cell[0][k] - f.cell[1][k];
        }
        cases.push(f);
    }
    for short in [4., 6., 8.] {
        let mut f = interior_fixture(0.);
        f.cell = [[short, 0., 0.], [1., 32., 0.], [0., 0., 32.]];
        f.positions.clear();
        for y in [1.5, 9.5, 17.5, 25.5] {
            for dz in [0., 1.] {
                f.positions.push([0.375 + y / 32., y, 2.5 + dz]);
            }
        }
        // At cutoff2.5 these first-axis dimensions are 1, 2, and 3.
        cases.push(f);
    }
    for f in cases {
        let s = periodic_oracle::system(&f);
        let nl = forced(&s, f.cutoff);
        let report = nl.triclinic_build_diagnostics().unwrap();
        assert_eq!(
            report.selected,
            NeighborBuildStrategy::FractionalStencil,
            "{report:?}"
        );
        let mut scalar = NeighborList::new(f.cutoff);
        scalar.build(&s);
        let mut cube: Vec<_> = periodic_oracle::fixed_cube(&f)
            .into_iter()
            .filter_map(|(i, j, _)| (i < j).then_some((i, j)))
            .collect();
        cube.sort_unstable();
        cube.dedup();
        assert_eq!(nl.pairs(), scalar.pairs());
        assert_eq!(nl.pairs(), cube);
        assert_eq!(nl.pairs().len(), 4);
        eprintln!("extended cell {:?}: {report:?}", f.cell);
    }
}

#[test]
fn stencil_skin_cache_survives_a_periodic_seam() {
    let f = interior_fixture(8.);
    let mut s = periodic_oracle::system(&f);
    let r = s.pos(0);
    s.set_pos(0, Vec3::new(31.75 + 8. * 3. / 64., r.0[1], r.0[2]));
    let mut nl = NeighborList::with_skin(2.5, 1.)
        .with_triclinic_strategy_for_diagnostics(NeighborBuildStrategy::FractionalStencil);
    nl.build(&s);
    assert_eq!(
        nl.triclinic_build_diagnostics().unwrap().selected,
        NeighborBuildStrategy::FractionalStencil
    );
    let original = nl.pairs().to_vec();
    let mut list_fixture = f.clone();
    list_fixture.cutoff = 3.5;
    list_fixture.positions = (0..s.n()).map(|i| s.pos(i).0).collect();
    let mut cube: Vec<_> = periodic_oracle::fixed_cube(&list_fixture)
        .into_iter()
        .filter_map(|(i, j, _)| (i < j).then_some((i, j)))
        .collect();
    cube.sort_unstable();
    cube.dedup();
    assert_eq!(original, cube);

    assert!(original.contains(&(0, 1)));
    s.set_pos(0, s.pos(0) + Vec3::new(0.4, 0., 0.));
    s.wrap_positions();
    nl.build(&s);
    assert_eq!(nl.rebuild_count(), 1);
    let mut scalar = NeighborList::new(2.5);
    scalar.build(&s);
    assert!(scalar.pairs().iter().all(|p| original.contains(p)));
}

#[test]
fn zero_skin_rebuilds_every_build() {
    let mut s = system(&[Vec3::zero(), Vec3::new(1., 0., 0.)]);
    let mut nl = NeighborList::new(2.)
        .with_triclinic_strategy_for_diagnostics(NeighborBuildStrategy::Scalar);
    nl.build(&s);
    nl.build(&s);
    assert_eq!(nl.rebuild_count(), 2);
    s.set_pos(1, Vec3::new(3., 0., 0.));
    nl.build(&s);
    assert_eq!(nl.rebuild_count(), 3);
    assert!(nl.pairs().is_empty());
}

#[test]
fn diagnostics_follow_successful_triclinic_cache_lifecycle() {
    let mut s = system(&[Vec3::zero(), Vec3::new(1., 0., 0.)]);
    let mut nl = NeighborList::with_skin(2., 0.3);
    nl.build(&s);
    let report = nl.triclinic_build_diagnostics().unwrap();
    nl.build(&s);
    assert_eq!(nl.rebuild_count(), 1);
    assert_eq!(nl.triclinic_build_diagnostics(), Some(report));
    nl = nl.with_triclinic_strategy_for_diagnostics(NeighborBuildStrategy::Scalar);
    assert!(nl.triclinic_build_diagnostics().is_none());
    nl.build(&s);
    assert_eq!(nl.rebuild_count(), 2);
    assert_eq!(nl.triclinic_build_diagnostics(), Some(report));
    s.box_size = SimBox::cubic(10.);
    nl.build(&s);
    assert!(nl.triclinic_build_diagnostics().is_none());
    s.box_size = cell(3.);
    nl.build(&s);
    assert!(nl.triclinic_build_diagnostics().is_some());
    s.resize(0);
    nl.build(&s);
    assert!(nl.triclinic_build_diagnostics().is_none());
}

fn snapshot(n: usize) -> LjSnapshot {
    LjSnapshot {
        protocol: SNAPSHOT_PROTOCOL.into(),
        units: SnapshotUnits::default(),
        atomic_numbers: vec![18; n],
        positions_angstrom: (0..n).map(|i| [i as f64, 0., 0.]).collect(),
        cell_angstrom: [[10., 0., 0.], [0., 10., 0.], [0., 0., 10.]],
        pbc: [true; 3],
        sigma_angstrom: 1.,
        epsilon_ev: 1.,
        cutoff_angstrom: 3.,
        skin_angstrom: 0.,
    }
}

#[test]
fn bounded_zero_skin_matches_scalar() {
    // The public snapshot protocol is orthorhombic. Private unit tests cover
    // the actual triclinic bounded builder without adding a public test hook.
    let input = snapshot(3);
    let bounded = evaluate_lj_snapshot(&input).unwrap();
    let mut s = System::new(SimBox::cubic(10.));
    s.resize(3);
    for (i, &p) in input.positions_angstrom.iter().enumerate() {
        s.set_pos(i, Vec3(p));
    }
    let mut nl = NeighborList::new(3.);
    let mut virial = Virial::default();
    let model = PairForceModel::new(LennardJones::new(1., 1., 3.).unwrap());
    let energy = model.compute(&mut s, &mut virial, &mut nl).unwrap();
    assert_eq!(energy, bounded.energy_ev);
    assert_eq!(nl.pairs(), &[(0, 1), (0, 2), (1, 2)]);
    assert_eq!(s.forces, bounded.forces_ev_per_angstrom.concat());
}

#[test]
fn bounded_pair_overflow_still_errors() {
    // 1415 atoms have 1,000,405 pairs, exceeding the public snapshot cap.
    // Coincident positions also prove cap checking precedes force evaluation.
    let mut input = snapshot(1415);
    input.positions_angstrom.fill([0.; 3]);
    let error = evaluate_lj_snapshot(&input).unwrap_err();
    assert!(error.to_string().contains("1000000 neighbor pairs"));
}

#[test]
fn same_bounding_box_changed_basis_rebuilds() {
    let mut s = system(&[Vec3::zero(), Vec3::new(3., 9., 0.)]);
    let mut nl = NeighborList::with_skin(1.1, 0.2);
    nl.build(&s);
    assert_eq!(nl.pairs(), &[(0, 1)]);
    let bounds = s.box_size.size;
    s.box_size = cell(-3.);
    assert_eq!(s.box_size.size, bounds);
    nl.build(&s);
    assert_eq!(nl.rebuild_count(), 2);
    assert!(nl.pairs().is_empty());
}

#[test]
fn periodic_crossing_skin_is_complete() {
    let mut s = system(&[Vec3::new(0.05, 0., 0.), Vec3::new(1.7, 0., 0.)]);
    let mut nl = NeighborList::with_skin(1.5, 0.6);
    nl.build(&s);
    assert_eq!(nl.pairs(), &[(0, 1)]); // Included by skin, outside force cutoff.
    s.set_pos(0, Vec3::new(-0.05, 0., 0.));
    s.set_pos(1, Vec3::new(1.41, 0., 0.));
    s.wrap_positions();
    nl.build(&s);
    assert_eq!(nl.rebuild_count(), 1);
    assert_eq!(nl.pairs(), &[(0, 1)]);
    s.set_pos(1, Vec3::new(1.3, 0., 0.));
    nl.build(&s);
    assert_eq!(nl.rebuild_count(), 2);
    assert_eq!(nl.pairs(), &[(0, 1)]);
}

#[test]
fn pair_scale_zero_zero_only_exclusion() {
    let s = system(&[Vec3::zero(), Vec3::new(1., 0., 0.), Vec3::new(2., 0., 0.)]);
    let scales = [
        PairScale {
            i: 1,
            j: 0,
            lj: 0.,
            coulomb: 0.,
        },
        PairScale {
            i: 2,
            j: 0,
            lj: 0.,
            coulomb: 0.7,
        },
        PairScale {
            i: 1,
            j: 2,
            lj: 1.3,
            coulomb: 0.,
        },
    ];
    let mut nl = NeighborList::with_pair_scales(3., 0., &scales)
        .unwrap()
        .with_triclinic_strategy_for_diagnostics(NeighborBuildStrategy::Scalar);
    nl.build(&s);
    assert_eq!(nl.pairs(), &[(0, 2), (1, 2)]);
    assert_eq!(nl.pair_scale(0, 2), (0., 0.7));
    assert_eq!(nl.pair_scale(2, 1), (1.3, 0.));
    assert_eq!(
        nl.triclinic_build_diagnostics().unwrap().candidate_checks,
        2
    );
}

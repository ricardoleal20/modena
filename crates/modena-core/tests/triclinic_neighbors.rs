#![cfg(feature = "neighbor-diagnostics")]

use modena_core::snapshot::{evaluate_lj_snapshot, LjSnapshot, SnapshotUnits, SNAPSHOT_PROTOCOL};
use modena_core::{
    ForceModel, LennardJones, NeighborBuildStrategy, NeighborList, PairForceModel, PairScale,
    SimBox, System, Vec3, Virial,
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
        // Task 5 stages the request seam; Task 6 supplies stencil execution.
        assert_eq!(report.selected, NeighborBuildStrategy::Scalar);
        assert_eq!(report.fallback, None);
        assert_eq!(report.candidate_checks, 6);
        assert_eq!(report.requested_workspace_bytes, 0);
    }
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

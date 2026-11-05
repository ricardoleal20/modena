use modena_core::snapshot::{evaluate_lj_snapshot, LjSnapshot, SnapshotUnits, SNAPSHOT_PROTOCOL};

fn pair() -> LjSnapshot {
    LjSnapshot {
        protocol: SNAPSHOT_PROTOCOL.into(),
        units: SnapshotUnits::default(),
        atomic_numbers: vec![18, 18],
        positions_angstrom: vec![[0.; 3], [4., 0., 0.]],
        cell_angstrom: [[20., 0., 0.], [0., 20., 0.], [0., 0., 20.]],
        pbc: [true; 3],
        sigma_angstrom: 3.4,
        epsilon_ev: 0.0103,
        cutoff_angstrom: 8.5,
        skin_angstrom: 0.3,
    }
}

#[test]
fn native_snapshot_is_dimensionally_consistent_and_has_no_mass_dependence() {
    let mut input = pair();
    let output = evaluate_lj_snapshot(&input).unwrap();
    // Hand-evaluated independent two-particle SF-LJ value (eV).
    assert!((output.energy_ev - (-0.008978481029880884)).abs() < 1e-14);
    assert!((output.forces_ev_per_angstrom[0][0] - 0.005608613407781233).abs() < 1e-14);
    assert_eq!(output.forces_ev_per_angstrom.len(), 2);
    assert_eq!(output.forces_ev_per_angstrom[0][1..], [0., 0.]);
    input.atomic_numbers = vec![118, 118];
    let same = evaluate_lj_snapshot(&input).unwrap();
    assert_eq!(same.energy_ev, output.energy_ev);
    assert_eq!(same.forces_ev_per_angstrom, output.forces_ev_per_angstrom);
    input.epsilon_ev = 0.;
    let zero = evaluate_lj_snapshot(&input).unwrap();
    assert_eq!(zero.energy_ev, 0.);
    assert!(zero
        .forces_ev_per_angstrom
        .iter()
        .flatten()
        .all(|&v| v == 0.));
}

#[test]
fn native_snapshot_rejects_nonfinite_direct_rust_inputs_and_resource_overflow() {
    for v in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        for field in 0..6 {
            let mut input = pair();
            match field {
                0 => input.sigma_angstrom = v,
                1 => input.epsilon_ev = v,
                2 => input.cutoff_angstrom = v,
                3 => input.skin_angstrom = v,
                4 => input.cell_angstrom[0][0] = v,
                _ => input.positions_angstrom[0][0] = v,
            }
            assert!(
                evaluate_lj_snapshot(&input).is_err(),
                "accepted field {field}: {v}"
            );
        }
    }
    let mut huge = pair();
    huge.cell_angstrom[0][0] = 1e200;
    assert!(evaluate_lj_snapshot(&huge).is_err());
    let mut tiny = pair();
    tiny.cutoff_angstrom = 1e-200;
    tiny.skin_angstrom = 0.;
    assert!(evaluate_lj_snapshot(&tiny).is_err());
    let mut overflow = pair();
    overflow.epsilon_ev = f64::MAX;
    assert!(evaluate_lj_snapshot(&overflow).is_err());
    let mut large = pair();
    large.atomic_numbers.resize(100001, 18);
    large.positions_angstrom.resize(100001, [0.; 3]);
    assert!(evaluate_lj_snapshot(&large).is_err());
}

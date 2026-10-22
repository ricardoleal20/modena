use modena_core::thermostat::{NoneThermostat, Thermostat, VelocityRescaling};
use modena_core::{
    Constraint, Constraints, ForceModel, Integrator, LangevinBAOAB, NeighborList, NoseHooverChain,
    Result, SimBox, System, Vec3, VelocityVerlet, Virial,
};

struct ZeroForce;
impl ForceModel for ZeroForce {
    fn compute(&self, _: &mut System, _: &mut Virial, _: &mut NeighborList) -> Result<f64> {
        Ok(0.0)
    }
}
fn gas() -> System {
    let mut s = System::new(SimBox::cubic(100.0));
    s.resize(4);
    s.velocities.fill(0.7);
    s
}
fn json_roundtrip(value: serde_json::Value) -> serde_json::Value {
    serde_json::from_str(&serde_json::to_string(&value).unwrap()).unwrap()
}

#[test]
fn csvr_json_state_reproduces_stochastic_continuation() {
    let mut system = gas();
    let mut thermostat = VelocityRescaling::new(1.0, 0.37, 42);
    for _ in 0..13 {
        thermostat.apply(&mut system, 0.01);
    }
    let state = json_roundtrip(thermostat.checkpoint_state().unwrap());
    let mut split = system.clone();
    let mut restored = VelocityRescaling::new(1.0, 0.37, 42);
    restored.restore_state(&state).unwrap();
    for _ in 0..17 {
        thermostat.apply(&mut system, 0.01);
        restored.apply(&mut split, 0.01);
    }
    assert_eq!(
        split.velocities, system.velocities,
        "CSVR RNG continuation differs"
    );
}

#[test]
fn langevin_json_state_reproduces_stochastic_continuation() {
    let mut system = gas();
    let mut integrator = LangevinBAOAB::new(1.23, 0.7, 42);
    let mut neighbors = NeighborList::new(2.5);
    for _ in 0..13 {
        integrator
            .step(&mut system, &ZeroForce, &mut neighbors, 0.01)
            .unwrap();
    }
    let state = json_roundtrip(integrator.checkpoint_state().unwrap());
    let mut split = system.clone();
    let mut restored = LangevinBAOAB::new(1.23, 0.7, 42);
    restored.restore_state(&state).unwrap();
    let mut split_neighbors = NeighborList::new(2.5);
    for _ in 0..17 {
        integrator
            .step(&mut system, &ZeroForce, &mut neighbors, 0.01)
            .unwrap();
        restored
            .step(&mut split, &ZeroForce, &mut split_neighbors, 0.01)
            .unwrap();
    }
    assert_eq!(
        split.velocities, system.velocities,
        "Langevin RNG continuation differs"
    );
    assert_eq!(split.positions, system.positions);
}

#[test]
fn nhc_json_state_preserves_extended_momenta_and_positions() {
    let mut system = gas();
    let mut integrator = NoseHooverChain::new(1.0, 0.4, 4);
    let mut neighbors = NeighborList::new(2.5);
    for _ in 0..13 {
        integrator
            .step(&mut system, &ZeroForce, &mut neighbors, 0.01)
            .unwrap();
    }
    let state = json_roundtrip(integrator.checkpoint_state().unwrap());
    let mut split = system.clone();
    let mut restored = NoseHooverChain::new(1.0, 0.4, 4);
    restored.restore_state(&state).unwrap();
    let mut split_neighbors = NeighborList::new(2.5);
    for _ in 0..17 {
        integrator
            .step(&mut system, &ZeroForce, &mut neighbors, 0.01)
            .unwrap();
        restored
            .step(&mut split, &ZeroForce, &mut split_neighbors, 0.01)
            .unwrap();
    }
    assert_eq!(
        split.velocities, system.velocities,
        "NHC continuation differs"
    );
    assert_eq!(
        restored.checkpoint_state().unwrap(),
        integrator.checkpoint_state().unwrap()
    );
}

#[test]
fn state_restore_rejects_different_kind_or_parameters() {
    let csvr = VelocityRescaling::new(1.0, 0.4, 42);
    let state = csvr.checkpoint_state().unwrap();
    assert!(VelocityRescaling::new(2.0, 0.4, 42)
        .restore_state(&state)
        .is_err());
    assert!(LangevinBAOAB::new(1.0, 1.0, 42)
        .restore_state(&state)
        .is_err());
    let nhc = NoseHooverChain::new(1.0, 0.4, 4);
    assert!(NoseHooverChain::new(1.0, 0.4, 3)
        .restore_state(&nhc.checkpoint_state().unwrap())
        .is_err());
}

#[test]
fn stateless_state_rejects_non_null() {
    assert_eq!(
        VelocityVerlet::new().checkpoint_state().unwrap(),
        serde_json::Value::Null
    );
    assert!(VelocityVerlet::new()
        .restore_state(&serde_json::json!({"bad":1}))
        .is_err());
    assert!(NoneThermostat
        .restore_state(&serde_json::json!({"bad":1}))
        .is_err());
}

#[test]
fn physical_dof_subtract_constraints_and_conserved_com_only() {
    let mut s = gas();
    assert_eq!(s.dof(), 12);
    s.set_constraint_count(4).unwrap();
    s.set_conserved_com_dof(3).unwrap();
    assert_eq!(s.dof(), 5);
    assert_eq!(s.constraint_count(), 4);
    assert_eq!(s.conserved_com_dof(), 3);
    assert!(s.set_constraint_count(20).is_err());
    assert!(s.set_conserved_com_dof(4).is_err());
}

#[test]
fn langevin_records_constraint_count_and_keeps_stochastic_com_dof() {
    let mut s = gas();
    s.set_conserved_com_dof(3).unwrap();
    s.set_pos(0, Vec3::new(30.0, 30.0, 30.0));
    s.set_pos(1, Vec3::new(31.0, 30.0, 30.0));
    let mut cons = Constraints::new();
    cons.bonds.push(Constraint { i: 0, j: 1, d: 1.0 });
    let mut integrator = LangevinBAOAB::with_constraints(1.0, 1.0, 42, cons);
    integrator
        .step(&mut s, &ZeroForce, &mut NeighborList::new(2.5), 0.01)
        .unwrap();
    assert_eq!(s.dof(), 11);
    assert_eq!(s.conserved_com_dof(), 0);
}

#[test]
fn csvr_one_dof_updates_without_gamma_shape_zero() {
    let mut s = System::new(SimBox::cubic(10.0));
    s.resize(1);
    s.set_constraint_count(2).unwrap();
    s.set_vel(0, Vec3::new(1.0, 0.0, 0.0));
    let mut thermostat = VelocityRescaling::new(1.0, 1.0, 42);
    let mut sum = 0.0;
    let mut sum2 = 0.0;
    for step in 0..41_000 {
        thermostat.apply(&mut s, 2.0);
        if step >= 1000 {
            let k = s.kinetic_energy();
            sum += k;
            sum2 += k * k;
        }
    }
    let mean = sum / 40_000.0;
    let variance = sum2 / 40_000.0 - mean * mean;
    assert!((mean - 0.5).abs() < 0.015, "one DOF canonical mean={mean}");
    assert!(
        (variance - 0.5).abs() < 0.04,
        "one DOF canonical variance={variance}"
    );
}

#[test]
fn state_restore_rejects_constraint_mismatch_and_malformed_internal_state() {
    let state = LangevinBAOAB::new(1.0, 1.0, 42).checkpoint_state().unwrap();
    let mut cons = Constraints::new();
    cons.bonds.push(Constraint { i: 0, j: 1, d: 1.0 });
    assert!(LangevinBAOAB::with_constraints(1.0, 1.0, 42, cons)
        .restore_state(&state)
        .is_err());
    let mut corrupted = state;
    corrupted["rng"]["word_pos"] = serde_json::json!("295147905179352825856"); // 2^68, outside stream
    assert!(LangevinBAOAB::new(1.0, 1.0, 42)
        .restore_state(&corrupted)
        .is_err());
    let mut nhc = NoseHooverChain::new(1.0, 0.4, 4);
    let before = nhc.checkpoint_state().unwrap();
    let mut corrupted = before.clone();
    corrupted["xi"] = serde_json::json!([0.0, 0.0]);
    assert!(nhc.restore_state(&corrupted).is_err());
    assert_eq!(
        nhc.checkpoint_state().unwrap(),
        before,
        "failed restore must preserve existing state"
    );
}

#[test]
fn checkpoint_rejects_nonfinite_constraint_parameters() {
    let mut constraints = Constraints::new();
    constraints.bonds.push(Constraint {
        i: 0,
        j: 1,
        d: f64::NAN,
    });
    assert!(
        LangevinBAOAB::with_constraints(1.0, 1.0, 42, constraints.clone())
            .checkpoint_state()
            .is_err()
    );
    assert!(NoseHooverChain::new(1.0, 0.4, 4)
        .with_constraints(constraints)
        .checkpoint_state()
        .is_err());
}

use modena_core::constraints::RigidWater;
use modena_core::{Constraint, Constraints, SimBox, System, Vec3};

fn dimer() -> (System, Constraints) {
    let mut s = System::new(SimBox::cubic(20.0));
    s.resize(2);
    s.set_pos(1, Vec3::new(1.2, 0.0, 0.0));
    let mut c = Constraints::new();
    c.bonds.push(Constraint { i: 0, j: 1, d: 1.0 });
    (s, c)
}
#[test]
fn shake_returns_error_when_iteration_budget_is_exhausted() {
    let (mut s, c) = dimer();
    let reference = s.positions.clone();
    assert!(c.shake(&mut s, &reference, 1e-12, 0).is_err());
}
#[test]
fn rattle_returns_error_when_iteration_budget_is_exhausted() {
    let (mut s, c) = dimer();
    s.set_vel(1, Vec3::new(1.0, 0.0, 0.0));
    assert!(c.rattle(&mut s, 1e-12, 0).is_err());
}
#[test]
fn projections_reject_invalid_indices_and_masses() {
    let (mut s, mut c) = dimer();
    let reference = s.positions.clone();
    c.bonds[0].j = 4;
    assert!(c.shake(&mut s, &reference, 1e-12, 100).is_err());
    c.bonds[0].j = 1;
    s.masses[1] = 0.0;
    assert!(c.shake(&mut s, &reference, 1e-12, 100).is_err());
    assert!(c.rattle(&mut s, 1e-12, 100).is_err());
}
#[test]
fn settle_errors_on_singular_or_invalid_water() {
    let mut s = System::new(SimBox::cubic(20.0));
    s.resize(3);
    let mut c = Constraints::new();
    c.waters.push(RigidWater {
        first: 0,
        d_oh: 1.0,
        d_hh: 1.6,
    });
    assert!(c.settle(&mut s, 1e-12).is_err());
    c.waters[0].first = 2;
    assert!(c.settle_velocities(&mut s, 1e-12).is_err());
}
#[test]
fn settle_projects_selected_water_and_leaves_unrelated_atoms_untouched() {
    let mut s = System::new(SimBox::cubic(20.0));
    s.resize(5);
    s.set_pos(0, Vec3::new(9.0, 9.0, 9.0));
    s.set_pos(4, Vec3::new(8.0, 8.0, 8.0));
    s.set_pos(1, Vec3::new(3.0, 3.0, 3.0));
    s.set_pos(2, Vec3::new(3.7, 3.7, 3.0));
    s.set_pos(3, Vec3::new(3.7, 2.3, 3.0));
    let mut c = Constraints::new();
    c.waters.push(RigidWater {
        first: 1,
        d_oh: 1.0,
        d_hh: 1.6,
    });
    c.settle(&mut s, 1e-12).unwrap();
    assert!(c.max_residual(&s) < 1e-12);
    assert_eq!(s.pos(0), Vec3::new(9.0, 9.0, 9.0));
    assert_eq!(s.pos(4), Vec3::new(8.0, 8.0, 8.0));
    s.set_vel(0, Vec3::new(1.0, 2.0, 3.0));
    s.set_vel(2, Vec3::new(0.3, 0.2, 0.1));
    c.settle_velocities(&mut s, 1e-12).unwrap();
    assert!(c.max_velocity_residual(&s) < 1e-12);
    assert_eq!(s.vel(0), Vec3::new(1.0, 2.0, 3.0));
}

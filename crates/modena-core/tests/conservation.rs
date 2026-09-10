//! Integration test: NVE energy & momentum conservation.
//!
//! The canonical correctness check for an MD integrator. A small LJ solid is
//! evolved under velocity-Verlet; total energy must drift by < 0.1% and linear
//! momentum must stay near zero. This is the test that makes the engine
//! defensible in a research paper.

use modena_core::{
    lattice, ForceModel, Integrator, LennardJones, NeighborList, PairForceModel, PairPotential,
    SimBox, System, Vec3, VelocityVerlet, Virial,
};

fn build_lj() -> (Box<dyn ForceModel>, f64) {
    let lj = LennardJones::new(1.0, 1.0, 2.5).unwrap();
    let cutoff = lj.cutoff();
    (Box::new(PairForceModel::new(lj)), cutoff)
}

#[test]
fn pair_energy_and_forces_are_correct() {
    // Two atoms at r = 1.5 in a 10³ box: PE must equal the single shifted-force
    // pair energy, and forces must be equal and opposite. Catches overcounting.
    let (model, cutoff) = build_lj();
    let mut sys = System::new(SimBox::cubic(10.0));
    sys.resize(2);
    sys.positions[0] = Vec3::new(0.0, 0.0, 0.0);
    sys.positions[1] = Vec3::new(1.5, 0.0, 0.0);

    let mut nl = NeighborList::new(cutoff);
    let mut virial = Virial::default();
    let pe = model.compute(&mut sys, &mut virial, &mut nl).unwrap();

    // Expected single-pair energy (the kernel halves the double-ended sum).
    let r = 1.5_f64;
    let (u_expected, _) = {
        let lj = LennardJones::new(1.0, 1.0, 2.5).unwrap();
        lj.energy_force(r * r, 0, 0)
    };
    assert!(
        (pe - u_expected).abs() < 1e-9,
        "PE={pe}, expected {u_expected}"
    );

    // Newton's third law: equal and opposite.
    let f0 = sys.forces[0];
    let f1 = sys.forces[1];
    assert!(
        (f0 + f1).norm() < 1e-9,
        "forces not equal/opposite: {f0:?} {f1:?}"
    );
    // Repulsive at r=1.5 (>r_min≈1.122? no, 1.5>r_min → attractive). Force on 0
    // should point toward +x (toward atom 1) → x component negative? At r=1.5
    // (>1.122), force is attractive: atom 0 pulled toward atom 1 → +x.
    assert!(
        f0.0[0] > 0.0,
        "expected attractive +x force on atom 0, got {f0:?}"
    );
}

#[test]
fn nve_energy_is_conserved() {
    let (model, cutoff) = build_lj();
    let mut system = lattice::fcc(1.6796, [2, 2, 2], "Ar"); // 32 atoms
    assert_eq!(system.n(), 32);

    // Energetic, deterministic initial velocities (amplitude large enough that
    // atoms actually translate and cross neighbor-list cell boundaries over the
    // run — this is what catches cutoff/force-discontinuity drift).
    for i in 0..system.n() {
        let fi = i as f64;
        let v = modena_core::Vec3::new(
            0.5 * fi.sin(),
            0.5 * (1.3 * fi).cos(),
            0.5 * (0.7 * fi).sin(),
        );
        system.velocities[i] = v;
    }
    // Remove COM drift.
    let mut com = modena_core::Vec3::zero();
    for v in &system.velocities {
        com += *v;
    }
    com *= 1.0 / system.n() as f64;
    for v in system.velocities.iter_mut() {
        *v -= com;
    }

    let mut neighbors = NeighborList::new(cutoff);
    let mut integrator = VelocityVerlet::new();

    // Initial forces so F(t) is cached for the first half-kick.
    use modena_core::Virial;
    let mut virial = Virial::default();
    let _ = model.compute(&mut system, &mut virial, &mut neighbors);

    let dt = 0.002;
    let nsteps = 2000;

    // Run and capture first & last total energy.
    let mut e_first: Option<f64> = None;
    let mut e_last: f64 = 0.0;
    for _ in 0..nsteps {
        let e = integrator
            .step(&mut system, &*model, &mut neighbors, dt)
            .expect("step");
        if e_first.is_none() {
            e_first = Some(e.total);
        }
        e_last = e.total;
    }

    let e0 = e_first.unwrap();
    let drift = (e_last - e0).abs() / e0.abs().max(1e-12);
    assert!(
        drift < 1e-3,
        "energy drift too large: {drift:.3e} (E0={e0:.6}, E_last={e_last:.6})"
    );

    // Momentum should remain ~zero (no external forces, PBC).
    let p = system.total_momentum();
    assert!(p.norm() < 1e-6, "linear momentum grew: |p| = {}", p.norm());
}

#[test]
fn nvt_langevin_stabilizes_temperature() {
    // Sanity check that the Langevin integrator thermalizes toward the target.
    use modena_core::{LangevinBAOAB, Virial};

    let (model, cutoff) = build_lj();
    let mut system = lattice::fcc(1.6796, [2, 2, 2], "Ar");
    // Hot start.
    for v in system.velocities.iter_mut() {
        *v = modena_core::Vec3::new(2.0, 2.0, 2.0);
    }
    let mut com = modena_core::Vec3::zero();
    for v in &system.velocities {
        com += *v;
    }
    com *= 1.0 / system.n() as f64;
    for v in system.velocities.iter_mut() {
        *v -= com;
    }

    let mut neighbors = NeighborList::new(cutoff);
    let mut integrator = LangevinBAOAB::new(1.0, 1.0, 7);
    let mut virial = Virial::default();
    let _ = model.compute(&mut system, &mut virial, &mut neighbors);

    let dt = 0.005;
    let mut last_t = system.temperature();
    for _ in 0..4000 {
        let _e = integrator
            .step(&mut system, &*model, &mut neighbors, dt)
            .expect("step");
        last_t = system.temperature();
    }
    // After many steps at gamma=1, T=1, the system should be within ~30% of 1.0.
    // (Loose bound — this is a stability smoke test, not a fluctuation test.)
    assert!(
        (last_t - 1.0).abs() < 0.5,
        "Langevin failed to thermalize: T = {last_t:.4}"
    );
}

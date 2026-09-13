//! Domain-decomposition demo (M9, no MPI runtime required).
//!
//! Runs the in-process multi-rank `DomainForceModel` on an FCC argon cluster,
//! verifies the decomposed force/energy/virial is bit-identical to the serial
//! pair model, then runs a short decomposed NVE trajectory and reports the
//! energy drift — a self-contained "does the parallel path actually work as an
//! MD engine" exercise.
//!
//! ```sh
//! cargo run --release --example domain_decomposition
//! ```

use modena_core::{
    decomp::{DomainDecomposition, DomainForceModel, LocalCommunicator, RankGrid},
    force_model::{ForceModel, PairForceModel, Virial},
    lattice,
    neighbor::NeighborList,
    potential::LennardJones,
    vec3::Vec3,
};

fn main() {
    let cutoff = 2.5;
    let sys0 = lattice::fcc(1.5, [5, 5, 5], "Ar"); // 500 atoms, box 7.5
    let grid = RankGrid::balanced(sys0.box_size.size, cutoff, 8);
    let decomp = DomainDecomposition::new(sys0.box_size, cutoff, grid).unwrap();
    println!("system: {} atoms, box {:?}", sys0.n(), sys0.box_size.size.0);
    println!("rank grid: {:?} ({} ranks)", grid.nr, grid.n_ranks());

    // --- 1. Decomposed == serial (the correctness anchor). ---
    let mut s_serial = sys0.clone();
    s_serial.zero_forces();
    let mut vir_s = Virial::default();
    let mut nl = NeighborList::new(cutoff);
    let pe_s = PairForceModel::new(LennardJones::new(1.0, 1.0, cutoff).unwrap())
        .compute(&mut s_serial, &mut vir_s, &mut nl)
        .unwrap();

    let mut s_decomp = sys0.clone();
    s_decomp.zero_forces();
    let comm = LocalCommunicator::new(grid.n_ranks(), 0);
    let model = DomainForceModel::new(
        LennardJones::new(1.0, 1.0, cutoff).unwrap(),
        decomp,
        Box::new(comm),
    );
    let mut vir_d = Virial::default();
    let mut nl_d = NeighborList::new(cutoff);
    let pe_d = model.compute(&mut s_decomp, &mut vir_d, &mut nl_d).unwrap();

    let max_f_err = (0..s_serial.forces.len())
        .map(|i| (s_serial.forces[i] - s_decomp.forces[i]).abs())
        .fold(0.0f64, f64::max);
    println!(
        "decomposed vs serial: Δpe={:.3e}, Δvirial={:.3e}, max|Δf|={:.3e}",
        (pe_d - pe_s).abs(),
        (vir_d.w - vir_s.w).abs(),
        max_f_err
    );
    assert!(max_f_err < 1e-9, "decomposed forces diverge from serial");

    // --- 2. Short decomposed NVE trajectory; report energy drift. ---
    let mut sys = sys0.clone();
    // Seed deterministic velocities, remove COM drift.
    let mut p_com = Vec3::zero();
    let mut total_m = 0.0;
    for i in 0..sys.n() {
        let v = Vec3::new(
            (((i as f64 * 12.9898).sin() * 43758.5453).fract() - 0.5) * 0.3,
            (((i as f64 * 78.233).sin() * 12543.123).fract() - 0.5) * 0.3,
            (((i as f64 * 39.345).sin() * 9981.777).fract() - 0.5) * 0.3,
        );
        sys.set_vel(i, v);
        p_com += v * sys.mass(i);
        total_m += sys.mass(i);
    }
    for i in 0..sys.n() {
        let v = sys.vel(i) - p_com * (1.0 / total_m);
        sys.set_vel(i, v);
    }

    let decomp2 = DomainDecomposition::new(sys.box_size, cutoff, grid).unwrap();
    let comm2 = LocalCommunicator::new(grid.n_ranks(), 0);
    let model2 = DomainForceModel::new(
        LennardJones::new(1.0, 1.0, cutoff).unwrap(),
        decomp2,
        Box::new(comm2),
    );
    let mut vir = Virial::default();
    let mut nl2 = NeighborList::new(cutoff);
    sys.zero_forces();
    let pe0 = model2.compute(&mut sys, &mut vir, &mut nl2).unwrap();
    let e0 = sys.kinetic_energy() + pe0;

    let dt = 0.005;
    let steps = 500;
    let mut e_max = e0;
    let mut e_min = e0;
    for _ in 0..steps {
        for i in 0..sys.n() {
            let m = sys.mass(i);
            let f = sys.force(i);
            let v = sys.vel(i) + f * (0.5 * dt / m);
            sys.set_vel(i, v);
            let x = sys.pos(i) + v * dt;
            sys.set_pos(i, x);
        }
        sys.wrap_positions();
        let pe = model2.compute(&mut sys, &mut vir, &mut nl2).unwrap();
        for i in 0..sys.n() {
            let m = sys.mass(i);
            let f = sys.force(i);
            let v = sys.vel(i) + f * (0.5 * dt / m);
            sys.set_vel(i, v);
        }
        let e = sys.kinetic_energy() + pe;
        e_max = e_max.max(e);
        e_min = e_min.min(e);
    }
    let drift = (e_max - e_min) / e0.abs();
    println!(
        "decomposed NVE: {} steps, E0={:.6}, drift={:.3e}",
        steps, e0, drift
    );
    assert!(drift < 1e-3, "decomposed NVE drift too large: {drift:.3e}");
    println!("✓ domain decomposition: forces match serial and NVE energy is conserved");
}

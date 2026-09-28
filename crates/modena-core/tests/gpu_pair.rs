//! GPU pair-kernel integration test (M9, `gpu` feature).
//!
//! Skips automatically when no GPU adapter is available (e.g. headless CI):
//! [`GpuPairForceModel::new`] returns `None` and the test exits early. Where an
//! adapter exists, it compares the GPU `f32` LJ result against the serial `f64`
//! half-list reference — they must agree to single precision.

#![cfg(feature = "gpu")]

use modena_core::{
    force_model::{ForceModel, PairForceModel, Virial},
    gpu_comm::{GpuGather, GpuPairForceModel},
    lattice,
    neighbor::NeighborList,
    potential::LennardJones,
    vec3::Vec3,
    System,
};

#[test]
fn gpu_pair_matches_serial_within_f32_tolerance() {
    let Some(gpu) = GpuPairForceModel::new(1.0, 1.0, 2.5) else {
        eprintln!("no GPU adapter available — skipping GPU pair test");
        return;
    };

    // Small FCC system so the f32 all-pairs kernel is cheap and deterministic.
    let mut sys = lattice::fcc(1.5, [4, 4, 4], "Ar"); // 256 atoms, box 6.0
    sys.zero_forces();
    let mut vir_gpu = Virial::default();
    let mut nl = NeighborList::new(2.5);
    let pe_gpu = gpu.compute(&mut sys, &mut vir_gpu, &mut nl).unwrap();

    // Serial f64 reference.
    let mut s = lattice::fcc(1.5, [4, 4, 4], "Ar");
    s.zero_forces();
    let mut vir_ref = Virial::default();
    let mut nl_ref = NeighborList::new(2.5);
    let pe_ref = PairForceModel::new(LennardJones::new(1.0, 1.0, 2.5).unwrap())
        .compute(&mut s, &mut vir_ref, &mut nl_ref)
        .unwrap();

    // Single-precision agreement: relative tolerance with a floor for near-zero
    // magnitudes. The GPU computes in f32, the CPU in f64.
    let close = |a: f64, b: f64| (a - b).abs() <= 1e-3 * a.abs().max(b.abs()).max(1.0);
    assert!(close(pe_gpu, pe_ref), "GPU pe {pe_gpu} vs serial {pe_ref}");
    assert!(
        close(vir_gpu.w, vir_ref.w),
        "GPU virial {} vs serial {}",
        vir_gpu.w,
        vir_ref.w
    );
    for i in 0..sys.n() {
        for k in 0..3 {
            let g = sys.forces[3 * i + k];
            let c = s.forces[3 * i + k];
            assert!(
                close(g, c),
                "force atom {i} axis {k}: GPU {g} vs serial {c}"
            );
        }
    }
}

/// Seed deterministic small velocities and remove the center-of-mass drift.
fn seed_velocities(sys: &mut System) {
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
}

#[test]
fn gpu_nve_energy_conservation() {
    // The real MD exercise: a velocity-Verlet NVE trajectory driven by the GPU
    // force model. Total energy must be conserved — proving the GPU kernel is a
    // *conservative* force field under integration, not just correct at one
    // snapshot. The GPU computes in f32, so the drift tolerance is looser than
    // the f64 CPU path but still far below any non-conservative failure.
    let Some(gpu) = GpuPairForceModel::new(1.0, 1.0, 2.5) else {
        eprintln!("no GPU adapter available — skipping GPU NVE test");
        return;
    };

    let mut sys = lattice::fcc(1.5, [4, 4, 4], "Ar"); // 256 atoms, box 6.0
    seed_velocities(&mut sys);

    let mut vir = Virial::default();
    let mut nl = NeighborList::new(2.5);
    sys.zero_forces();
    let pe0 = gpu.compute(&mut sys, &mut vir, &mut nl).unwrap();
    let e0 = sys.kinetic_energy() + pe0;

    let dt = 0.004;
    let steps = 300;
    let mut e_max = e0;
    let mut e_min = e0;
    for _ in 0..steps {
        // Velocity-Verlet: half-kick, drift, force, half-kick.
        for i in 0..sys.n() {
            let m = sys.mass(i);
            let f = sys.force(i);
            let v = sys.vel(i) + f * (0.5 * dt / m);
            sys.set_vel(i, v);
            let x = sys.pos(i) + v * dt;
            sys.set_pos(i, x);
        }
        sys.wrap_positions();
        sys.zero_forces();
        let pe = gpu.compute(&mut sys, &mut vir, &mut nl).unwrap();
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
    let drift = (e_max - e_min) / e0.abs().max(1e-12);
    assert!(
        drift < 5e-2,
        "GPU NVE energy drift too large: {drift:.3e} (f32 kernel)"
    );
}

#[test]
fn gpu_pair_scales_to_larger_system() {
    // A second system size (864 atoms, box 9.0) — guards against a fixed-size
    // assumption in the buffer/dispatch math and confirms the f32 agreement
    // holds as N grows (the all-pairs kernel must dispatch enough workgroups).
    let Some(gpu) = GpuPairForceModel::new(1.0, 1.0, 2.5) else {
        eprintln!("no GPU adapter available — skipping GPU scaling test");
        return;
    };

    let mut sys = lattice::fcc(1.5, [6, 6, 6], "Ar"); // 864 atoms, box 9.0
    sys.zero_forces();
    let mut vir_gpu = Virial::default();
    let mut nl = NeighborList::new(2.5);
    let pe_gpu = gpu.compute(&mut sys, &mut vir_gpu, &mut nl).unwrap();

    let mut s = lattice::fcc(1.5, [6, 6, 6], "Ar");
    s.zero_forces();
    let mut vir_ref = Virial::default();
    let mut nl_ref = NeighborList::new(2.5);
    let pe_ref = PairForceModel::new(LennardJones::new(1.0, 1.0, 2.5).unwrap())
        .compute(&mut s, &mut vir_ref, &mut nl_ref)
        .unwrap();

    let close = |a: f64, b: f64| (a - b).abs() <= 1e-3 * a.abs().max(b.abs()).max(1.0);
    assert!(close(pe_gpu, pe_ref), "GPU pe {pe_gpu} vs serial {pe_ref}");
    for i in 0..sys.n() {
        for k in 0..3 {
            let g = sys.forces[3 * i + k];
            let c = s.forces[3 * i + k];
            assert!(
                close(g, c),
                "force atom {i} axis {k}: GPU {g} vs serial {c}"
            );
        }
    }
}

#[test]
fn gpu_cell_list_matches_all_pairs_and_serial() {
    // The cell-list (CSR) kernel must agree with the all-pairs oracle and the
    // serial f64 reference on a thermal (non-lattice) snapshot, through a
    // neighbor-list rebuild (skin 0.3) — exercising the CSR build/upload path.
    let (Some(gpu_ap), Some(gpu_csr)) = (
        GpuPairForceModel::new(1.0, 1.0, 2.5),
        GpuPairForceModel::with_gather(1.0, 1.0, 2.5, GpuGather::CellList),
    ) else {
        eprintln!("no GPU adapter available — skipping GPU cell-list test");
        return;
    };

    let mut sys = lattice::fcc(1.5, [4, 4, 4], "Ar"); // 256 atoms
    seed_velocities(&mut sys);
    // A few free-flight steps to decorrelate from the perfect lattice (a
    // lattice snapshot can hide range/order bugs behind symmetry).
    for i in 0..sys.n() {
        let x = sys.pos(i) + sys.vel(i) * 0.05;
        sys.set_pos(i, x);
    }
    sys.wrap_positions();

    let mut s_ap = sys.clone();
    s_ap.zero_forces();
    let mut vir_ap = Virial::default();
    let mut nl_ap = NeighborList::with_skin(2.5, 0.3);
    let pe_ap = gpu_ap.compute(&mut s_ap, &mut vir_ap, &mut nl_ap).unwrap();

    let mut s_csr = sys.clone();
    s_csr.zero_forces();
    let mut vir_csr = Virial::default();
    let mut nl_csr = NeighborList::with_skin(2.5, 0.3);
    let pe_csr = gpu_csr
        .compute(&mut s_csr, &mut vir_csr, &mut nl_csr)
        .unwrap();

    let mut s = sys.clone();
    s.zero_forces();
    let mut vir_ref = Virial::default();
    let mut nl_ref = NeighborList::new(2.5);
    let pe_ref = PairForceModel::new(LennardJones::new(1.0, 1.0, 2.5).unwrap())
        .compute(&mut s, &mut vir_ref, &mut nl_ref)
        .unwrap();

    let close = |a: f64, b: f64| (a - b).abs() <= 1e-3 * a.abs().max(b.abs()).max(1.0);
    assert!(close(pe_csr, pe_ap), "CSR pe {pe_csr} vs all-pairs {pe_ap}");
    assert!(close(pe_csr, pe_ref), "CSR pe {pe_csr} vs serial {pe_ref}");
    assert!(
        close(vir_csr.w, vir_ap.w),
        "CSR virial {} vs all-pairs {}",
        vir_csr.w,
        vir_ap.w
    );
    for i in 0..sys.n() {
        for k in 0..3 {
            let g = s_csr.forces[3 * i + k];
            let a = s_ap.forces[3 * i + k];
            let c = s.forces[3 * i + k];
            assert!(
                close(g, a),
                "CSR vs all-pairs force atom {i} axis {k}: {g} vs {a}"
            );
            assert!(
                close(g, c),
                "CSR vs serial force atom {i} axis {k}: {g} vs {c}"
            );
        }
    }
}

#[test]
fn gpu_cell_list_nve_conservation_and_rebuilds() {
    // NVE through the cell-list kernel with a Verlet skin: the list rebuilds
    // during the trajectory (displacement > skin/2), exercising the
    // fingerprint → CSR re-upload path, and the energy must stay conserved to
    // the same f32 tolerance as the all-pairs kernel.
    let Some(gpu) = GpuPairForceModel::with_gather(1.0, 1.0, 2.5, GpuGather::CellList) else {
        eprintln!("no GPU adapter available — skipping GPU cell-list NVE test");
        return;
    };

    let mut sys = lattice::fcc(1.5, [4, 4, 4], "Ar"); // 256 atoms
    seed_velocities(&mut sys);

    let mut vir = Virial::default();
    let mut nl = NeighborList::with_skin(2.5, 0.3);
    sys.zero_forces();
    let pe0 = gpu.compute(&mut sys, &mut vir, &mut nl).unwrap();
    let e0 = sys.kinetic_energy() + pe0;

    let dt = 0.004;
    let steps = 300;
    let mut e_max = e0;
    let mut e_min = e0;
    for step in 0..steps {
        for i in 0..sys.n() {
            let m = sys.mass(i);
            let f = sys.force(i);
            let v = sys.vel(i) + f * (0.5 * dt / m);
            sys.set_vel(i, v);
            let x = sys.pos(i) + v * dt;
            sys.set_pos(i, x);
        }
        sys.wrap_positions();
        // Force a deterministic rebuild mid-run so the fingerprint → CSR
        // re-upload path is exercised regardless of thermal displacements.
        if step == 150 {
            nl.force_rebuild(&sys);
        }
        sys.zero_forces();
        let pe = gpu.compute(&mut sys, &mut vir, &mut nl).unwrap();
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
    assert!(
        nl.rebuild_count() >= 2,
        "expected neighbor-list rebuilds during the NVE run, got {}",
        nl.rebuild_count()
    );
    let drift = (e_max - e_min) / e0.abs().max(1e-12);
    assert!(
        drift < 5e-2,
        "cell-list GPU NVE energy drift too large: {drift:.3e} (f32 kernel)"
    );
}

#[test]
fn gpu_cell_list_empty_and_far_pairs() {
    // Edge cases: a system where no pair is within the cutoff (empty j-ranges
    // — the j_list buffer is the 4-byte dummy) and one where a list built for
    // a smaller cutoff than the model must be rejected.
    let Some(gpu) = GpuPairForceModel::with_gather(1.0, 1.0, 2.5, GpuGather::CellList) else {
        eprintln!("no GPU adapter available — skipping GPU cell-list edge test");
        return;
    };

    // Two atoms far apart in a big box: zero pairs.
    let mut sys = System::new(modena_core::SimBox::new(Vec3::new(20.0, 20.0, 20.0)));
    sys.resize(2);
    sys.set_pos(0, Vec3::new(1.0, 1.0, 1.0));
    sys.set_pos(1, Vec3::new(10.0, 10.0, 10.0));
    sys.zero_forces();
    let mut vir = Virial::default();
    let mut nl = NeighborList::with_skin(2.5, 0.3);
    let pe = gpu.compute(&mut sys, &mut vir, &mut nl).unwrap();
    assert!(pe.abs() < 1e-12, "far pair must give zero pe, got {pe}");
    assert!(vir.w.abs() < 1e-12);
    for f in &sys.forces {
        assert!(f.abs() < 1e-12);
    }

    // A list whose cutoff does not cover the model cutoff must be an error,
    // not a silently wrong force field.
    let mut nl_short = NeighborList::with_skin(2.0, 0.1);
    let err = gpu.compute(&mut sys, &mut vir, &mut nl_short).unwrap_err();
    assert!(
        err.to_string().contains("neighbor list"),
        "expected cutoff-mismatch error, got {err}"
    );
}

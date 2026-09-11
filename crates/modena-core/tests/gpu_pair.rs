//! GPU pair-kernel integration test (M9, `gpu` feature).
//!
//! Skips automatically when no GPU adapter is available (e.g. headless CI):
//! [`GpuPairForceModel::new`] returns `None` and the test exits early. Where an
//! adapter exists, it compares the GPU `f32` LJ result against the serial `f64`
//! half-list reference — they must agree to single precision.

#![cfg(feature = "gpu")]

use modena_core::{
    force_model::{ForceModel, PairForceModel, Virial},
    gpu_comm::GpuPairForceModel,
    lattice,
    neighbor::NeighborList,
    potential::LennardJones,
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

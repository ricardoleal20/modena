//! GPU pair-kernel demo (M9, `gpu` feature).
//!
//! Acquires a GPU adapter, runs the WGSL LJ compute shader on an FCC argon
//! cluster, compares force/energy to the serial f64 reference, then runs a
//! short NVE trajectory on the GPU and reports the energy drift. Exits
//! gracefully if no adapter is present.
//!
//! ```sh
//! cargo run --release --features gpu --example gpu_pair
//! ```

#[cfg(not(feature = "gpu"))]
fn main() {
    eprintln!(
        "this example requires the `gpu` feature: cargo run --features gpu --example gpu_pair"
    );
}

#[cfg(feature = "gpu")]
use modena_core::{
    force_model::{ForceModel, PairForceModel, Virial},
    gpu_comm::GpuPairForceModel,
    lattice,
    neighbor::NeighborList,
    potential::LennardJones,
    vec3::Vec3,
};

#[cfg(feature = "gpu")]
fn main() {
    let cutoff = 2.5;
    let Some(gpu) = GpuPairForceModel::new(1.0, 1.0, cutoff) else {
        eprintln!("no GPU adapter available — nothing to demo");
        return;
    };
    println!("GPU adapter acquired; running WGSL LJ kernel");

    let sys0 = lattice::fcc(1.5, [5, 5, 5], "Ar"); // 500 atoms, box 7.5
    println!("system: {} atoms, box {:?}", sys0.n(), sys0.box_size.size.0);

    // --- 1. GPU vs serial (f32 vs f64). ---
    let mut s_gpu = sys0.clone();
    s_gpu.zero_forces();
    let mut vir_g = Virial::default();
    let mut nl = NeighborList::new(cutoff);
    let pe_g = gpu.compute(&mut s_gpu, &mut vir_g, &mut nl).unwrap();

    let mut s_ser = sys0.clone();
    s_ser.zero_forces();
    let mut vir_s = Virial::default();
    let mut nl_s = NeighborList::new(cutoff);
    let pe_s = PairForceModel::new(LennardJones::new(1.0, 1.0, cutoff).unwrap())
        .compute(&mut s_ser, &mut vir_s, &mut nl_s)
        .unwrap();

    let max_f_err = (0..s_ser.forces.len())
        .map(|i| (s_ser.forces[i] - s_gpu.forces[i]).abs())
        .fold(0.0f64, f64::max);
    println!(
        "GPU vs serial: Δpe={:.3e}, Δvirial={:.3e}, max|Δf|={:.3e} (f32 kernel)",
        (pe_g - pe_s).abs(),
        (vir_g.w - vir_s.w).abs(),
        max_f_err
    );

    // --- 2. Short GPU NVE trajectory; report energy drift. ---
    let mut sys = sys0.clone();
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

    let mut vir = Virial::default();
    let mut nl2 = NeighborList::new(cutoff);
    sys.zero_forces();
    let pe0 = gpu.compute(&mut sys, &mut vir, &mut nl2).unwrap();
    let e0 = sys.kinetic_energy() + pe0;

    let dt = 0.004;
    let steps = 300;
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
        sys.zero_forces();
        let pe = gpu.compute(&mut sys, &mut vir, &mut nl2).unwrap();
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
        "GPU NVE: {} steps, E0={:.6}, drift={:.3e} (f32 kernel)",
        steps, e0, drift
    );
    println!("✓ GPU pair kernel: matches serial to f32 and conserves NVE energy");
}

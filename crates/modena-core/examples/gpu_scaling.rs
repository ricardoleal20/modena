//! GPU vs CPU pair-kernel scaling benchmark (M9, `gpu` feature).
//!
//! Times `GpuPairForceModel::compute` (f32 WGSL shader, O(N^2) all-pairs) vs
//! `PairForceModel::compute` (f64 rayon half-list, O(N) with cutoff) over a
//! fixed step count at several system sizes. Reports wall time and speedup at
//! each size. The GPU kernel computes in f32 while the CPU path uses f64 —
//! the speedup comparison is the point.
//!
//! ```sh
//! cargo run --release --features gpu --example gpu_scaling
//! ```
//!
//! Prints CSV: `n_atoms,cpu_wall_s,cpu_us_per_step,gpu_wall_s,gpu_us_per_step,speedup`

#[cfg(not(feature = "gpu"))]
fn main() {
    eprintln!(
        "this example requires the `gpu` feature: cargo run --features gpu --example gpu_scaling"
    );
}

#[cfg(feature = "gpu")]
use modena_core::{
    force_model::{ForceModel, PairForceModel, Virial},
    gpu_comm::GpuPairForceModel,
    lattice,
    neighbor::NeighborList,
    potential::LennardJones,
};

/// Choose roughly cubic FCC replicas to reach ~`target` atoms (4 per cell).
#[cfg(feature = "gpu")]
fn fcc_replicas_for_target(target: usize) -> [u32; 3] {
    let cells = (target as f64 / 4.0).cbrt().round() as u32;
    let cells = cells.max(1);
    [cells, cells, cells]
}

/// Time `steps` force evaluations of `model` on `sys`, returning (wall_s, us_per_step).
/// Forces and virial are zeroed before each step (the trait contract is accumulate).
#[cfg(feature = "gpu")]
fn time_model(
    model: &dyn ForceModel,
    sys: &mut modena_core::system::System,
    steps: usize,
    cutoff: f64,
) -> (f64, f64) {
    let mut virial = Virial::default();
    let mut nl = NeighborList::new(cutoff);

    // Warmup (untimed).
    sys.zero_forces();
    let _ = model.compute(sys, &mut virial, &mut nl).unwrap();

    let t0 = std::time::Instant::now();
    for _ in 0..steps {
        sys.zero_forces();
        virial = Virial::default();
        let _ = model.compute(sys, &mut virial, &mut nl).unwrap();
    }
    let wall = t0.elapsed().as_secs_f64();
    (wall, wall / steps as f64 * 1.0e6)
}

#[cfg(feature = "gpu")]
fn main() {
    let cutoff = 2.5;
    let a = 1.6796; // rho ~ 0.8442
    let steps = 100;

    let Some(gpu) = GpuPairForceModel::new(1.0, 1.0, cutoff) else {
        eprintln!("no GPU adapter available — nothing to benchmark");
        return;
    };
    eprintln!("GPU adapter acquired; running WGSL LJ kernel vs CPU rayon");

    let cpu_model = PairForceModel::new(LennardJones::new(1.0, 1.0, cutoff).unwrap());

    // ~10^3, 10^4, 5x10^4 atoms (drop largest if too slow).
    let sizes: &[usize] = &[1_000, 10_000, 50_000];

    println!("n_atoms,cpu_wall_s,cpu_us_per_step,gpu_wall_s,gpu_us_per_step,speedup");

    for &target in sizes {
        let replicas = fcc_replicas_for_target(target);
        let sys0 = lattice::fcc(a, replicas, "Ar");
        let n = sys0.n();
        eprintln!("  testing {} atoms (target {})...", n, target);

        let mut sys_cpu = sys0.clone();
        let (cpu_wall, cpu_us) = time_model(&cpu_model, &mut sys_cpu, steps, cutoff);

        let mut sys_gpu = sys0.clone();
        let (gpu_wall, gpu_us) = time_model(&gpu, &mut sys_gpu, steps, cutoff);

        let speedup = cpu_wall / gpu_wall;
        println!(
            "{},{:.6},{:.2},{:.6},{:.2},{:.2}",
            n, cpu_wall, cpu_us, gpu_wall, gpu_us, speedup
        );
    }
}

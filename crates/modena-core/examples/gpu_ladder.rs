//! GPU-d decision-ladder driver (committee G3, `gpu` feature).
//!
//! Runs full velocity-Verlet NVE integration steps — not micro-benchmarks —
//! for one of four engines at one ladder size, and prints per-run plus
//! summary CSV for the ladder table in `docs/PERF_RESULTS.md`:
//!
//! - `cpu` — the shipped `PairForceModel` (f64, rayon, default threads),
//!   driven through the *same* `VelocityVerlet` integrator as the runner, so
//!   CPU vs GPU differ only in the force model.
//! - `gpu` — the WGSL all-pairs LJ kernel (`GpuGather::AllPairs`). Every
//!   `compute` blocks on the staging-buffer map
//!   (`device.poll(PollType::wait_indefinitely())`), so the per-step
//!   `Instant` window is fence-synced: the GPU work is complete before the
//!   clock reads.
//! - `gpu-csr` — the cell-list gather kernel (`GpuGather::CellList`): CPU
//!   Verlet bins → contiguous per-atom j-ranges (CSR) → one GPU thread per
//!   i, no atomics. Also fence-synced.
//! - `gather` — the all-pairs kernel with a 0.01 cutoff: the j-loop still
//!   streams all N partner positions (gather + minimum image + distance) but
//!   no pair passes the cutoff, so the LJ math is skipped. Its time bounds
//!   what a cell-list gather could save (the cell-list replaces the N²
//!   gathers with ~55N).
//!
//! Initial state matches the matched decks bit-for-bit (same FCC lattice,
//! same `init_velocities` seed 87287, dt 0.005, reduced units), so the CPU
//! and GPU drift rows are the same trajectory under the two engines.
//!
//! ```sh
//! cargo run --release --features gpu --example gpu_ladder -- gpu 2048 2000 5
//! ```
//!
//! Usage: `gpu_ladder <cpu|gpu|gpu-csr|gather> <atoms> <steps> <runs>` (one
//! warm-up run is always executed and discarded). Supported atom counts:
//! 864, 2048, 5000, 10000, 20000 (rectangular FCC blocks; see
//! `replicas_for`).

#[cfg(not(feature = "gpu"))]
fn main() {
    eprintln!(
        "this example requires the `gpu` feature: cargo run --features gpu --example gpu_ladder"
    );
}

#[cfg(feature = "gpu")]
use modena_core::{
    conservation::{ConservationReport, ConservationVerdict},
    force_model::{ForceModel, PairForceModel, Virial},
    gpu_comm::{GpuGather, GpuPairForceModel},
    integrator::{Integrator, VelocityVerlet},
    lattice,
    neighbor::NeighborList,
    potential::LennardJones,
    runner::init_velocities,
};

/// Ladder sizes → FCC replica counts (4 atoms/cell). 5000/10000/20000 are not
/// 4·n³, so they use rectangular blocks with the smallest aspect ratio that
/// factorizes exactly; the smallest box edge (5·1.6796 = 8.398) still exceeds
/// 2·cutoff = 5.0, so minimum image is well defined.
#[cfg(feature = "gpu")]
fn replicas_for(atoms: usize) -> [u32; 3] {
    match atoms {
        864 => [6, 6, 6],
        2048 => [8, 8, 8],
        5000 => [5, 10, 25],
        10000 => [10, 10, 25],
        20000 => [10, 20, 25],
        _ => {
            eprintln!("unsupported ladder size {atoms} (supported: 864 2048 5000 10000 20000)");
            std::process::exit(2);
        }
    }
}

/// One NVE run: fresh system at the deck seed, `steps` velocity-Verlet steps
/// through `model`, per-step timings, and a conservation report. Returns
/// (wall_seconds, per_step_median_ms, conservation, energy band / |E0|).
#[cfg(feature = "gpu")]
fn run_once(
    model: &dyn ForceModel,
    replicas: [u32; 3],
    steps: usize,
) -> (f64, f64, ConservationReport, f64) {
    const A: f64 = 1.6796; // rho = 0.8442
    const CUTOFF: f64 = 2.5;
    const SKIN: f64 = 0.3; // matches the decks' run.neighbor_skin default
    const DT: f64 = 0.005;
    const SEED: u64 = 87287;

    let mut sys = lattice::fcc(A, replicas, "Ar");
    init_velocities(&mut sys, 1.0, SEED);

    let mut neighbors = NeighborList::with_skin_exclusions(CUTOFF, SKIN, &[]);
    let mut vv = VelocityVerlet::new();

    // Initial force evaluation so F(t) is cached for the first half-kick
    // (mirrors `Simulation::from_config`). Untimed.
    sys.forces.fill(0.0);
    let mut virial = Virial::default();
    let pe0 = model
        .compute(&mut sys, &mut virial, &mut neighbors)
        .unwrap();
    let e0 = sys.kinetic_energy() + pe0;
    let mut report = ConservationReport::pending(e0);

    let mut step_ms = Vec::with_capacity(steps);
    let mut e_min = e0;
    let mut e_max = e0;
    let t0 = std::time::Instant::now();
    for _ in 0..steps {
        // Fence-synced on the GPU path: `compute` returns only after the
        // staging-buffer map completes, i.e. after the GPU work is done.
        let s0 = std::time::Instant::now();
        let e = vv.step(&mut sys, model, &mut neighbors, DT).unwrap();
        step_ms.push(s0.elapsed().as_secs_f64() * 1.0e3);
        report.observe(e.total);
        e_min = e_min.min(e.total);
        e_max = e_max.max(e.total);
    }
    let wall = t0.elapsed().as_secs_f64();
    step_ms.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let step_median = step_ms[steps / 2];
    (wall, step_median, report, (e_max - e_min) / e0.abs())
}

#[cfg(feature = "gpu")]
fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() != 5 {
        eprintln!("usage: gpu_ladder <cpu|gpu|gather> <atoms> <steps> <runs>");
        std::process::exit(2);
    }
    let engine = args[1].clone();
    let atoms: usize = args[2].parse().unwrap();
    let steps: usize = args[3].parse().unwrap();
    let runs: usize = args[4].parse().unwrap_or(1).max(1);

    let replicas = replicas_for(atoms);
    let cpu_model = || -> Box<dyn ForceModel> {
        Box::new(PairForceModel::new(
            LennardJones::new(1.0, 1.0, 2.5).unwrap(),
        ))
    };
    let gpu_model = |cutoff, gather| -> Box<dyn ForceModel> {
        match GpuPairForceModel::with_gather(1.0, 1.0, cutoff, gather) {
            Some(m) => Box::new(m),
            None => {
                eprintln!("no GPU adapter available — nothing to ladder");
                std::process::exit(3);
            }
        }
    };
    let model: Box<dyn ForceModel> = match engine.as_str() {
        "cpu" => {
            eprintln!("engine: MODENA CPU (f64, rayon default threads), {atoms} atoms");
            cpu_model()
        }
        "gpu" => {
            eprintln!("engine: MODENA GPU all-pairs (WGSL LJ, f32+Kahan), {atoms} atoms");
            gpu_model(2.5, GpuGather::AllPairs)
        }
        "gpu-csr" => {
            eprintln!(
                "engine: MODENA GPU cell-list gather (CSR j-ranges, f32+Kahan), {atoms} atoms"
            );
            gpu_model(2.5, GpuGather::CellList)
        }
        "gather" => {
            // Same j-loop (gather + minimum image + distance), LJ math skipped
            // for every pair: bounds the cell-list-gather saving.
            eprintln!("engine: GPU gather probe (cutoff 0.01, no LJ math), {atoms} atoms");
            gpu_model(0.01, GpuGather::AllPairs)
        }
        _ => {
            eprintln!("unknown engine {engine} (cpu | gpu | gpu-csr | gather)");
            std::process::exit(2);
        }
    };

    println!("run,engine,atoms,steps,wall_s,mean_ms_per_step,median_ms_per_step,drift_final,band,verdict");
    let mut means = Vec::with_capacity(runs);
    let mut meds = Vec::with_capacity(runs);
    let mut summary_drift = None;
    for run in 0..=runs {
        // Run 0 is the warm-up (buffers, shader JIT, caches); discarded.
        let (wall, step_median, report, band) = run_once(model.as_ref(), replicas, steps);
        if run == 0 {
            eprintln!("warm-up run discarded ({wall:.3} s)");
            continue;
        }
        means.push(wall * 1.0e3 / steps as f64);
        meds.push(step_median);
        let (drift_s, verdict): (String, String) = if engine == "gather" {
            // Free flight (no pair passes the cutoff): drift is meaningless.
            ("n/a".into(), "n/a".into())
        } else {
            summary_drift = Some((report.relative_drift, band));
            (
                format!("{:.3e}", report.relative_drift),
                format!("{:?}", report.verdict),
            )
        };
        println!(
            "{},{},{},{},{:.4},{:.4},{:.4},{},{:.3e},{}",
            run,
            engine,
            atoms,
            steps,
            wall,
            means[run - 1],
            step_median,
            drift_s,
            band,
            verdict,
        );
    }

    means.sort_by(|a, b| a.partial_cmp(b).unwrap());
    meds.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let median = means[means.len() / 2];
    let min = means[0];
    let max = means[means.len() - 1];
    let med_of_meds = meds[meds.len() / 2];
    let (drift_s, band_s) = match summary_drift {
        Some((d, b)) => (format!("{d:.3e}"), format!("{b:.3e}")),
        None => ("n/a".into(), "n/a".into()),
    };
    let verdict_s = match summary_drift {
        Some((d, _)) => format!("{:?}", ConservationVerdict::from_drift(d)),
        None => "n/a".into(),
    };
    println!(
        "summary,{engine},{atoms},{runs},{min_ms:.4},{median_ms:.4},{max_ms:.4},{med_step_ms:.4},{drift},{band},{verdict}",
        min_ms = min,
        median_ms = median,
        max_ms = max,
        med_step_ms = med_of_meds,
        drift = drift_s,
        band = band_s,
        verdict = verdict_s,
    );
}

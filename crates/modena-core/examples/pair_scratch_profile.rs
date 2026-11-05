//! Isolated allocation/static-numerics profile. Timings here are instrumented,
//! never the primary performance comparison (use bench/pair_scratch.py).
#[path = "../tests/support/pair_scratch_alloc.rs"]
mod allocation;
use allocation::*;
use modena_core::{ForceModel, LennardJones, NeighborList, PairForceModel, Virial};
use std::sync::atomic::Ordering;

#[global_allocator]
static ALLOCATOR: CountingAllocator = CountingAllocator;

fn main() {
    let args: Vec<_> = std::env::args().collect();
    assert_eq!(
        args.len(),
        4,
        "usage: pair_scratch_profile ATOMS THREADS PASSES"
    );
    let n: usize = args[1].parse().unwrap();
    let threads: usize = args[2].parse().unwrap();
    let passes: usize = args[3].parse().unwrap();
    assert!(threads > 0 && passes > 0);
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .unwrap();
    let model = PairForceModel::new(LennardJones::new(1.0, 1.0, 2.5).unwrap());
    let mut s = system(n);
    assert_eq!(s.n(), n);
    // Break FCC symmetry for a meaningful per-component comparison.
    for (i, x) in s.positions.iter_mut().enumerate() {
        *x += 0.013 * (i as f64 * 0.7).sin();
    }
    let mut neighbors = NeighborList::with_skin(2.5, 0.3);
    let mut v = Virial::default();
    pool.install(|| model.compute(&mut s, &mut v, &mut neighbors).unwrap());
    let pairs = neighbors.n_pairs();
    reset(3 * n * size_of::<f64>());
    let live_before = LIVE_BYTES.load(Ordering::SeqCst);
    let start = std::time::Instant::now();
    let mut u = 0.0;
    for _ in 0..passes {
        s.zero_forces();
        v.w = 0.0;
        u = pool.install(|| model.compute(&mut s, &mut v, &mut neighbors).unwrap());
    }
    let elapsed = start.elapsed().as_secs_f64();
    let allocations = FORCE_ALLOCS.load(Ordering::SeqCst);
    let reallocations = FORCE_REALLOCS.load(Ordering::SeqCst);
    let bytes = ALLOCATED_BYTES.load(Ordering::SeqCst);
    let live_after = LIVE_BYTES.load(Ordering::SeqCst);
    let high_water = HIGH_WATER.load(Ordering::SeqCst);
    FORCE_BYTES.store(0, Ordering::SeqCst);
    let mut reference = s.clone();
    reference.zero_forces();
    let mut rv = Virial::default();
    let ru = model
        .compute_deterministic(&mut reference, &mut rv, &mut neighbors)
        .unwrap();
    let close = |a: f64, b: f64, atol: f64| (a - b).abs() <= atol + 1e-12 * b.abs();
    let exact = u == ru && v.w == rv.w && s.forces == reference.forces;
    let passed = close(u, ru, 1e-8)
        && close(v.w, rv.w, 1e-8)
        && s.forces
            .iter()
            .zip(&reference.forces)
            .all(|(&a, &b)| close(a, b, 1e-10))
        && (threads != 1 || exact);
    let max_force_error = s
        .forces
        .iter()
        .zip(&reference.forces)
        .map(|(a, b)| (a - b).abs())
        .fold(0.0, f64::max);
    let before_drop = LIVE_BYTES.load(Ordering::SeqCst);
    drop(model);
    let released = before_drop - LIVE_BYTES.load(Ordering::SeqCst);
    println!(
        "{}",
        serde_json::json!({"atoms":n,"threads":threads,"passes":passes,"pairs":pairs,
        "parallel_eligible":pairs>=100_000 && threads>1,"force_bytes":3*n*8,
        "warm_force_allocations":allocations,"warm_force_reallocations":reallocations,
        "warm_force_allocated_bytes":bytes,"live_bytes_before":live_before,"live_bytes_after":live_after,
        "live_bytes_high_water":high_water,"model_drop_released_bytes":released,
        "instrumented_seconds":elapsed,"energy":u,"virial":v.w,"serial_energy":ru,"serial_virial":rv.w,
        "max_force_error":max_force_error,"one_thread_exact":exact,"static_gate_passed":passed})
    );
    assert!(passed, "fixed static numerical gate failed");
}

#[path = "support/pair_scratch_alloc.rs"]
mod allocation;
use allocation::*;
use modena_core::{ForceModel, LennardJones, NeighborList, PairForceModel, Virial};
use std::sync::atomic::Ordering;

#[global_allocator]
static ALLOCATOR: CountingAllocator = CountingAllocator;

#[test]
fn warm_parallel_pass_reuses_force_buffers() {
    // Keep allocator-sensitive checks in one isolated integration-test process.
    for threads in [1, 4, 8] {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .unwrap();
        let model = PairForceModel::new(LennardJones::new(1.0, 1.0, 2.5).unwrap());
        let mut s = system(5000);
        let mut neighbors = NeighborList::with_skin(2.5, 0.3);
        let mut virial = Virial::default();
        pool.install(|| model.compute(&mut s, &mut virial, &mut neighbors).unwrap());
        assert!(
            neighbors.n_pairs() >= 100_000,
            "must exercise parallel branch"
        );
        s.zero_forces();
        virial.w = 0.0;
        reset(3 * s.n() * size_of::<f64>());
        pool.install(|| model.compute(&mut s, &mut virial, &mut neighbors).unwrap());
        let observed = FORCE_ALLOCS.load(Ordering::SeqCst) + FORCE_REALLOCS.load(Ordering::SeqCst);
        FORCE_BYTES.store(0, Ordering::SeqCst);
        eprintln!(
            "threads={threads} atoms={} pairs={} force_allocations={observed}",
            s.n(),
            neighbors.n_pairs()
        );
        assert_eq!(
            observed, 0,
            "warm pair pass allocated an N-sized force buffer"
        );
    }
}

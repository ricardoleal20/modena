//! Cell-list neighbor-build benchmark.
//!
//! Measures `NeighborList::build` for FCC LJ systems at three sizes. The build
//! is O(N) in atoms; this isolates it from the force kernel.

use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion};
use modena_core::{lattice, NeighborList};

fn build_system(n_per_side: u32) -> modena_core::System {
    // a = 1.5444 → ρ ≈ 1.085 (dense LJ solid). Larger replicas → larger N.
    lattice::fcc(1.5444, [n_per_side, n_per_side, n_per_side], "Ar")
}

fn bench_build(c: &mut Criterion) {
    let mut group = c.benchmark_group("neighbor_build");
    for n in [4u32, 6, 8] {
        let system = build_system(n);
        group.throughput(criterion::Throughput::Elements(system.n() as u64));
        group.bench_with_input(BenchmarkId::from_parameter(n), &system, |b, sys| {
            b.iter(|| {
                let mut nl = NeighborList::new(2.5);
                nl.build(black_box(sys));
                black_box(nl);
            });
        });
    }
    group.finish();
}

criterion_group!(benches, bench_build);
criterion_main!(benches);

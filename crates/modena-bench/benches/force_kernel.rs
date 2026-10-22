//! Force clearing + model compute + neighbor displacement check. Rebuild isolated.
mod common;
use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion};
use modena_core::{NeighborList, Virial};
fn bench_force(c: &mut Criterion) {
    for (name, mut system, model, cutoff) in common::workloads() {
        let mut nl = NeighborList::with_skin(cutoff, 0.3);
        nl.build(&system);
        for rebuild in [false, true] {
            let mut group = c.benchmark_group(if rebuild {
                "force_with_rebuild"
            } else {
                "force_cached_neighbors"
            });
            group.throughput(criterion::Throughput::Elements(system.n() as u64));
            group.bench_function(BenchmarkId::new(&name, system.n()), |b| {
                b.iter(|| {
                    system.zero_forces();
                    if rebuild {
                        nl.force_rebuild(&system);
                    }
                    let u = model
                        .compute(&mut system, &mut Virial::default(), &mut nl)
                        .unwrap();
                    black_box(u);
                })
            });
            group.finish();
        }
    }
}
criterion_group!(benches, bench_force);
criterion_main!(benches);

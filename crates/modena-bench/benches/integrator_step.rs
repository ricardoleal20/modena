//! One steady velocity-Verlet step; setup clones system+cached neighbors untimed.
//! Includes kicks, force clearing and neighbor displacement checks; no forced rebuild.
mod common;
use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use modena_core::{Integrator, NeighborList, Vec3, VelocityVerlet, Virial};
use std::hint::black_box;
fn bench_step(c: &mut Criterion) {
    let mut group = c.benchmark_group("integrator_cached_neighbors");
    for (name, mut system, model, cutoff) in common::workloads() {
        for i in 0..system.n() {
            system.set_vel(i, Vec3::new(0.01 * (i as f64).sin(), 0.0, 0.0));
        }
        let mut nl = NeighborList::with_skin(cutoff, 0.3);
        model
            .compute(&mut system, &mut Virial::default(), &mut nl)
            .unwrap();
        group.throughput(criterion::Throughput::Elements(system.n() as u64));
        group.bench_function(BenchmarkId::new(name, system.n()), |b| {
            b.iter_batched(
                || {
                    let s = system.clone();
                    let mut local_nl = NeighborList::with_skin(cutoff, 0.3);
                    local_nl.build(&s);
                    (s, local_nl)
                },
                |(mut s, mut local_nl)| {
                    let e = VelocityVerlet::new()
                        .step(&mut s, &*model, &mut local_nl, 0.0001)
                        .unwrap();
                    black_box(e);
                },
                criterion::BatchSize::SmallInput,
            )
        });
    }
    group.finish();
}
criterion_group!(benches, bench_step);
criterion_main!(benches);

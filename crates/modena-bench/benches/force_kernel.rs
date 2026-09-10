//! Force-kernel benchmark: one `ForceModel::compute` call (neighbor build +
//! rayon-parallel pair evaluation) on FCC LJ systems.

use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion};
use modena_core::{
    lattice, ForceModel, LennardJones, NeighborList, PairForceModel, PairPotential, System, Virial,
};

fn build_system(n_per_side: u32) -> System {
    lattice::fcc(1.5444, [n_per_side, n_per_side, n_per_side], "Ar")
}

fn bench_force(c: &mut Criterion) {
    let lj = LennardJones::new(1.0, 1.0, 2.5).unwrap();
    let cutoff = lj.cutoff();
    let model = PairForceModel::new(lj);

    let mut group = c.benchmark_group("force_kernel");
    for n in [4u32, 6, 8] {
        let system = build_system(n);
        group.throughput(criterion::Throughput::Elements(system.n() as u64));
        group.bench_with_input(BenchmarkId::from_parameter(n), &system, |b, sys| {
            b.iter_batched(
                || sys.clone(),
                |mut s| {
                    let mut nl = NeighborList::new(cutoff);
                    let mut virial = Virial::default();
                    let pe = model.compute(&mut s, &mut virial, &mut nl).unwrap();
                    black_box(pe);
                },
                criterion::BatchSize::SmallInput,
            );
        });
    }
    group.finish();
}

criterion_group!(benches, bench_force);
criterion_main!(benches);

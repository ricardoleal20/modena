//! Integrator-step benchmark: one velocity-Verlet step (two half-kicks + drift
//! + one force evaluation). This is the per-step cost that determines steps/s.

use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion};
use modena_core::{
    lattice, ForceModel, Integrator, LennardJones, NeighborList, PairForceModel, PairPotential,
    System, Vec3, VelocityVerlet, Virial,
};

fn build_system(n_per_side: u32) -> System {
    lattice::fcc(1.5444, [n_per_side, n_per_side, n_per_side], "Ar")
}

fn bench_step(c: &mut Criterion) {
    let lj = LennardJones::new(1.0, 1.0, 2.5).unwrap();
    let cutoff = lj.cutoff();
    let model = PairForceModel::new(lj);
    let mut integrator = VelocityVerlet::new();
    let dt = 0.005;

    let mut group = c.benchmark_group("integrator_step");
    for n in [4u32, 6, 8] {
        let mut system = build_system(n);
        // Seed a non-trivial velocity so forces are exercised realistically.
        for i in 0..system.n() {
            let f = i as f64;
            system.velocities[i] =
                Vec3::new(0.3 * f.sin(), 0.3 * (1.7 * f).cos(), 0.3 * (0.9 * f).sin());
        }
        let mut nl = NeighborList::new(cutoff);
        let mut virial = Virial::default();
        let _ = model.compute(&mut system, &mut virial, &mut nl); // cache F(t)

        group.throughput(criterion::Throughput::Elements(system.n() as u64));
        group.bench_with_input(BenchmarkId::from_parameter(n), &system, |b, sys| {
            b.iter_batched(
                || sys.clone(),
                |mut s| {
                    let mut local_nl = NeighborList::new(cutoff);
                    let e = integrator.step(&mut s, &model, &mut local_nl, dt).unwrap();
                    black_box(e);
                },
                criterion::BatchSize::SmallInput,
            );
        });
    }
    group.finish();
}

criterion_group!(benches, bench_step);
criterion_main!(benches);

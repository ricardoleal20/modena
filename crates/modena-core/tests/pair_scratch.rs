#[path = "support/pair_scratch_alloc.rs"]
mod allocation;
use allocation::*;
use modena_core::{
    Bond, BondedForceModel, BondedTopology, CompositeForceModel, ForceModel, LennardJones,
    NeighborList, PairForceModel, PairPotential, PairScale, System, Virial,
};
use std::sync::atomic::Ordering;

#[global_allocator]
static ALLOCATOR: CountingAllocator = CountingAllocator;

#[test]
fn warm_parallel_pass_reuses_force_buffers() {
    for changed in [
        (-0.0, 0.0, &[0.0][..]),
        (0.0, -0.0, &[0.0][..]),
        (0.0, 0.0, &[-0.0][..]),
    ] {
        assert!(
            !exact_state((0.0, 0.0, &[0.0]), changed),
            "exact E/F/W gate must reject opposite signed zeros"
        );
    }
    // Keep allocator-sensitive checks in one isolated integration-test process.
    let model = PairForceModel::new(LennardJones::new(1.0, 1.0, 2.5).unwrap());
    for (threads, atoms) in [
        (1, 5000),
        (4, 5000),
        (8, 10000),
        (11, 5000),
        (4, 5000),
        (4, 10000),
    ] {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .unwrap();
        let mut s = system(atoms);
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
    numerical_resize_and_composition();
    concurrent_and_reentrant_calls();
}

fn lj() -> LennardJones {
    LennardJones::new(1.0, 1.0, 2.5).unwrap()
}
fn perturbed(n: usize) -> System {
    let mut s = system(n);
    for (i, x) in s.positions.iter_mut().enumerate() {
        *x += 0.013 * (i as f64 * 0.7).sin();
    }
    s
}
fn evaluate(model: &dyn ForceModel, source: &System, scales: &[PairScale]) -> (f64, f64, Vec<f64>) {
    let mut s = source.clone();
    s.forces.fill(2.5);
    let mut v = Virial { w: 7.0 };
    let mut n = NeighborList::with_pair_scales(2.5, 0.3, scales).unwrap();
    let u = model.compute(&mut s, &mut v, &mut n).unwrap();
    (u, v.w, s.forces)
}
fn close(actual: &(f64, f64, Vec<f64>), expected: &(f64, f64, Vec<f64>), exact: bool) {
    if exact {
        assert!(exact_state(
            (actual.0, actual.1, &actual.2),
            (expected.0, expected.1, &expected.2)
        ));
        return;
    }
    for (a, b) in [(actual.0, expected.0), (actual.1, expected.1)] {
        assert!((a - b).abs() <= 1e-8 + 1e-12 * b.abs(), "scalar {a} vs {b}");
    }
    for (i, (&a, &b)) in actual.2.iter().zip(&expected.2).enumerate() {
        assert!(
            (a - b).abs() <= 1e-10 + 1e-12 * b.abs(),
            "force[{i}] {a} vs {b}"
        );
    }
}
struct Generic(LennardJones);
impl PairPotential for Generic {
    fn cutoff(&self) -> f64 {
        self.0.cutoff()
    }
    fn cutoff2(&self) -> f64 {
        self.0.cutoff2()
    }
    fn energy_force(&self, r2: f64, i: u32, j: u32) -> (f64, f64) {
        self.0.energy_force(r2, i, j)
    }
}
fn numerical_resize_and_composition() {
    let typed = PairForceModel::new(lj());
    let generic = PairForceModel::from_boxed(Box::new(Generic(lj())));
    let scales = [
        PairScale {
            i: 0,
            j: 1,
            lj: 0.25,
            coulomb: 0.0,
        },
        PairScale {
            i: 0,
            j: 2,
            lj: 0.0,
            coulomb: 0.0,
        },
    ];
    for (threads, n) in [
        (1, 5000),
        (4, 5000),
        (8, 10000),
        (3, 5000),
        (11, 10000),
        (2, 2048),
        (4, 5000),
        (8, 864),
        (4, 10000),
        (1, 5000),
    ] {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .unwrap();
        let s = perturbed(n);
        let mut reference = s.clone();
        reference.forces.fill(2.5);
        let mut v = Virial { w: 7.0 };
        let mut list = NeighborList::with_pair_scales(2.5, 0.3, &scales).unwrap();
        let u = typed
            .compute_deterministic(&mut reference, &mut v, &mut list)
            .unwrap();
        let expected = (u, v.w, reference.forces);
        for model in [&typed, &generic] {
            for _ in 0..3 {
                close(
                    &pool.install(|| evaluate(model, &s, &scales)),
                    &expected,
                    threads == 1,
                );
            }
        }
    }
    let composite = CompositeForceModel::new(typed);
    let s = perturbed(5000);
    let expected = evaluate(&generic, &s, &[]);
    close(&evaluate(&composite, &s, &[]), &expected, false);
    let mut empty = system(864);
    empty.resize(0);
    assert_eq!(evaluate(&composite, &empty, &[]), (0.0, 7.0, vec![]));
    close(&evaluate(&composite, &s, &[]), &expected, false);

    // Actual two-provider composition, including preexisting F/W. Build the
    // reference by independently summing contributions once each.
    let topology = BondedTopology {
        bonds: vec![Bond {
            i: 0,
            j: 1,
            r0: 1.0,
            k: 3.0,
        }],
        ..BondedTopology::default()
    };
    let bonded = evaluate(&BondedForceModel::new(topology.clone()), &s, &[]);
    let sum = (
        expected.0 + bonded.0,
        expected.1 + bonded.1 - 7.0,
        expected
            .2
            .iter()
            .zip(&bonded.2)
            .map(|(a, b)| a + b - 2.5)
            .collect(),
    );
    let combined = composite.with_bonded(BondedForceModel::new(topology));
    close(&evaluate(&combined, &s, &[]), &sum, false);
}

fn concurrent_and_reentrant_calls() {
    use std::sync::atomic::AtomicBool;
    use std::sync::{Arc, OnceLock, Weak};
    let shared = Arc::new(PairForceModel::new(lj()));
    let s = perturbed(5000);
    let expected = evaluate(&*shared, &s, &[]);
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(4)
        .build()
        .unwrap();
    pool.install(|| {
        rayon::scope(|scope| {
            for _ in 0..8 {
                let m = &shared;
                let s = &s;
                let expected = &expected;
                scope.spawn(move |_| {
                    for _ in 0..3 {
                        close(&evaluate(&**m, s, &[]), expected, false);
                    }
                });
            }
        })
    });
    struct Reentrant {
        potential: LennardJones,
        entered: AtomicBool,
        model: Arc<OnceLock<Weak<PairForceModel>>>,
        source: System,
        expected: (f64, f64, Vec<f64>),
    }
    impl PairPotential for Reentrant {
        fn cutoff(&self) -> f64 {
            2.5
        }
        fn cutoff2(&self) -> f64 {
            6.25
        }
        fn energy_force(&self, r2: f64, i: u32, j: u32) -> (f64, f64) {
            if !self.entered.swap(true, Ordering::SeqCst) {
                let m = self.model.get().unwrap().upgrade().unwrap();
                // Outer call owns scratch while this nested same-pool call
                // evaluates the same shared model. Blocking locks deadlock.
                rayon::join(
                    || close(&evaluate(&*m, &self.source, &[]), &self.expected, false),
                    || {},
                );
            }
            self.potential.energy_force(r2, i, j)
        }
    }
    let slot = Arc::new(OnceLock::new());
    let reentrant = Arc::new(PairForceModel::new(Reentrant {
        potential: lj(),
        entered: AtomicBool::new(false),
        model: slot.clone(),
        source: s.clone(),
        expected: expected.clone(),
    }));
    slot.set(Arc::downgrade(&reentrant)).unwrap();
    let (send, recv) = std::sync::mpsc::channel();
    std::thread::spawn(move || {
        pool.install(|| close(&evaluate(&*reentrant, &s, &[]), &expected, false));
        send.send(()).unwrap();
    });
    recv.recv_timeout(std::time::Duration::from_secs(20))
        .expect("nested Rayon shared-model call deadlocked");
}

use modena_core::{
    EwaldForceModel, ForceModel, NeighborList, SimBox, SpmeForceModel, System, Vec3, Virial,
};
use std::alloc::{GlobalAlloc, Layout, System as AllocSystem};
use std::sync::atomic::{AtomicUsize, Ordering};

struct CountingAllocator;
static ALLOCATIONS: AtomicUsize = AtomicUsize::new(0);
unsafe impl GlobalAlloc for CountingAllocator {
    unsafe fn alloc(&self, l: Layout) -> *mut u8 {
        ALLOCATIONS.fetch_add(1, Ordering::Relaxed);
        unsafe { AllocSystem.alloc(l) }
    }
    unsafe fn dealloc(&self, p: *mut u8, l: Layout) {
        unsafe { AllocSystem.dealloc(p, l) }
    }
    unsafe fn realloc(&self, p: *mut u8, l: Layout, n: usize) -> *mut u8 {
        ALLOCATIONS.fetch_add(1, Ordering::Relaxed);
        unsafe { AllocSystem.realloc(p, l, n) }
    }
}
#[global_allocator]
static ALLOCATOR: CountingAllocator = CountingAllocator;

fn system(edge: f64) -> System {
    let mut s = System::new(SimBox::cubic(edge));
    s.resize(3);
    for (i, p) in [
        Vec3::new(1.1, 2.3, 3.4),
        Vec3::new(2.7, 3.2, 3.8),
        Vec3::new(5.2, 5.8, 6.0),
    ]
    .into_iter()
    .enumerate()
    {
        s.set_pos(i, p);
    }
    s.charges = vec![1.0, -0.7, 0.4];
    s
}
fn evaluate(m: &dyn ForceModel, s: &System) -> (f64, f64, Vec<f64>) {
    let mut t = s.clone();
    t.zero_forces();
    let mut v = Virial::default();
    let u = m
        .compute(&mut t, &mut v, &mut NeighborList::with_skin(4.2, 0.3))
        .unwrap();
    (u, v.w, t.forces)
}
#[test]
fn spme_reuses_owned_work_buffers_and_shared_models_keep_box_states_separate() {
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(1)
        .build()
        .unwrap();
    let model = SpmeForceModel::new(0.8, 4.2, 6.4).unwrap();
    let mut s = system(9.0);
    let mut n = NeighborList::with_skin(4.2, 0.3);
    let mut v = Virial::default();
    pool.install(|| model.compute(&mut s, &mut v, &mut n).unwrap());
    s.zero_forces();
    v.w = 0.0;
    let before = ALLOCATIONS.load(Ordering::SeqCst);
    pool.install(|| model.compute(&mut s, &mut v, &mut n).unwrap());
    let warm = ALLOCATIONS.load(Ordering::SeqCst) - before;
    eprintln!("SPME warm compute allocations={warm}");
    assert_eq!(warm,0,"stable mesh/box single-thread SPME must reuse grid, plans, FFT scratch and influence factors");

    let pool4 = rayon::ThreadPoolBuilder::new()
        .num_threads(4)
        .build()
        .unwrap();
    let parallel_model = SpmeForceModel::new(0.8, 4.2, 6.4).unwrap();
    let mut large = system(20.0);
    let mut large_n = NeighborList::with_skin(4.2, 0.3);
    let mut large_v = Virial::default();
    pool4.install(|| {
        parallel_model
            .compute(&mut large, &mut large_v, &mut large_n)
            .unwrap()
    });
    large.zero_forces();
    large_v.w = 0.0;
    let before = ALLOCATIONS.load(Ordering::SeqCst);
    pool4.install(|| {
        parallel_model
            .compute(&mut large, &mut large_v, &mut large_n)
            .unwrap()
    });
    let parallel_allocations = ALLOCATIONS.load(Ordering::SeqCst) - before;
    eprintln!("SPME warm parallel compute allocations={parallel_allocations}");
    assert_eq!(
        parallel_allocations, 0,
        "worker scratch must also persist on parallel FFT meshes"
    );

    // The box changes but the 5-smooth mesh is unchanged; influence factors
    // still need invalidation. A change in kappa with the same mesh does too.
    let mut model = model;
    for (edge, kappa) in [(9.01, 0.8), (9.01, 0.9), (9.0, 0.8)] {
        model.kappa = kappa;
        let s = system(edge);
        let reference = SpmeForceModel::new(kappa, 4.2, 6.4).unwrap();
        assert_eq!(evaluate(&model, &s), evaluate(&reference, &s));
    }
    // Same per-instance workspace may be shared concurrently by simulations.
    // Serialized evaluation must hold the workspace's box key through use.
    for spme in [true, false] {
        let shared: std::sync::Arc<dyn ForceModel> = if spme {
            std::sync::Arc::new(SpmeForceModel::new(0.8, 4.2, 6.4).unwrap())
        } else {
            std::sync::Arc::new(EwaldForceModel::new(0.8, 4.2, 6.4).unwrap())
        };
        let expected: Vec<_> = [9.0, 11.0]
            .into_iter()
            .map(|edge| {
                let reference: Box<dyn ForceModel> = if spme {
                    Box::new(SpmeForceModel::new(0.8, 4.2, 6.4).unwrap())
                } else {
                    Box::new(EwaldForceModel::new(0.8, 4.2, 6.4).unwrap())
                };
                evaluate(&*reference, &system(edge))
            })
            .collect();
        std::thread::scope(|scope| {
            for (idx, edge) in [9.0, 11.0].into_iter().enumerate() {
                let m = &shared;
                let expected = &expected[idx];
                scope.spawn(move || {
                    let s = system(edge);
                    for _ in 0..12 {
                        assert_eq!(&evaluate(&**m, &s), expected);
                    }
                });
            }
        });
    }
}

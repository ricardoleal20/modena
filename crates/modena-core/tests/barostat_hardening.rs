use modena_core::{
    Barostat, ForceModel, LennardJones, MonteCarloBarostat, NeighborList, PairForceModel, SimBox,
    System, Vec3, Virial,
};

#[test]
fn ideal_gas_log_volume_moves_sample_gamma_volume_distribution() {
    let mut sum = 0.0;
    let mut sq = 0.0;
    let mut count = 0.0;
    for seed in [13, 417, 9091] {
        let mut s = System::new(SimBox::cubic(5.0_f64.cbrt()));
        s.resize(4);
        let mut b = MonteCarloBarostat::new(1.0, 1.0, 0.8, seed);
        for i in 0..120000 {
            b.attempt(&mut s, 0.0, |_| 0.0);
            if i >= 20000 {
                let v = s.box_size.volume();
                sum += v;
                sq += v * v;
                count += 1.0;
            }
        }
    }
    // Independent NPT ideal gas p(V) ∝ V^N exp(-βPV): Gamma(N+1, 1/βP).
    let mean = sum / count;
    let variance = sq / count - mean * mean;
    println!("MC ideal gas: samples={count} mean={mean:.12} variance={variance:.12}; expected mean=5 variance=5");
    assert!((mean - 5.0).abs() < 0.12, "mean={mean}, expected 5");
    assert!(
        (variance - 5.0).abs() < 0.35,
        "variance={variance}, expected 5"
    );
}

#[test]
fn accepted_and_rejected_model_moves_preserve_cache_integrity() {
    let model = PairForceModel::new(LennardJones::new(1.0, 1.0, 2.5).unwrap());
    let mut sys = System::new(SimBox::cubic(6.0));
    sys.resize(2);
    sys.set_pos(0, Vec3::new(1.0, 1.0, 1.0));
    sys.set_pos(1, Vec3::new(2.2, 1.0, 1.0));
    let mut neighbors = NeighborList::with_skin(2.5, 0.7);
    let mut v = Virial::default();
    let mut pe = model.compute(&mut sys, &mut v, &mut neighbors).unwrap();
    let mut b = MonteCarloBarostat::new(0.02, 1.0, 0.5, 51);
    let mut accepted = 0;
    let mut rejected = 0;
    for _ in 0..100 {
        let old = sys.clone();
        let old_pe = pe;
        let old_w = v.w;
        if b.apply_with_model(&mut sys, &model, &mut neighbors, &mut pe, &mut v, 0.001)
            .unwrap()
        {
            accepted += 1;
        } else {
            rejected += 1;
            assert_eq!(sys.positions, old.positions);
            assert_eq!(sys.velocities, old.velocities);
            assert_eq!(sys.forces, old.forces);
            assert_eq!(sys.box_size, old.box_size);
            assert_eq!(pe, old_pe);
            assert_eq!(v.w, old_w);
        }
        let mut fresh = sys.clone();
        fresh.zero_forces();
        let mut w = Virial::default();
        let reference = model
            .compute(&mut fresh, &mut w, &mut NeighborList::new(2.5))
            .unwrap();
        assert_eq!(fresh.forces, sys.forces);
        assert_eq!(reference, pe);
        assert_eq!(w.w, v.w);
        // Reuse the actual post-rejection list, including a large skin.
        fresh.zero_forces();
        w = Virial::default();
        assert_eq!(
            model.compute(&mut fresh, &mut w, &mut neighbors).unwrap(),
            pe
        );
        assert_eq!(fresh.forces, sys.forces);
    }
    assert!(
        accepted > 0 && rejected > 0,
        "accepted={accepted} rejected={rejected}"
    );
}

#[test]
fn interval_and_rng_state_continue_identically() {
    let model = PairForceModel::new(LennardJones::new(1.0, 0.0, 2.5).unwrap());
    let mut a = System::new(SimBox::cubic(5.0_f64.cbrt()));
    a.resize(4);
    let mut b = a.clone();
    let mut first = MonteCarloBarostat::new(1.0, 1.0, 0.5, 91)
        .with_interval(3)
        .unwrap();
    let mut n = NeighborList::new(2.5);
    let mut pe = 0.0;
    let mut w = Virial::default();
    for _ in 0..7 {
        first
            .apply_with_model(&mut a, &model, &mut n, &mut pe, &mut w, 0.1)
            .unwrap();
    }
    assert_eq!(first.attempted, 2);
    b.clone_from(&a);
    let state = first.checkpoint_state().unwrap();
    let mut second = MonteCarloBarostat::new(1.0, 1.0, 0.5, 91)
        .with_interval(3)
        .unwrap();
    second.restore_state(&state).unwrap();
    let mut nn = NeighborList::new(2.5);
    let mut pp = pe;
    let mut ww = w;
    for _ in 0..23 {
        let aa = first
            .apply_with_model(&mut a, &model, &mut n, &mut pe, &mut w, 0.1)
            .unwrap();
        let bb = second
            .apply_with_model(&mut b, &model, &mut nn, &mut pp, &mut ww, 0.1)
            .unwrap();
        assert_eq!(aa, bb);
        assert_eq!(a.box_size, b.box_size);
        assert_eq!(a.positions, b.positions);
    }
    assert_eq!(first.attempted, 10);
    assert_eq!(
        first.checkpoint_state().unwrap(),
        second.checkpoint_state().unwrap()
    );
}

#[test]
fn failed_trial_rolls_back_mechanics_and_neighbor_state() {
    struct Fail;
    impl ForceModel for Fail {
        fn compute(
            &self,
            s: &mut System,
            w: &mut Virial,
            _n: &mut NeighborList,
        ) -> modena_core::Result<f64> {
            s.forces.fill(123.0);
            w.w = 999.0;
            Err(modena_core::ModenaError::Physics(
                "deliberate trial failure".into(),
            ))
        }
    }
    let mut s = System::new(SimBox::cubic(3.0));
    s.resize(2);
    s.set_pos(1, Vec3::new(1.0, 0.0, 0.0));
    s.forces.fill(0.7);
    let old = s.clone();
    let mut n = NeighborList::new(2.5);
    n.build(&s);
    let pairs = n.pairs().to_vec();
    let mut p = 1.0;
    let mut w = Virial { w: 2.0 };
    let mut b = MonteCarloBarostat::new(1.0, 1.0, 0.1, 42);
    assert!(b
        .apply_with_model(&mut s, &Fail, &mut n, &mut p, &mut w, 0.1)
        .is_err());
    assert_eq!(s.positions, old.positions);
    assert_eq!(s.forces, old.forces);
    assert_eq!(s.box_size, old.box_size);
    assert_eq!(p, 1.0);
    assert_eq!(w.w, 2.0);
    assert_eq!(n.pairs(), pairs);
}

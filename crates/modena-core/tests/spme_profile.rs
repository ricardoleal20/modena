//! Profiling must preserve the same analytic mesh functional and additive contract.
use modena_core::stress::TensorAvailability;
use modena_core::{
    force_model::{ForceModel, Virial},
    neighbor::{NeighborList, PairScale},
    spme::{SpmeForceModel, SpmeProfile},
    system::{SimBox, System},
};
use std::time::Duration;

fn fields(p: &SpmeProfile) -> [Duration; 12] {
    [
        p.neighbor,
        p.real_space,
        p.pair_correction,
        p.workspace_setup,
        p.influence_setup,
        p.spread,
        p.fft_forward,
        p.influence,
        p.fft_inverse,
        p.gather,
        p.self_background,
        p.total,
    ]
}
fn close(a: f64, b: f64) {
    assert!((a - b).abs() <= 1e-12 * b.abs().max(1.), "{a} != {b}");
}
fn special() -> (System, Vec<PairScale>) {
    let mut s = System::new(SimBox::from_array([9., 10., 11.]));
    s.resize(3);
    s.positions = vec![1.1, 2.3, 3.4, 2.7, 3.2, 3.8, 6.2, 5.8, 7.];
    s.charges = vec![1., -0.7, 0.4];
    (
        s,
        vec![
            PairScale {
                i: 0,
                j: 1,
                lj: 1.,
                coulomb: 0.,
            },
            PairScale {
                i: 0,
                j: 2,
                lj: 1.,
                coulomb: 0.4,
            },
        ],
    )
}
fn ionic() -> (System, Vec<PairScale>) {
    let mut s = modena_core::lattice::rocksalt(2.2, [2; 3], "Na", "Cl", 1., 1., 1., -1.);
    for (i, r) in s.positions.chunks_exact_mut(3).enumerate() {
        for (axis, v) in r.iter_mut().enumerate() {
            let z = ((i as u64 + 1) * 1_664_525 + (axis as u64 + 1) * 1_013_904_223 + 19 * 69_069)
                % 65_536;
            *v = (v.rem_euclid(4.4) + 0.05 * 2.2 * (2. * z as f64 / 65_535. - 1.)).rem_euclid(4.4);
        }
    }
    s.charges[1] = -0.5;
    (s, vec![])
}
fn parity(mut initial: System, scales: Vec<PairScale>) {
    let rc = 0.4
        * initial
            .box_size
            .size
            .0
            .into_iter()
            .fold(f64::INFINITY, f64::min);
    let mut ordinary = SpmeForceModel::new(1. / 2.2, rc, 12. / 2.2)
        .unwrap()
        .with_pair_scales(scales.clone())
        .unwrap();
    let mut profiled = SpmeForceModel::new(1. / 2.2, rc, 12. / 2.2)
        .unwrap()
        .with_pair_scales(scales.clone())
        .unwrap();
    let mut nl = NeighborList::with_pair_scales(rc, 0.3, &scales).unwrap();
    let mut np = NeighborList::with_pair_scales(rc, 0.3, &scales).unwrap();
    for step in 0..7 {
        match step {
            2 => initial.positions[0] += 0.001,
            3 => {
                ordinary.kappa *= 1.1;
                profiled.kappa *= 1.1;
            }
            4 => {
                initial.box_size = SimBox::from_array(initial.box_size.size.0.map(|x| x * 1.01));
            }
            5 => {
                ordinary.kspace_cutoff *= 1.3;
                profiled.kspace_cutoff *= 1.3;
            }
            _ => {}
        }
        let mut a = initial.clone();
        let mut b = initial.clone();
        let mut zero = initial.clone();
        zero.zero_forces();
        for (i, f) in a.forces.iter_mut().enumerate() {
            *f = 0.125 + (i as f64 * 0.3).sin();
        }
        b.forces.clone_from(&a.forces);
        let preload = a.forces.clone();
        let mut wa = Virial { w: 1.75 };
        let mut wb = Virial { w: 1.75 };
        let ua = ordinary.compute(&mut a, &mut wa, &mut nl).unwrap();
        let (ub, p) = profiled.compute_profiled(&mut b, &mut wb, &mut np).unwrap();
        close(ua, ub);
        close(wa.w, wb.w);
        for (&x, &y) in a.forces.iter().zip(&b.forces) {
            close(x, y);
        }
        let mut wz = Virial::default();
        let (uz, _) = profiled
            .compute_profiled(&mut zero, &mut wz, &mut np)
            .unwrap();
        close(ub, uz);
        close(wb.w, wz.w + 1.75);
        for ((&x, &y), &pre) in b.forces.iter().zip(&zero.forces).zip(&preload) {
            close(x, y + pre);
        }
        // Tensor and legacy entry points must run the same physics exactly once.
        let mut evaluated = initial.clone();
        evaluated.forces.clone_from(&preload);
        let ea = ordinary.evaluate(&mut evaluated, &mut nl).unwrap();
        let mut evaluated_profiled = initial.clone();
        evaluated_profiled.forces.clone_from(&preload);
        let (eb, tensor_profile) = profiled
            .evaluate_profiled(&mut evaluated_profiled, &mut np)
            .unwrap();
        assert_eq!(ea, eb);
        close(ea.energy, ua);
        close(ea.scalar_virial + 1.75, wa.w);
        let TensorAvailability::Complete(tensor) = ea.tensor else {
            panic!("SPME must provide its full tensor");
        };
        close(tensor.trace(), ea.scalar_virial);
        assert_eq!(evaluated.forces, a.forces);
        assert_eq!(evaluated_profiled.forces, a.forces);
        // Fresh models catch stale shared influence factors after cell/kappa/grid changes.
        let fresh = SpmeForceModel::new(ordinary.kappa, rc, ordinary.kspace_cutoff)
            .unwrap()
            .with_pair_scales(scales.clone())
            .unwrap();
        let mut fresh_system = initial.clone();
        fresh_system.forces.clone_from(&preload);
        let fresh_e = fresh
            .evaluate(
                &mut fresh_system,
                &mut NeighborList::with_pair_scales(rc, 0.3, &scales).unwrap(),
            )
            .unwrap();
        assert_eq!(fresh_e, eb);
        assert_eq!(fresh_system.forces, a.forces);
        assert_eq!(fields(&tensor_profile).len(), 12);
        assert_eq!(fields(&p).len(), 12); // No timing thresholds or phase-sum equality.
        assert_eq!(nl.rebuild_count(), np.rebuild_count());
        if step == 1 {
            assert_eq!(np.rebuild_count(), 1);
        }
    }
}
#[test]
fn charged_scaled_three_ion_profile_parity() {
    let (s, p) = special();
    parity(s, p);
}
#[test]
fn perturbed_charged_sixty_four_ion_profile_parity() {
    let (s, p) = ionic();
    assert_eq!(s.n(), 64);
    parity(s, p);
}

#[test]
fn parallel_fft_tensor_profile_and_reused_workspace_match_one_worker() {
    let (mut source, scales) = special();
    source.box_size = SimBox::from_array([20.; 3]);
    let model = SpmeForceModel::new(0.8, 3.5, 7.5)
        .unwrap()
        .with_pair_scales(scales.clone())
        .unwrap();
    assert_eq!(
        modena_core::spme::mesh_for_box(&source.box_size, 7.5),
        [48; 3]
    );
    let mut expected = None;
    // 48³ crosses the production parallel-FFT threshold. Reuse the same
    // workspace across pools so growing worker scratch is also exercised.
    for threads in [1, 4] {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .unwrap();
        let mut s = source.clone();
        s.forces.fill(0.125);
        let mut nl = NeighborList::with_pair_scales(3.5, 0.3, &scales).unwrap();
        let e = pool.install(|| model.evaluate(&mut s, &mut nl).unwrap());
        let TensorAvailability::Complete(t) = e.tensor else {
            panic!("missing tensor");
        };
        close(t.trace(), e.scalar_virial);
        let mut p = source.clone();
        p.forces.fill(0.125);
        let (pe, profile) = pool.install(|| model.evaluate_profiled(&mut p, &mut nl).unwrap());
        assert_eq!(e, pe);
        assert_eq!(s.forces, p.forces);
        assert_eq!(fields(&profile).len(), 12);
        if let Some((ref old_e, ref old_f)) = expected {
            assert_eq!(&e, old_e);
            assert_eq!(&s.forces, old_f);
        } else {
            expected = Some((e, s.forces));
        }
    }
}

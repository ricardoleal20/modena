//! Ensemble & sampling integration tests (M5).
//!
//! Three correctness checks: (1) the Nosé–Hoover chain conserves its extended
//! Hamiltonian H* = KE + U + Σ p_ξk²/(2Q_k) + N_f kT ξ_1 + Σ kT ξ_k (the shadow
//! quantity whose drift exposes a wrong Suzuki–Yoshida decomposition or a sign
//! error in the chain propagation); (2) replica exchange runs the parallel-
//! tempering ladder end-to-end with finite energies and non-zero swap acceptance;
//! (3) a Berendsen barostat configured via TOML reaches `Simulation::run` and
//! actually moves the box volume (exercises `build_barostat` + `run_streaming`).

use modena_core::{
    lattice, ForceModel, Integrator, LennardJones, NeighborList, NoseHooverChain, PairForceModel,
    PairPotential, ReplicaExchange, System, Vec3, Virial,
};

fn build_lj() -> (Box<dyn ForceModel>, f64) {
    let lj = LennardJones::new(1.0, 1.0, 2.5).unwrap();
    let cutoff = lj.cutoff();
    (Box::new(PairForceModel::new(lj)), cutoff)
}

/// Maxwellian velocities at `temperature`, COM drift removed, rescaled to the
/// exact target (mirrors the runner's `init_velocities`).
fn thermalize(system: &mut System, temperature: f64, seed: u64) {
    use rand::rngs::StdRng;
    use rand::SeedableRng;
    use rand_distr::{Distribution, StandardNormal};
    let mut rng = StdRng::seed_from_u64(seed);
    let n = system.n();
    let normal = StandardNormal;
    let mut p_com = Vec3::zero();
    let mut total_m = 0.0;
    for i in 0..n {
        let m = system.masses[i];
        let s = (temperature / m).sqrt();
        let gx: f64 = normal.sample(&mut rng);
        let gy: f64 = normal.sample(&mut rng);
        let gz: f64 = normal.sample(&mut rng);
        let v = Vec3::new(s * gx, s * gy, s * gz);
        system.set_vel(i, v);
        p_com += v * m;
        total_m += m;
    }
    let v_com = p_com * (1.0 / total_m);
    for i in 0..n {
        system.set_vel(i, system.vel(i) - v_com);
    }
    let ke = system.kinetic_energy();
    let dof = system.dof() as f64;
    if ke > 0.0 && dof > 0.0 {
        let current_t = 2.0 * ke / dof;
        if current_t > 0.0 {
            let lambda = (temperature / current_t).sqrt();
            for i in 0..n {
                system.set_vel(i, system.vel(i) * lambda);
            }
        }
    }
}

#[test]
fn nhc_extended_hamiltonian_is_conserved() {
    let (model, cutoff) = build_lj();
    let mut system = lattice::fcc(1.6796, [2, 2, 2], "Ar"); // 32 atoms
    assert_eq!(system.n(), 32);
    thermalize(&mut system, 1.0, 42);

    let mut neighbors = NeighborList::new(cutoff);
    let mut virial = Virial::default();
    let pe0 = model
        .compute(&mut system, &mut virial, &mut neighbors)
        .unwrap();

    let mut nhc = NoseHooverChain::new(1.0, 0.5, 4);
    let dt = 0.0025; // small dt → tight shadow-Hamiltonian conservation
    let steps = 4000;

    let h0 = nhc.extended_hamiltonian(&system, pe0);
    let mut h_min = h0;
    let mut h_max = h0;
    let mut t_sum = 0.0;
    for _ in 0..steps {
        let e = nhc.step(&mut system, &*model, &mut neighbors, dt).unwrap();
        let h = nhc.extended_hamiltonian(&system, e.potential);
        h_min = h_min.min(h);
        h_max = h_max.max(h);
        t_sum += e.temperature;
    }
    let drift = (h_max - h_min).abs() / h0.abs();
    let mean_t = t_sum / steps as f64;
    // The NHC shadow Hamiltonian should be conserved to well under 1% over a
    // short run; a sign/order error in the chain propagation drifts secularly.
    assert!(
        drift < 1e-2,
        "NHC extended-Hamiltonian drift too large: {drift:.3e} (H0={h0:.6}, range=[{h_min:.6},{h_max:.6}])"
    );
    // The chain must actually thermostat: mean T within 30% of the 1.0 target.
    assert!(
        (mean_t - 1.0).abs() < 0.3,
        "NHC mean temperature {mean_t:.4} far from target 1.0"
    );
}

#[test]
fn replica_exchange_runs_and_swaps() {
    let (model, _cutoff) = build_lj();
    let system = lattice::fcc(1.6796, [2, 2, 2], "Ar"); // 32 atoms

    let temperatures = vec![0.5, 0.8, 1.2, 2.0];
    let mut rex = ReplicaExchange::new(&system, temperatures, 0.005, 1.0, 42).unwrap();
    let stats = rex
        .run(&*model, 500, 25, 0.005)
        .expect("replica-exchange run must complete");

    assert!(
        stats.swaps_attempted > 0,
        "expected swap attempts, got {}",
        stats.swaps_attempted
    );
    assert_eq!(
        stats.final_energies.len(),
        4,
        "one final energy per replica"
    );
    for &e in &stats.final_energies {
        assert!(e.is_finite(), "non-finite replica energy: {e}");
    }
    // A healthy parallel-tempering ladder accepts *some* swaps (a zero acceptance
    // rate means the ladder is broken or the criterion is inverted). We don't
    // require a high rate — just that the path is not fully blocked.
    assert!(
        stats.swaps_accepted > 0,
        "no swaps accepted (attempted={}, accepted={}) — ladder may be broken",
        stats.swaps_attempted,
        stats.swaps_accepted
    );
}

#[test]
fn replica_exchange_uses_actual_cutoff_above_two_point_five() {
    let mut s = System::new(modena_core::SimBox::cubic(12.0));
    s.resize(2);
    s.set_pos(1, Vec3::new(3.2, 0.0, 0.0));
    let pot = LennardJones::new(1.0, 1.0, 3.7).unwrap();
    let model = PairForceModel::new(pot.clone());
    let mut rex = ReplicaExchange::new(&s, vec![1e-4, 2e-4], 1e-4, 0.0, 42).unwrap();
    let stats = rex.run(&model, 2, 10, 1e-4).unwrap();
    for (s, u) in rex.systems().into_iter().zip(stats.final_energies) {
        let d = s.box_size.min_image(s.pos(0) - s.pos(1));
        assert!(d.norm() > 2.5 && d.norm() < 3.7);
        let (expected, fr) = pot.energy_force(d.norm2(), 0, 0);
        assert!(expected.abs() > 1e-5);
        assert!(
            (u - expected).abs() < 1e-14,
            "replica neighbor list must include r > 2.5 pairs"
        );
        for k in 0..3 {
            assert!((s.force(0).0[k] - d.0[k] * fr).abs() < 1e-14);
        }
    }
}

#[test]
fn runner_applies_berendsen_barostat_and_moves_volume() {
    use modena_core::{InputConfig, Simulation};
    use std::path::Path;

    // NVT (NHC) + a Berendsen barostat targeting a high pressure (P0 = 5.0) →
    // the box must compress over the run. This exercises `build_barostat` and
    // the `barostat.apply(...)` call in `run_streaming`.
    let toml = r#"
units = "reduced"

[structure]
type = "fcc"
lattice_constant = 1.6796
replicas = [2, 2, 2]
element = "Ar"

[model]
type = "lennard_jones"
sigma = 1.0
epsilon = 1.0
cutoff = 2.5

[integrator]
type = "nose_hoover_chain"
dt = 0.005
target_temperature = 1.0
tau = 0.5
chain_len = 4

[barostat]
type = "berendsen"
target_pressure = 5.0
tau = 0.5
compressibility = 1.0

[run]
steps = 400
initial_temperature = 1.0
seed = 42
"#;
    let config = InputConfig::from_str(toml).expect("config must parse");
    let base_dir = Path::new(".");
    let mut sim = Simulation::from_config(&config, base_dir).expect("sim must build");

    let v0 = sim.system().box_size.volume();
    let summary = sim.run().expect("run must complete");
    let v1 = sim.system().box_size.volume();

    // A meaningful pressure-coupling step changes the volume. With P0 = 5.0
    // (high) the box should compress; either way a no-op barostat (the
    // regression we're guarding against) would leave v1 == v0 exactly.
    assert!(
        (v1 - v0).abs() > 1e-9,
        "barostat did not move the box: v0={v0}, v1={v1}"
    );
    // Sanity: the run produced a conservation report.
    assert!(summary.steps > 0);
}

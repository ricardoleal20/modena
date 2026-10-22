//! Independent rigid-body regressions: free rotation cannot dissipate, and
//! three nonlinear rotational DOF have canonical K_rot~Gamma(3/2,T).
use modena_core::constraints::RigidWater;
use modena_core::{
    Constraints, ForceModel, Integrator, LangevinBAOAB, NeighborList, Result, SimBox, System, Vec3,
    VelocityVerlet, Virial,
};

struct ZeroForce;
impl ForceModel for ZeroForce {
    fn compute(&self, _: &mut System, _: &mut Virial, _: &mut NeighborList) -> Result<f64> {
        Ok(0.0)
    }
}

fn equilateral_rotor() -> (System, Constraints) {
    let mut s = System::new(SimBox::cubic(100.0));
    s.resize(3);
    let center = Vec3::splat(50.0);
    for i in 0..3 {
        let angle = 2.0 * std::f64::consts::PI * i as f64 / 3.0;
        let relative = Vec3::new(angle.cos(), angle.sin(), 0.0);
        s.set_pos(i, center + relative);
        s.set_vel(i, Vec3::new(-relative.0[1], relative.0[0], 0.0));
    }
    let length = 3.0_f64.sqrt();
    let c = Constraints {
        bonds: vec![],
        waters: vec![RigidWater {
            first: 0,
            d_oh: length,
            d_hh: length,
        }],
    };
    (s, c)
}

#[test]
fn rigid_water_free_rotor_conserves_kinetic_energy_across_timesteps() {
    let mut worst = 0.0_f64;
    for dt in [0.01_f64, 0.005, 0.001] {
        for baoab in [true, false] {
            let (mut s, c) = equilateral_rotor();
            let mut integrator: Box<dyn Integrator> = if baoab {
                Box::new(LangevinBAOAB::with_constraints(0.0, 1.0, 7, c.clone()))
            } else {
                Box::new(VelocityVerlet::with_constraints(c.clone()))
            };
            let initial = s.kinetic_energy();
            let mut neighbors = NeighborList::new(2.5);
            let mut max_drift = 0.0_f64;
            for _ in 0..(10.0 / dt).round() as usize {
                integrator
                    .step(&mut s, &ZeroForce, &mut neighbors, dt)
                    .unwrap();
                max_drift = max_drift.max((s.kinetic_energy() / initial - 1.0).abs());
            }
            assert!(c.max_residual(&s) < 1e-10);
            eprintln!("rotor baoab={baoab} dt={dt} initial_K={initial:.15} final_K={:.15} max_relative_drift={max_drift:e}",s.kinetic_energy());
            worst = worst.max(max_drift);
        }
    }
    // Rotation about a principal axis of this free rigid triangle keeps K
    // exactly constant. A timestep-proportional dissipative drift is forbidden.
    assert!(
        worst < 1e-7,
        "free rigid rotor lost kinetic energy, worst relative drift={worst:e}"
    );
}

fn rotational_energy(system: &System) -> f64 {
    let total_mass: f64 = system.masses.iter().sum();
    let center_velocity = system.total_momentum() * (1.0 / total_mass);
    (0..system.n())
        .map(|i| 0.5 * system.masses[i] * (system.vel(i) - center_velocity).norm2())
        .sum()
}

#[test]
fn rigid_water_rotational_modes_have_canonical_equipartition() {
    const SAMPLES: usize = 100_000;
    const SEEDS: [u64; 3] = [7, 42, 91];
    for (gamma, dt) in [(0.1_f64, 0.04_f64), (2.0, 0.01)] {
        let (mut sum, mut sum2) = (0.0, 0.0);
        for seed in SEEDS {
            let (mut s, c) = equilateral_rotor();
            s.velocities.fill(0.0);
            let mut integrator = LangevinBAOAB::with_constraints(gamma, 1.0, seed, c);
            let mut neighbors = NeighborList::new(2.5);
            let burn = (20.0 / (gamma * dt)).ceil() as usize; // twenty friction times
            for step in 0..burn + SAMPLES {
                integrator
                    .step(&mut s, &ZeroForce, &mut neighbors, dt)
                    .unwrap();
                if step >= burn {
                    let k = rotational_energy(&s);
                    sum += k;
                    sum2 += k * k;
                }
            }
        }
        let samples = (SAMPLES * SEEDS.len()) as f64;
        let mean = sum / samples;
        let variance = sum2 / samples - mean * mean;
        // Gaussian tangent momenta give three rotational quadratic modes:
        // independent canonical mean=variance=3/2 at T=1. Account for the
        // slowest OU energy correlation c^2=exp(-2 gamma dt). The centered
        // squared-energy estimator's canonical variance is 13.5; using the
        // same correlation factor is conservative (its second mode is faster).
        let correlation = (-2.0 * gamma * dt).exp();
        let inflation = (1.0 + correlation) / (1.0 - correlation);
        let mean_tolerance = 5.0 * (1.5 * inflation / samples).sqrt() + 0.01;
        let variance_tolerance = 5.0 * (13.5 * inflation / samples).sqrt() + 0.01;
        eprintln!("rigid rotation gamma={gamma} dt={dt} samples={samples} mean={mean} variance={variance} canonical_mean=1.5 canonical_variance=1.5 mean_tolerance={mean_tolerance} variance_tolerance={variance_tolerance}");
        assert!(
            (mean - 1.5).abs() < mean_tolerance,
            "rotational equipartition failed at gamma={gamma},dt={dt}: mean={mean}"
        );
        assert!(
            (variance - 1.5).abs() < variance_tolerance,
            "rotational canonical variance failed at gamma={gamma},dt={dt}: variance={variance}"
        );
    }
}

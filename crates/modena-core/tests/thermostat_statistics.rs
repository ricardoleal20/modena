//! Independent canonical-distribution checks of the public CSVR thermostat.

use modena_core::{SimBox, System, Thermostat, VelocityRescaling};

fn system_with_kinetic_energy(atoms: usize, kinetic_energy: f64) -> System {
    let mut system = System::new(SimBox::cubic(10.0));
    system.resize(atoms);
    // Unequal masses exercise KE-weighted rescaling without removing modes.
    for (i, mass) in system.masses.iter_mut().enumerate() {
        *mass = 1.0 + 0.25 * i as f64;
    }
    let squared_speed = 2.0 * kinetic_energy / system.masses.iter().sum::<f64>();
    system.velocities.fill((squared_speed / 3.0).sqrt());
    system
}

#[test]
fn bussi_one_step_matches_noncentral_chi_squared_moments() {
    // The exact transition is s * chi_squared(Nf, a*K/s), independently
    // of how the implementation generates its random variables.
    const SAMPLES: usize = 50_000;
    let kinetic_energy = 2.0;
    let target_kinetic_energy = 3.5;
    for ndof in [1_usize, 2, 3, 6, 60] {
        let mut system = system_with_kinetic_energy(ndof.div_ceil(3), kinetic_energy);
        // Isolate the thermostat's DOF contract; constraint projection is
        // covered separately. At Nf=1 and 2 only active coordinates move.
        system.set_constraint_count(3 * system.n() - ndof).unwrap();
        system.velocities[ndof..].fill(0.0);
        let normalization = (kinetic_energy / system.kinetic_energy()).sqrt();
        for velocity in &mut system.velocities {
            *velocity *= normalization;
        }
        let initial_velocities = system.velocities.clone();
        let ndof = ndof as f64;
        for dt in [0.05_f64, 1.0, 8.0] {
            let a = (-dt).exp();
            let scale = (1.0 - a) * target_kinetic_energy / ndof;
            let expected_mean = a * kinetic_energy + (1.0 - a) * target_kinetic_energy;
            let expected_variance = 2.0 * ndof * scale * scale + 4.0 * a * kinetic_energy * scale;
            let fourth_cumulant =
                48.0 * ndof * scale.powi(4) + 192.0 * a * kinetic_energy * scale.powi(3);
            let mean_tolerance = 7.0 * (expected_variance / SAMPLES as f64).sqrt();
            let variance_tolerance =
                7.0 * ((fourth_cumulant + 2.0 * expected_variance.powi(2)) / SAMPLES as f64).sqrt();
            for seed in [7, 31] {
                let mut thermostat =
                    VelocityRescaling::new(2.0 * target_kinetic_energy / ndof, 1.0, seed);
                let mut mean = 0.0;
                let mut squared_deviations = 0.0;
                for sample in 1..=SAMPLES {
                    system.velocities.copy_from_slice(&initial_velocities);
                    thermostat.apply(&mut system, dt);
                    let new_kinetic_energy = system.kinetic_energy();
                    let delta = new_kinetic_energy - mean;
                    mean += delta / sample as f64;
                    squared_deviations += delta * (new_kinetic_energy - mean);
                }
                let variance = squared_deviations / (SAMPLES - 1) as f64;
                assert!(
                    (mean - expected_mean).abs() < mean_tolerance,
                    "Nf={ndof}, dt/tau={dt}, seed={seed}: mean {mean}, expected {expected_mean} +/- {mean_tolerance}"
                );
                assert!(
                    (variance - expected_variance).abs() < variance_tolerance,
                    "Nf={ndof}, dt/tau={dt}, seed={seed}: variance {variance}, expected {expected_variance} +/- {variance_tolerance}"
                );
            }
        }
    }
}

#[test]
fn bussi_stationary_kinetic_energy_has_canonical_mean_and_variance() {
    const SAMPLES: usize = 30_000;
    let temperature = 1.7;
    // Cover weak coupling, dt=tau, and almost independent bath updates.
    for (atoms, dt, tau) in [(1, 0.025, 0.25), (4, 0.05, 0.05), (16, 0.05, 0.01)] {
        for seed in [17, 101, 9929] {
            let mut system = system_with_kinetic_energy(atoms, 0.1);
            let ndof = system.dof() as f64;
            let expected_mean = 0.5 * ndof * temperature;
            let expected_variance = 0.5 * ndof * temperature * temperature;
            let mut thermostat = VelocityRescaling::new(temperature, tau, seed);
            // Relax for 20 tau, then separate samples by at least 5 tau.
            // The exact KE autocorrelation is exp(-t/tau), below 0.007 here.
            let burn_in = (20.0_f64 * tau / dt).ceil() as usize;
            let stride = (5.0_f64 * tau / dt).ceil().max(1.0) as usize;
            for _ in 0..burn_in {
                thermostat.apply(&mut system, dt);
            }
            let mut mean = 0.0;
            let mut squared_deviations = 0.0;
            for sample in 1..=SAMPLES {
                for _ in 0..stride {
                    thermostat.apply(&mut system, dt);
                }
                let kinetic_energy = system.kinetic_energy();
                let delta = kinetic_energy - mean;
                mean += delta / sample as f64;
                squared_deviations += delta * (kinetic_energy - mean);
            }
            let variance = squared_deviations / (SAMPLES - 1) as f64;
            // Seven standard errors, with a 2% allowance for residual
            // correlation. Gamma(Nf/2,T) has excess kurtosis 12/Nf.
            let mean_tolerance = 7.0 * (1.02 * expected_variance / SAMPLES as f64).sqrt();
            let variance_tolerance =
                7.0 * expected_variance * (1.02 * (2.0 + 12.0 / ndof) / SAMPLES as f64).sqrt();
            assert!(
                (mean - expected_mean).abs() < mean_tolerance,
                "Nf={ndof}, dt={dt}, tau={tau}, seed={seed}: mean KE {mean}, canonical {expected_mean} +/- {mean_tolerance}"
            );
            assert!(
                (variance - expected_variance).abs() < variance_tolerance,
                "Nf={ndof}, dt={dt}, tau={tau}, seed={seed}: KE variance {variance}, canonical {expected_variance} +/- {variance_tolerance}"
            );
        }
    }
}

#[test]
fn bussi_decorrelated_rescaling_reverses_velocities_with_equal_probability() {
    const SAMPLES: usize = 50_000;
    let mut system = system_with_kinetic_energy(1, 1.5);
    // The single-DOF limit must produce a signed Maxwellian velocity.
    system.set_constraint_count(2).unwrap();
    system.velocities[1..].fill(0.0);
    system.velocities[0] = 1.0;
    let initial_velocities = system.velocities.clone();
    let mut thermostat = VelocityRescaling::new(1.0, 1.0, 238);
    let mut reversals = 0;
    let mut sum = 0.0;
    let mut sum2 = 0.0;
    for _ in 0..SAMPLES {
        system.velocities.copy_from_slice(&initial_velocities);
        thermostat.apply(&mut system, 100.0);
        reversals += usize::from(system.velocities[0] < 0.0);
        sum += system.velocities[0];
        sum2 += system.velocities[0] * system.velocities[0];
    }
    let fraction = reversals as f64 / SAMPLES as f64;
    assert!(
        (fraction - 0.5).abs() < 0.02,
        "complete bath decorrelation reversed {fraction} of velocities; expected 0.5"
    );
    assert!((sum / SAMPLES as f64).abs() < 0.035);
    assert!((sum2 / SAMPLES as f64 - 1.0).abs() < 0.05);
}

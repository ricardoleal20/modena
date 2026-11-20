#[path = "../examples/stress_validation/mod.rs"]
mod validation;
use validation::protocol::CaseStatus;
use validation::statistics::estimate_replicates;

#[test]
fn actual_rotor_pressure_campaign_passes_analytic_gates() {
    let dir = tempfile::tempdir().unwrap();
    let gates = validation::fixtures::run_deterministic(dir.path()).unwrap();
    if let Some(path) = std::env::var_os("B07_ENGINEERING_ROWS") {
        let path = std::path::Path::new(&path);
        validation::protocol::write_new(&path.join("engineering-gates.json"), &gates).unwrap();
        std::fs::copy(
            dir.path().join("deterministic.jsonl"),
            path.join("engineering-deterministic.jsonl"),
        )
        .unwrap();
    }
    for gate in &gates {
        assert_eq!(gate.status, CaseStatus::Passed, "{gate:?}");
    }
    assert!(gates
        .iter()
        .any(|g| g.name == "vv_rotor_mixed_time_negative_control" && g.observed[0] < -1e-4));
    assert!(gates
        .iter()
        .any(|g| g.name.contains("nonprincipal") && g.name.contains("interval")));
}

#[test]
fn replica_uncertainty_and_positive_pair_blocking() {
    let series = |repeat: usize| {
        (0..8)
            .map(|seed| {
                (0..8000)
                    .map(|j| {
                        seed as f64
                            + if (j / repeat).is_multiple_of(2) {
                                -1.
                            } else {
                                1.
                            }
                    })
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>()
    };
    let fast = estimate_replicates(&series(1), 0.1).unwrap();
    let slow = estimate_replicates(&series(8), 0.1).unwrap();
    assert!((fast.mean - 3.5).abs() < 1e-12);
    assert!((fast.se - (0.75_f64).sqrt()).abs() < 1e-12);
    assert!(slow.block_samples > fast.block_samples && slow.n_eff < fast.n_eff);
    assert!(slow.ci95[0] < slow.mean && slow.ci95[1] > slow.mean);
    let width = |repeat| {
        let rows = series(repeat);
        let t = validation::statistics::first_half_tau(&rows).unwrap();
        let ci = validation::statistics::infer(&rows, (10. * t).ceil() as usize, t)
            .unwrap()
            .bootstrap_ci95
            .unwrap();
        ci[1] - ci[0]
    };
    assert!(width(8) > width(1));
    let mut invalid = series(8);
    invalid[3][40] = f64::NAN;
    assert!(estimate_replicates(&invalid, 0.1).is_err());
}

#[test]
fn prepared_snapshots_and_sampling_clocks_match_the_registered_roster() {
    use validation::ensemble::sample_step;
    let f = validation::matched::fixtures().unwrap();
    let base = f.iter().find(|f| f.id == "interacting32").unwrap();
    assert_eq!(base.positions.len(), 96);
    let mut prior = &base.velocities;
    for seed in 2001..=2004 {
        let state = f
            .iter()
            .find(|f| f.id == format!("interacting32-{seed}"))
            .unwrap();
        assert_eq!(base.positions, state.positions);
        assert_ne!(prior, &state.velocities);
        assert!(state.system().unwrap().total_momentum().norm() > 0.);
        prior = &state.velocities;
    }
    assert_eq!(
        (1..=6).map(|k| sample_step(0, k)).collect::<Vec<_>>(),
        [3, 5, 8, 10, 13, 15]
    );
    assert_eq!(
        validation::matched::roster()
            .iter()
            .map(|j| j.steps)
            .sum::<u64>(),
        297000
    );
    let mut bad = base.clone();
    bad.masses[1] = 0.;
    assert!(bad.system().is_err());
    let mut bad = base.clone();
    bad.positions.pop();
    assert!(bad.system().is_err());
}

#[test]
fn weighted_fit_retains_shared_seed_uncertainty() {
    use validation::protocol::Estimate;
    use validation::statistics::{fit, Inference};
    let dt = [0.04_f64, 0.02, 0.01, 0.005];
    let rows: Vec<_> = dt
        .iter()
        .map(|h| {
            let mu = 2. + 5000. * h * h;
            let means: Vec<_> = (0..8)
                .map(|seed| mu + (seed as f64 - 3.5) * 0.001)
                .collect();
            let se = validation::statistics::replica_se(&means).unwrap();
            Inference {
                estimate: Estimate {
                    mean: mu,
                    se,
                    ci95: [mu - se, mu + se],
                    n_eff: 1500.,
                    block_samples: 5,
                    retained_blocks: 160,
                },
                seed_means: means,
                bootstrap_ci95: None,
                tau_samples: 0.5,
                enough: true,
            }
        })
        .collect();
    let f = fit(&rows, &dt, 2.).unwrap();
    assert!((f.intercept - 2.).abs() < 1e-12 && (f.slope - 5000.).abs() < 1e-8);
    assert!((f.se - rows[0].estimate.se).abs() < 1e-13);
    assert!(f.weak_order_resolved && f.drop_paired_se < 1e-13);
    assert!(f.weighted_lack_of_fit < 1e-15);
    let mut sign_change = rows;
    for (r, bias) in sign_change.iter_mut().zip([16., 4., -1., -0.25]) {
        let shift = 2. + bias - r.estimate.mean;
        r.estimate.mean += shift;
        for seed_mean in &mut r.seed_means {
            *seed_mean += shift;
        }
    }
    let incompatible = fit(&sign_change, &dt, 2.).unwrap();
    assert!(
        !incompatible.weak_order_resolved,
        "a resolved sign reversal contradicts the other adjacent ratios"
    );
    let record = serde_json::to_value(&incompatible).unwrap();
    assert_eq!(
        record["resolved_bias_sign_changes"],
        serde_json::json!([[0.02, 0.01, 4., -1.]])
    );
}
#[test]
fn gate_publication_retains_partial_bytes_without_replacing_a_completed_file() {
    use serde::ser::{Error, SerializeSeq};
    struct Interrupted;
    impl serde::Serialize for Interrupted {
        fn serialize<S: serde::Serializer>(&self, s: S) -> Result<S::Ok, S::Error> {
            let mut seq = s.serialize_seq(None)?;
            seq.serialize_element(&17)?;
            Err(S::Error::custom("engineering interrupted serializer"))
        }
    }
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("ensemble-gates.json");
    assert!(validation::protocol::write_new(&path, &Interrupted).is_err());
    assert!(
        !path.exists(),
        "incomplete bytes must not be published as a completed gate file"
    );
    let partial = dir.path().join("ensemble-gates.json.partial");
    assert_eq!(std::fs::read(&partial).unwrap(), b"[17");
    let done = dir.path().join("water-400-gates.json");
    validation::protocol::write_new(&done, &[1, 2]).unwrap();
    assert!(validation::protocol::write_new(&done, &[3]).is_err());
    assert_eq!(std::fs::read(&done).unwrap(), b"[1,2]\n");
}
#[test]
fn finite_checks_precede_reductions() {
    for bad in [f64::NAN, f64::INFINITY, -f64::INFINITY] {
        assert!(validation::protocol::gate("observed", &[1., bad], &[1., 0.], &[1., 1.]).is_err());
        assert!(validation::protocol::gate("expected", &[1., 0.], &[1., bad], &[1., 1.]).is_err());
        assert!(validation::protocol::gate("tolerance", &[1., 0.], &[1., 0.], &[1., bad]).is_err());
    }
}

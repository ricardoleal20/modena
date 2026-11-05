use std::path::Path;

use modena_core::{
    ConservationApplicability, ConservationVerdict, EnergyErrorThresholds, EnergyMonitor,
    InputConfig, RunEvent, Simulation,
};
use serde_json::Value;

const BASE: &str = r#"
[structure]
type = "fcc"
lattice_constant = 2.4
replicas = [1,1,1]
[model]
type = "lennard_jones"
sigma = 1.0
epsilon = 1.0
cutoff = 2.5
[integrator]
type = "velocity_verlet"
dt = 0.001
[run]
steps = 16
initial_temperature = 1.0
seed = 51
log_every = 1
"#;

fn run(text: &str) -> Value {
    let config = InputConfig::from_str(text).unwrap();
    let summary = Simulation::from_config_without_outputs(&config, Path::new("."))
        .unwrap()
        .run()
        .unwrap();
    serde_json::to_value(summary.conservation).unwrap()
}

#[test]
fn runner_samples_every_step_in_reduced_time_independently_of_log_cadence() {
    let dense = run(BASE);
    let sparse = run(&BASE.replace("log_every = 1", "log_every = 100"));
    let d = &dense["diagnostics"];
    assert_eq!(
        d["samples"], 17,
        "initial state plus all 16 integration steps"
    );
    assert_eq!(dense["applicability"], "nve");
    assert_eq!(d["normalization"], "initial_kinetic_energy");
    assert!((d["reference_energy"].as_f64().unwrap() - 4.5).abs() < 1e-12);
    assert!((d["duration"].as_f64().unwrap() - 0.016).abs() < 1e-15);
    assert_eq!(d["energy_units"], "reduced");
    assert_eq!(d["time_units"], "reduced");
    assert_eq!(
        dense, sparse,
        "logging must not alter the numerical diagnostic"
    );
}

#[test]
fn thermostatted_and_pressure_coupled_runs_do_not_receive_an_nve_verdict() {
    let nvt = BASE.replace(
        "type = \"velocity_verlet\"",
        "type = \"langevin_baoab\"\ngamma = 1.0\ntarget_temperature = 1.0",
    );
    let npt = format!(
        "{nvt}\n[barostat]\ntype = \"monte_carlo\"\ntarget_pressure = 1.0\ntarget_temperature = 1.0\ninterval = 2\n"
    );
    for text in [nvt, npt] {
        let report = run(&text);
        assert_eq!(report["applicability"], "energy_exchange");
        assert_eq!(report["verdict"], "not_applicable");
        assert_eq!(report["diagnostics"]["samples"], 17);
        assert!(report["diagnostics"]["delta_energy"].is_number());
    }
}

#[test]
fn cold_nve_reports_a_missing_kinetic_reference_instead_of_inventing_a_scale() {
    let report = run(&BASE.replace("initial_temperature = 1.0", "initial_temperature = 0.0"));
    assert_eq!(report["applicability"], "nve");
    assert_eq!(report["verdict"], "unavailable");
    assert!(report["diagnostics"]["reference_energy"].is_null());
    assert!(report["diagnostics"]["normalized_max_abs_deviation"].is_null());
    assert!(report["diagnostics"]["energy_range"].is_number());
}

#[test]
fn normalized_metrics_are_uncalibrated_until_a_caller_selects_error_thresholds() {
    assert_eq!(run(BASE)["verdict"], "uncalibrated");
    let report = run(&format!(
        "{BASE}\n[run.energy_error_thresholds]\nwarning = 0.1\ndrifting = 0.2\n"
    ));
    assert_eq!(report["verdict"], "healthy");
    assert_eq!(report["diagnostics"]["advisory_thresholds"]["warning"], 0.1);
    for (warning, drifting) in [(0.0, 0.1), (-1.0, 0.1), (0.2, 0.1), (0.1, 0.1)] {
        let text = format!(
            "{BASE}\n[run.energy_error_thresholds]\nwarning = {warning}\ndrifting = {drifting}\n"
        );
        assert!(InputConfig::from_str(&text).is_err());
    }
}

#[test]
fn online_diagnostic_matches_an_analytic_linear_energy_series() {
    let mut monitor = EnergyMonitor::new(-10.0, 2, 2.0, ConservationApplicability::Nve).unwrap();
    for (time, energy) in [(1.0, -9.75), (2.0, -9.5), (3.0, -9.25)] {
        monitor.observe(time, energy).unwrap();
    }
    let report = monitor.report();
    let d = report.diagnostics.unwrap();
    assert_eq!(d.samples, 4);
    assert_eq!(d.delta_energy, 0.75);
    assert_eq!(d.delta_per_atom, 0.375);
    assert_eq!(d.energy_range, 0.75);
    assert_eq!(d.max_abs_deviation, 0.75);
    assert!((d.rms_fluctuation - 0.078125_f64.sqrt()).abs() < 1e-14);
    assert!((d.slope.unwrap() - 0.25).abs() < 1e-14);
    assert_eq!(d.normalized_change, Some(0.375));
    assert_eq!(d.normalized_max_abs_deviation, Some(0.375));
}

#[test]
fn adding_an_energy_constant_does_not_change_the_new_metrics_or_verdict() {
    let evaluate = |offset| {
        let mut monitor = EnergyMonitor::new(1.0 + offset, 2, 1.0, ConservationApplicability::Nve)
            .unwrap()
            .with_thresholds(EnergyErrorThresholds {
                warning: 1e-4,
                drifting: 1e-3,
            })
            .unwrap();
        for (t, delta) in [
            (0.25, 1.0 / 2048.0),
            (0.5, -1.0 / 4096.0),
            (0.75, 1.0 / 4096.0),
        ] {
            monitor.observe(t, 1.0 + offset + delta).unwrap();
        }
        monitor.report()
    };
    let base = evaluate(0.0);
    let shifted = evaluate(1024.0);
    assert_eq!(base.diagnostics, shifted.diagnostics);
    assert_eq!(base.verdict, shifted.verdict);
    assert_eq!(base.verdict, ConservationVerdict::Warning);
    // Preserve the old number for consumers that explicitly need it; it is not the verdict basis.
    assert_ne!(base.relative_drift, shifted.relative_drift);
}

#[test]
fn returning_to_initial_energy_does_not_hide_a_bounded_excursion() {
    let mut monitor = EnergyMonitor::new(1.0, 2, 1.0, ConservationApplicability::Nve)
        .unwrap()
        .with_thresholds(EnergyErrorThresholds {
            warning: 0.1,
            drifting: 0.2,
        })
        .unwrap();
    for (t, e) in [(1.0, 1.5), (2.0, 1.0), (3.0, 1.5), (4.0, 1.0)] {
        monitor.observe(t, e).unwrap();
    }
    let report = monitor.report();
    let d = report.diagnostics.unwrap();
    assert_eq!(d.delta_energy, 0.0);
    assert_eq!(d.energy_range, 0.5);
    assert_eq!(d.slope, Some(0.0));
    assert_eq!(d.normalized_max_abs_deviation, Some(0.5));
    assert_eq!(report.verdict, ConservationVerdict::Drifting);
}

#[test]
fn monitor_rejects_invalid_samples_without_corrupting_previous_statistics() {
    let mut monitor = EnergyMonitor::new(1.0, 2, 1.0, ConservationApplicability::Nve).unwrap();
    let initial = monitor.report();
    assert_eq!(
        serde_json::to_value(initial.verdict).unwrap(),
        "insufficient_data"
    );
    for (t, e) in [
        (0.0, 1.0),
        (-1.0, 1.0),
        (1.0, f64::NAN),
        (f64::INFINITY, 1.0),
    ] {
        assert!(monitor.observe(t, e).is_err());
        assert_eq!(monitor.report(), initial);
    }
    assert!(EnergyMonitor::new(0.0, 0, 1.0, ConservationApplicability::Nve).is_err());
    assert!(EnergyMonitor::new(f64::NAN, 2, 1.0, ConservationApplicability::Nve).is_err());
    assert!(EnergyMonitor::new(0.0, 2, -1.0, ConservationApplicability::Nve).is_err());
}

#[test]
fn thermo_events_carry_full_rate_statistics_even_when_output_is_sparse() {
    let config = InputConfig::from_str(&BASE.replace("log_every = 1", "log_every = 8")).unwrap();
    let mut sim = Simulation::from_config_without_outputs(&config, Path::new(".")).unwrap();
    let mut samples = Vec::new();
    let summary = sim
        .run_streaming(None, |event| {
            if let RunEvent::Thermo {
                step, conservation, ..
            } = event
            {
                samples.push((*step, *conservation));
            }
        })
        .unwrap();
    assert_eq!(
        samples.iter().map(|(step, _)| *step).collect::<Vec<_>>(),
        [0, 8, 16]
    );
    assert_eq!(samples[0].1.diagnostics.unwrap().samples, 1);
    assert_eq!(samples[1].1.diagnostics.unwrap().samples, 9);
    assert_eq!(samples[2].1, summary.conservation);
}

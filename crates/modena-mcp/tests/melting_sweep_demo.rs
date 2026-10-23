//! Deterministic integration test of the agent workflow — the paper's
//! reproducibility demo, run without a live LLM.
//!
//! Simulates an agent finding structure in a temperature sweep:
//!   build_system → estimate_run_cost → launch_sweep (T ladder) →
//!   poll get_sweep_results → detect_anomaly.
//!
//! Every step is a *typed tool call* over the shared `RunManager` core — no
//! shell, no script generation, no text parsing. This is the headline
//! capability of Feature A.

use modena_core::InputConfig;
use modena_mcp::{RunManager, RunState, SweepField};

/// A reduced LJ liquid config (FCC start, BAOAB-Langevin NVT). The agent
/// sweeps `target_temperature`.
fn base_config() -> InputConfig {
    let toml = r#"
units = "reduced"
[structure]
type = "fcc"
lattice_constant = 1.6796
replicas = [3, 3, 3]
element = "Ar"
[model]
type = "lennard_jones"
sigma = 1.0
epsilon = 1.0
cutoff = 2.5
[integrator]
type = "langevin_baoab"
dt = 0.005
gamma = 1.0
target_temperature = 1.0
seed = 7
[run]
steps = 500
initial_temperature = 1.0
seed = 42
"#;
    toml::from_str(toml).expect("valid config")
}

#[test]
fn agent_sweep_workflow_is_reproducible() {
    let mgr = RunManager::default();

    // 1. The agent first inspects the structure it intends to simulate.
    let built = mgr.build_system(&base_config().structure).unwrap();
    assert_eq!(built.atoms, 108);
    assert!((built.density - 0.8442).abs() < 1e-3);

    // 2. It budgets the run before launching.
    let cost = mgr.estimate_run_cost(&base_config()).unwrap();
    assert!(cost.est_wall_s > 0.0);
    assert!(cost.est_rss_mb > 0.0);

    // 3. It launches a temperature ladder as one typed call.
    let sweep = mgr
        .launch_sweep(
            base_config(),
            SweepField::LangevinTargetTemperature,
            vec![0.7, 0.9, 1.1, 1.3, 1.5],
        )
        .unwrap();
    assert_eq!(sweep.run_ids.len(), 5);

    // 4. It waits for every run to finish (a real agent polls get_run_status).
    for rid in &sweep.run_ids {
        mgr.join_run(*rid).unwrap();
    }

    // 5. It collects the aggregated sweep results.
    let results = mgr.get_sweep_results(sweep.sweep_id).unwrap();
    assert_eq!(results.runs.len(), 5);
    for row in &results.runs {
        assert_eq!(
            row.status.state,
            RunState::Completed,
            "run {} did not complete: {:?}",
            row.run_id,
            row.status
        );
        assert!(row.status.latest_thermo.is_some());
        // The sweep is Langevin NVT, so the conservation verdict is *not*
        // physically meaningful here (the thermostat exchanges energy). The
        // manager correctly marks these as non-NVE so drift anomalies are
        // suppressed for them.
        assert!(
            !row.status.is_nve,
            "Langevin sweep run {} should not be flagged NVE",
            row.run_id
        );
    }

    // 6. It asks the engine for structured anomaly findings.
    let report = mgr.detect_anomaly(None, Some(sweep.sweep_id)).unwrap();
    // The transition detector should surface exactly one possible_phase_transition
    // finding (the largest |d<PE>/dT| interval), since all 5 runs completed.
    let transitions: Vec<_> = report
        .anomalies
        .iter()
        .filter(|a| a.typ == "possible_phase_transition")
        .collect();
    assert_eq!(
        transitions.len(),
        1,
        "expected one transition finding, got {}: {:?}",
        transitions.len(),
        report.anomalies
    );
    // The finding carries a structured, agent-citable location + evidence.
    assert!(transitions[0]
        .location
        .starts_with("langevin_target_temperature"));
    assert!(transitions[0].evidence["interval"].is_array());
    assert!(transitions[0].evidence["n_completed"].as_u64() == Some(5));
}

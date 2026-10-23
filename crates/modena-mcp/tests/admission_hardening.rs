use std::io::Cursor;
use std::path::PathBuf;
use std::sync::{Arc, Barrier};
use std::thread;
use std::time::{Duration, Instant};

use modena_core::{InputConfig, StructureSpec};
use modena_mcp::{call_tool, server, McpError, RunManager, RunState, SweepField};
use serde_json::{json, Value};

fn config(steps: usize) -> InputConfig {
    serde_json::from_value(json!({
        "structure": {"type": "fcc", "lattice_constant": 2.4, "replicas": [2, 2, 2]},
        "model": {"type": "lennard_jones", "sigma": 1.0, "epsilon": 1.0, "cutoff": 2.5},
        "integrator": {"type": "velocity_verlet", "dt": 0.001},
        "run": {"steps": steps, "initial_temperature": 1.0, "log_every": 1000}
    }))
    .unwrap()
}

fn stop_runs(mgr: &RunManager, ids: impl IntoIterator<Item = u64>) {
    for id in ids {
        if mgr.get_run_status(id).is_ok() {
            mgr.cancel_run(id).unwrap();
            mgr.join_run(id).unwrap();
        }
    }
}

fn assert_bad_tool_input_keeps_server_alive(name: &str, arguments: Value) {
    let input = [
        json!({"jsonrpc": "2.0", "id": 1, "method": "tools/call",
            "params": {"name": name, "arguments": arguments}}),
        json!({"jsonrpc": "2.0", "id": 2, "method": "ping"}),
    ]
    .map(|v| v.to_string())
    .join("\n");
    let mgr = RunManager::default();
    let mut output = Vec::new();
    server::serve(&mgr, Cursor::new(input), &mut output).unwrap();
    let responses: Vec<Value> = String::from_utf8(output)
        .unwrap()
        .lines()
        .map(|line| serde_json::from_str(line).unwrap())
        .collect();
    assert_eq!(responses.len(), 2);
    assert!(responses[0]["result"]["content"][0]["text"]
        .as_str()
        .unwrap()
        .starts_with("invalid tool input:"));
    assert_eq!(responses[1]["result"], json!({}));
}

#[test]
fn fcc_count_overflow_is_a_recoverable_tool_error() {
    assert_bad_tool_input_keeps_server_alive(
        "build_system",
        json!({"structure": {
            "type": "fcc", "lattice_constant": 1.0, "replicas": [65536, 65536, 1]
        }}),
    );
}

#[test]
fn rocksalt_count_overflow_is_a_recoverable_tool_error() {
    assert_bad_tool_input_keeps_server_alive(
        "build_system",
        json!({"structure": {
            "type": "rock_salt", "lattice_constant": 1.0, "replicas": [65536, 65536, 1],
            "cation": "Na", "anion": "Cl", "cation_charge": 1.0, "anion_charge": -1.0
        }}),
    );
}

#[test]
fn lattice_boxes_must_have_finite_edges_and_volume() {
    for structure in [
        json!({"type": "fcc", "lattice_constant": 1e308, "replicas": [2, 1, 1]}),
        json!({"type": "fcc", "lattice_constant": 1e200, "replicas": [1, 1, 1]}),
        json!({"type": "rock_salt", "lattice_constant": 1e308, "replicas": [2, 1, 1],
            "cation": "Na", "anion": "Cl", "cation_charge": 1.0, "anion_charge": -1.0}),
    ] {
        let spec: StructureSpec = serde_json::from_value(structure).unwrap();
        assert!(RunManager::default().build_system(&spec).is_err());
    }
}

#[test]
fn concurrent_callers_share_one_atomic_admission_limit() {
    let mgr = Arc::new(RunManager::new(1, PathBuf::from(".")));
    let barrier = Arc::new(Barrier::new(12));
    let callers: Vec<_> = (0..12)
        .map(|_| {
            let mgr = Arc::clone(&mgr);
            let barrier = Arc::clone(&barrier);
            thread::spawn(move || {
                barrier.wait();
                mgr.run_async(config(100_000_000))
            })
        })
        .collect();
    let results: Vec<_> = callers.into_iter().map(|h| h.join().unwrap()).collect();
    let admitted: Vec<_> = results
        .iter()
        .filter_map(|r| r.as_ref().ok().copied())
        .collect();
    stop_runs(&mgr, admitted.iter().copied());
    assert_eq!(admitted.len(), 1, "admitted run ids: {admitted:?}");
    assert!(results
        .iter()
        .filter_map(|r| r.as_ref().err())
        .all(|e| matches!(e, McpError::PoolFull(1))));
}

#[test]
fn later_invalid_sweep_values_do_not_launch_earlier_runs() {
    let mgr = RunManager::new(2, PathBuf::from("."));
    let result = mgr.launch_sweep(
        config(100_000_000),
        SweepField::InitialTemperature,
        vec![0.8, -1.0],
    );
    let launched = mgr.get_run_status(1).is_ok();
    stop_runs(&mgr, 1..=2);
    assert!(result.is_err());
    assert!(!launched, "a rejected sweep started its first run");
}

#[test]
fn every_sweep_config_is_validated_before_launching() {
    let mgr = RunManager::new(2, PathBuf::from("."));
    let result = mgr.launch_sweep(config(100_000_000), SweepField::Density, vec![0.8, 1e-320]);
    let launched = mgr.get_run_status(1).is_ok();
    stop_runs(&mgr, 1..=2);
    assert!(
        result.is_err(),
        "a density producing an infinite box was accepted"
    );
    assert!(!launched, "a rejected sweep started its first run");
}

#[test]
fn insufficient_capacity_rejects_the_whole_sweep() {
    let mgr = RunManager::new(2, PathBuf::from("."));
    let existing = mgr.run_async(config(100_000_000)).unwrap();
    let deadline = Instant::now() + Duration::from_secs(10);
    while mgr.get_run_status(existing).unwrap().state == RunState::Queued {
        assert!(
            Instant::now() < deadline,
            "initial run did not finish setup"
        );
        thread::yield_now();
    }
    assert_eq!(
        mgr.get_run_status(existing).unwrap().state,
        RunState::Running
    );
    let result = mgr.launch_sweep(
        config(100_000_000),
        SweepField::InitialTemperature,
        vec![0.8, 1.2],
    );
    let extra_run = mgr.get_run_status(existing + 1).is_ok();
    stop_runs(&mgr, existing..=existing + 2);
    assert!(matches!(result, Err(McpError::PoolFull(2))));
    assert!(!extra_run, "a rejected sweep consumed a remaining slot");
}

#[test]
fn terminal_runs_release_capacity_for_the_next_admission() {
    let mgr = RunManager::new(1, PathBuf::from("."));
    let mut invalid = config(1);
    invalid.run.log_every = 0;
    let failed = mgr.run_async(invalid).unwrap();
    mgr.join_run(failed).unwrap();
    assert_eq!(mgr.get_run_status(failed).unwrap().state, RunState::Failed);
    for _ in 0..2 {
        let id = mgr.run_async(config(1)).unwrap();
        mgr.join_run(id).unwrap();
        assert_eq!(mgr.get_run_status(id).unwrap().state, RunState::Completed);
    }
}

#[test]
fn cost_estimate_count_overflow_is_a_recoverable_tool_error() {
    let mut input = serde_json::to_value(config(1000)).unwrap();
    input["structure"]["replicas"] = json!([u32::MAX, u32::MAX, u32::MAX]);
    assert_bad_tool_input_keeps_server_alive("estimate_run_cost", json!({"config": input}));
}

#[test]
fn cost_estimate_rejects_unsupported_counts_and_invalid_boxes() {
    for structure in [
        json!({"type": "fcc", "lattice_constant": 1.0, "replicas": [65536, 65536, 1]}),
        json!({"type": "fcc", "lattice_constant": 1.0, "replicas": [0, 1, 1]}),
        json!({"type": "fcc", "lattice_constant": 1e308, "replicas": [2, 1, 1]}),
        json!({"type": "rock_salt", "lattice_constant": 1.0, "replicas": [65536, 65536, 1],
            "cation": "Na", "anion": "Cl", "cation_charge": 1.0, "anion_charge": -1.0}),
    ] {
        let mut input = config(1000);
        input.structure = serde_json::from_value(structure).unwrap();
        let result = call_tool(
            "estimate_run_cost",
            &json!({"config": input}),
            &RunManager::default(),
        );
        assert!(
            result.is_err(),
            "invalid geometry received an estimate: {result:?}"
        );
    }
}

#[test]
fn cost_estimate_reports_unreadable_xyz_instead_of_zero_atoms() {
    let mut input = config(1000);
    input.structure = serde_json::from_value(json!({
        "type": "xyz", "box_size": [10.0, 10.0, 10.0],
        "file": std::env::temp_dir().join(format!("modena-missing-cost-input-{}.xyz", std::process::id()))
    }))
    .unwrap();
    let result = call_tool(
        "estimate_run_cost",
        &json!({"config": input}),
        &RunManager::default(),
    );
    assert!(
        result.is_err(),
        "an unreadable structure received an estimate: {result:?}"
    );
}

#[test]
fn cost_estimate_rejects_empty_xyz() {
    let path = std::env::temp_dir().join(format!(
        "modena-empty-cost-input-{}-{}.xyz",
        std::process::id(),
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos()
    ));
    std::fs::write(&path, "0\nempty frame\n").unwrap();
    let mut input = config(1000);
    input.structure = serde_json::from_value(json!({
        "type": "xyz", "box_size": [10.0, 10.0, 10.0], "file": path
    }))
    .unwrap();
    let result = call_tool(
        "estimate_run_cost",
        &json!({"config": input}),
        &RunManager::default(),
    );
    std::fs::remove_file(path).unwrap();
    assert!(
        result.is_err(),
        "an empty structure received an estimate: {result:?}"
    );
}

#[test]
fn cost_estimate_preserves_the_valid_budget_response() {
    let result = call_tool(
        "estimate_run_cost",
        &json!({"config": config(1000)}),
        &RunManager::default(),
    )
    .unwrap();
    assert_eq!(result["atoms"], 32);
    assert_eq!(result["steps"], 1000);
    assert!((result["est_wall_s"].as_f64().unwrap() - 0.0016).abs() < 1e-12);
    assert!((result["est_rss_mb"].as_f64().unwrap() - 5.00384).abs() < 1e-12);
}

use modena_mcp::{call_tool, RunManager};
use serde_json::{json, Value};
fn config() -> Value {
    json!({"structure":{"type":"fcc","lattice_constant":2.4,"replicas":[1,1,1]},"model":{"type":"lennard_jones","sigma":1.0,"epsilon":1.0,"cutoff":2.5},"integrator":{"type":"velocity_verlet","dt":0.001},"run":{"steps":8,"initial_temperature":1.0,"log_every":1}})
}
#[test]
fn barostat_run_is_not_nve_through_request_path() {
    let mgr = RunManager::default();
    let mut c = config();
    c["barostat"] =
        json!({"type":"berendsen","target_pressure":0.1,"tau":1.0,"compressibility":0.01});
    let r = call_tool("run_async", &json!({"config":c}), &mgr).unwrap();
    let id = r["run_id"].as_u64().unwrap();
    mgr.join_run(id).unwrap();
    let st = call_tool("get_run_status", &json!({"run_id":id}), &mgr).unwrap();
    assert_eq!(st["state"], "completed");
    assert_eq!(st["is_nve"], false);
}
#[test]
fn actual_means_counts_and_real_units_are_exposed() {
    let mgr = RunManager::default();
    let mut c = config();
    c["units"] = json!("real");
    c["structure"]["mass"] = json!(23.0);
    c["model"]["sigma"] = json!(2.0);
    c["model"]["epsilon"] = json!(0.5);
    let r = call_tool("run_async", &json!({"config":c}), &mgr).unwrap();
    let id = r["run_id"].as_u64().unwrap();
    mgr.join_run(id).unwrap();
    let st = call_tool("get_run_status", &json!({"run_id":id}), &mgr).unwrap();
    let events = call_tool("get_run_events", &json!({"run_id":id,"since_step":0}), &mgr).unwrap();
    let samples = events.as_array().unwrap();
    let expected = samples
        .iter()
        .map(|s| s["energies"]["potential"].as_f64().unwrap())
        .sum::<f64>()
        / samples.len() as f64;
    println!(
        "actual real event mean: expected={expected:.12}, reported={}, samples={}, window={}",
        st["mean_potential_energy"], st["sample_count"], st["sample_window"]
    );
    assert_eq!(st["thermo_units"], "real");
    assert_eq!(st["metadata"]["output_units"], "real");
    assert_eq!(st["sample_count"], samples.len());
    assert!((st["mean_potential_energy"].as_f64().unwrap() - expected).abs() < 1e-12);
    assert_eq!(st["sample_window"]["includes_initial"], false);
}
#[test]
fn discovery_and_probe_surface_supported_domain() {
    let mgr = RunManager::default();
    let cap = call_tool("get_capabilities", &json!({}), &mgr);
    assert!(cap.is_ok(), "missing capabilities {cap:?}");

    let schema = modena_mcp::tool_list()
        .into_iter()
        .find(|t| t.name == "run_probe")
        .unwrap()
        .input_schema
        .to_string();
    for field in [
        "lj_scales",
        "coulomb_scales",
        "exceptions",
        "water",
        "spme",
        "monte_carlo",
        "restart_mode",
        "model.types",
    ] {
        assert!(schema.contains(field), "missing {field}");
    }
    let mut c = config();
    c["structure"]["replicas"] = json!([2, 2, 2]);
    c["constraints"] = json!({"bonds":[{"atoms":[0,1],"length":1.697056274847714}]});
    let result = call_tool("run_probe", &json!({"config":c,"max_steps":2}), &mgr).unwrap();
    assert_eq!(
        result["metadata"]["pressure_status"],
        "diagnostic_without_constraint_impulse_stress"
    );
    assert_eq!(result["thermo_units"], "reduced");
}

#[test]
fn canonical_temperature_sweep_rejects_mc_coupling_explicitly() {
    let mgr = RunManager::default();
    let mut c = config();
    c["thermostat"] = json!({"type":"velocity_rescaling","target_temperature":1.0});
    c["barostat"] = json!({"type":"monte_carlo","target_pressure":0.1,"target_temperature":1.0});
    let r = call_tool(
        "launch_sweep",
        &json!({"base_config":c,"field":"thermostat_target_temperature","values":[0.8,1.2]}),
        &mgr,
    );
    assert!(
        r.is_err(),
        "MC coupled sweep must reject synchronously, received {r:?}"
    );
    assert!(r.unwrap_err().to_string().contains("MC"));
}
#[test]
fn accumulated_mean_survives_event_fifo_and_transition_uses_it() {
    let mgr = RunManager::default();
    let mut c = config();
    c["run"]["steps"] = json!(4200);
    c["run"]["log_every"] = json!(1);
    let sweep = call_tool(
        "launch_sweep",
        &json!({"base_config":c,"field":"initial_temperature","values":[0.7,1.1]}),
        &mgr,
    )
    .unwrap();
    for id in sweep["run_ids"].as_array().unwrap() {
        mgr.join_run(id.as_u64().unwrap()).unwrap();
    }
    let id = sweep["sweep_id"].as_u64().unwrap();
    let results = call_tool("get_sweep_results", &json!({"sweep_id":id}), &mgr).unwrap();
    for row in results["runs"].as_array().unwrap() {
        assert_eq!(row["status"]["sample_count"], 4200);
        let events = call_tool("get_run_events", &json!({"run_id":row["run_id"]}), &mgr).unwrap();
        assert_eq!(events.as_array().unwrap().len(), 4096);
    }
    println!("accumulated sweep rows: {}", results["runs"]);
    let mean0 = results["runs"][0]["status"]["mean_potential_energy"]
        .as_f64()
        .unwrap();
    let mean1 = results["runs"][1]["status"]["mean_potential_energy"]
        .as_f64()
        .unwrap();
    let report = call_tool("detect_anomaly", &json!({"sweep_id":id}), &mgr).unwrap();
    let transition = report["anomalies"]
        .as_array()
        .unwrap()
        .iter()
        .find(|a| a["type"] == "possible_phase_transition")
        .unwrap();
    assert!(
        (transition["evidence"]["slope_magnitude"].as_f64().unwrap() - (mean1 - mean0).abs() / 0.4)
            .abs()
            < 1e-12
    );
}

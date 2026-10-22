use modena_mcp::{call_tool, RunManager};
use serde_json::json;
#[test]
fn discovery_request_lists_every_production_pair_family_with_schema_parameters() {
    let mgr = RunManager::default();
    let result = call_tool("list_potentials", &json!({}), &mgr).unwrap();
    let entries = result.as_array().unwrap();
    assert_eq!(entries.len(), 5);
    for (name, params) in [
        (
            "lennard_jones",
            vec!["sigma", "epsilon", "cutoff", "types", "overrides", "mixing"],
        ),
        (
            "buckingham",
            vec!["a", "b", "rho", "cutoff", "types", "overrides"],
        ),
        (
            "morse",
            vec!["d", "alpha", "r0", "cutoff", "types", "overrides"],
        ),
        (
            "soft_sphere",
            vec!["sigma", "epsilon", "n", "cutoff", "types", "overrides"],
        ),
        ("tabulated", vec!["file", "interpolation"]),
    ] {
        let entry = entries
            .iter()
            .find(|e| e["name"] == name)
            .unwrap_or_else(|| panic!("missing {name}"));
        let mut actual: Vec<_> = entry["params"]
            .as_array()
            .unwrap()
            .iter()
            .map(|p| p.as_str().unwrap())
            .collect();
        let mut expected = params;
        actual.sort_unstable();
        expected.sort_unstable();
        assert_eq!(actual, expected);
    }
}
#[test]
fn exact_tiny_alloy_completes_real_mcp_probe_and_async_lifecycle() {
    let mgr = RunManager::default();
    let c: modena_core::InputConfig = toml::from_str(include_str!(
        "../../modena-core/tests/fixtures/final_alloy.toml"
    ))
    .unwrap();
    let config = serde_json::to_value(c).unwrap();
    let built = call_tool(
        "build_system",
        &json!({"structure":config["structure"]}),
        &mgr,
    )
    .unwrap();
    assert_eq!(built["atoms"], 4);
    let probe = call_tool("run_probe", &json!({"config":config,"max_steps":2}), &mgr).unwrap();
    assert_eq!(probe["steps"], 2);
    assert_eq!(probe["metadata"]["kinetic_dof"], 9);
    let run = call_tool("run_async", &json!({"config":config}), &mgr).unwrap();
    let id = run["run_id"].as_u64().unwrap();
    mgr.join_run(id).unwrap();
    let status = call_tool("get_run_status", &json!({"run_id":id}), &mgr).unwrap();
    assert_eq!(status["state"], "completed");
    assert_eq!(status["step"], 2);
    let mut invalid = config.clone();
    invalid["structure"]["species"][0]["concentration"] = json!(2.0);
    let failed = call_tool("run_async", &json!({"config":invalid}), &mgr).unwrap();
    let failed_id = failed["run_id"].as_u64().unwrap();
    mgr.join_run(failed_id).unwrap();
    let status = call_tool("get_run_status", &json!({"run_id":failed_id}), &mgr).unwrap();
    assert_eq!(status["state"], "failed");
    assert!(status["error"].as_str().unwrap().contains("concentrations"));
    assert!(call_tool("build_system",&json!({"structure":{"type":"random_alloy","lattice_constant":3.0,"replicas":[4294967295u32,4294967295u32,4294967295u32],"species":[{"element":"X","mass":1.0,"concentration":1.0}]}}),&mgr).is_err());
}

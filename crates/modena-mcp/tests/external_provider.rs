use modena_mcp::{call_tool, tool_list, RunManager};
use serde_json::{json, Value};

// Exercise the actual newline-delimited JSON-RPC process, including MCP error wrapping.
fn rpc_tool(name: &str, args: Value) -> Value {
    use std::io::Write;
    use std::process::{Command, Stdio};
    let mut child = Command::new(env!("CARGO_BIN_EXE_modena-mcp"))
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .unwrap();
    let request = json!({"jsonrpc":"2.0", "id":1, "method":"tools/call",
        "params":{"name":name,"arguments":args}});
    writeln!(child.stdin.take().unwrap(), "{request}").unwrap();
    let output = child.wait_with_output().unwrap();
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let response: Value = serde_json::from_slice(&output.stdout).unwrap();
    assert_eq!(response["id"], 1);
    response["result"].clone()
}

fn config(python: &str, version: &str, hash: &str) -> Value {
    json!({"units":"real",
        "structure":{"type":"fcc","lattice_constant":3.61,"replicas":[1,1,1],"mass":63.546},
        "model":{"type":"ase_emt","python":python,"ase_version":version,"model_sha256":hash,
                 "atomic_numbers":[29,29,29,29],"composition":"total","timeout_ms":120000},
        "integrator":{"type":"velocity_verlet","dt":0.0001},
        "run":{"steps":2,"initial_temperature":10.0,"log_every":100}})
}

#[test]
fn external_tool_estimate_is_explicitly_unsupported_without_launching_python() {
    let mgr = RunManager::default();
    let cfg = config("/nonexistent/modena-python", "3.26.0", &"a".repeat(64));
    let error = call_tool("estimate_run_cost", &json!({"config":cfg}), &mgr)
        .unwrap_err()
        .to_string();
    assert!(
        error.contains("unsupported") && error.contains("ase_emt"),
        "{error}"
    );
    let result = rpc_tool("estimate_run_cost", json!({"config":cfg}));
    assert_eq!(result["isError"], true);
    assert!(result["content"][0]["text"]
        .as_str()
        .unwrap()
        .contains("unsupported cost estimate for ase_emt"));
}

#[test]
fn discovery_schema_exposes_exact_external_pins_and_scope() {
    let tools = tool_list();
    let schema = &tools
        .iter()
        .find(|t| t.name == "run_probe")
        .unwrap()
        .input_schema;
    let models = schema["definitions"]["config"]["properties"]["model"]["oneOf"]
        .as_array()
        .unwrap();
    let model = models
        .iter()
        .find(|m| m["properties"]["type"]["const"] == "ase_emt")
        .expect("external model must be discoverable through actual tool schema");
    assert_eq!(model["additionalProperties"], false);
    assert_eq!(
        model["required"],
        json!([
            "type",
            "ase_version",
            "model_sha256",
            "atomic_numbers",
            "composition"
        ])
    );
    let props = model["properties"].as_object().unwrap();
    let mut fields: Vec<_> = props.keys().map(String::as_str).collect();
    fields.sort_unstable();
    assert_eq!(
        fields,
        vec![
            "asap_cutoff",
            "ase_version",
            "atomic_numbers",
            "composition",
            "model_sha256",
            "pbc",
            "python",
            "timeout_ms",
            "type"
        ]
    );
    assert_eq!(props["pbc"]["minItems"], 3);
    assert_eq!(props["pbc"]["maxItems"], 3);
    assert_eq!(props["pbc"]["items"]["const"], true);
    assert_eq!(props["asap_cutoff"]["const"], false);
    assert_eq!(props["composition"]["const"], "total");
    assert_eq!(
        props["atomic_numbers"]["items"]["enum"],
        json!([13, 28, 29, 46, 47, 78, 79])
    );
    let caps = call_tool("get_capabilities", &json!({}), &RunManager::default()).unwrap();
    let external = &caps["external_providers"][0];
    assert_eq!(external["name"], "ase_emt");
    assert_eq!(external["optional"], true);
    assert_eq!(external["cost_estimate"], "unsupported");
    assert_eq!(external["neighbors"], "provider_managed");
    assert_eq!(external["material_validation"], false);
    assert!(external["cancellation"]
        .as_str()
        .unwrap()
        .contains("timeout"));
}

#[test]
fn actual_ase_tool_json_probe_and_async_status_preserve_discovered_identity() {
    let Ok(python) = std::env::var("MODENA_ASE_PYTHON") else {
        eprintln!("SKIP: set MODENA_ASE_PYTHON for actual ASE MCP reference");
        return;
    };
    let identity = modena_core::describe_ase_emt(&python, 120000).unwrap();
    let cfg = config(&python, &identity.ase_version, &identity.model_sha256);
    let wire = rpc_tool("run_probe", json!({"config":cfg,"max_steps":2}));
    assert!(wire.get("isError").is_none(), "{wire}");
    let wire_probe: Value =
        serde_json::from_str(wire["content"][0]["text"].as_str().unwrap()).unwrap();
    assert_eq!(
        wire_probe["metadata"]["provider"],
        serde_json::to_value(&identity).unwrap()
    );
    let mgr = RunManager::default();
    let probe = call_tool("run_probe", &json!({"config":cfg,"max_steps":2}), &mgr).unwrap();
    assert_eq!(
        probe["metadata"]["provider"],
        serde_json::to_value(&identity).unwrap()
    );
    assert_eq!(
        probe["metadata"]["pressure_status"],
        "unconstrained_provider_stress"
    );
    let run = call_tool("run_async", &json!({"config":cfg}), &mgr).unwrap();
    let id = run["run_id"].as_u64().unwrap();
    // Observe populated live metadata before the final state; provider startup can be slow.
    let start = std::time::Instant::now();
    loop {
        let status = call_tool("get_run_status", &json!({"run_id":id}), &mgr).unwrap();
        if !status["metadata"].is_null() {
            assert_eq!(
                status["metadata"]["provider"],
                serde_json::to_value(&identity).unwrap()
            );
            break;
        }
        assert_ne!(status["state"], "failed", "{status}");
        assert!(start.elapsed().as_secs() < 120, "{status}");
        std::thread::sleep(std::time::Duration::from_millis(10));
    }
    mgr.join_run(id).unwrap();
    let status = call_tool("get_run_status", &json!({"run_id":id}), &mgr).unwrap();
    assert_eq!(status["state"], "completed");
    assert_eq!(
        status["metadata"]["provider"],
        serde_json::to_value(identity).unwrap()
    );
}

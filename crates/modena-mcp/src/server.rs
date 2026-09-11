//! The MCP stdio server: reads newline-delimited JSON-RPC from a reader,
//! dispatches to [`RunManager`](crate::RunManager) tools, and writes responses
//! to a writer. Generic over `BufRead` + `Write` so it is unit-testable with
//! `Cursor`.

use std::io::{BufRead, Write};

use crate::protocol::{
    CallToolResult, InitializeResult, Response, RpcError, ServerInfo, ToolDescriptor,
    PROTOCOL_VERSION,
};
use crate::tools::{call_tool, tool_list};
use crate::RunManager;

/// Run the server until the input stream is exhausted (EOF). Each line is one
/// JSON-RPC message. Responses (for messages carrying an id) are written one
/// per line to `output`, which is flushed after every write.
pub fn serve<R: BufRead, W: Write>(
    mgr: &RunManager,
    input: R,
    mut output: W,
) -> std::io::Result<()> {
    for line in input.lines() {
        let line = line?;
        if line.trim().is_empty() {
            continue;
        }
        let response = handle_line(mgr, &line);
        if let Some(resp) = response {
            serde_json::to_writer(&mut output, &resp)?;
            writeln!(output)?;
            output.flush()?;
        }
    }
    Ok(())
}

/// Handle one raw JSON line. Returns `None` for notifications (no id).
fn handle_line(mgr: &RunManager, line: &str) -> Option<Response> {
    let req: crate::protocol::Request = match serde_json::from_str(line) {
        Ok(r) => r,
        Err(e) => {
            // Parse error: respond with null id per JSON-RPC.
            return Some(Response::err(
                serde_json::Value::Null,
                RpcError::new(RpcError::PARSE_ERROR, format!("parse error: {e}")),
            ));
        }
    };

    // Notifications (no id) get no response.
    let id = match req.id.clone() {
        Some(id) if !id.is_null() => id,
        _ => {
            // Still handle side-effectful notifications (none currently), then drop.
            return None;
        }
    };

    let result = dispatch(mgr, &req.method, &req.params);
    match result {
        Ok(value) => Some(Response::ok(id, value)),
        Err(rpc_err) => Some(Response::err(id, rpc_err)),
    }
}

fn dispatch(
    mgr: &RunManager,
    method: &str,
    params: &serde_json::Value,
) -> Result<serde_json::Value, RpcError> {
    match method {
        "initialize" => Ok(serde_json::to_value(InitializeResult {
            protocol_version: PROTOCOL_VERSION,
            capabilities: crate::protocol::Capabilities {
                tools: crate::protocol::ToolsCapability {},
            },
            server_info: ServerInfo {
                name: "modena-mcp",
                version: env!("CARGO_PKG_VERSION"),
            },
        })
        .unwrap()),
        "ping" => Ok(serde_json::json!({})),
        "tools/list" => {
            let tools: Vec<ToolDescriptor> = tool_list()
                .into_iter()
                .map(|t| ToolDescriptor {
                    name: t.name,
                    description: t.description,
                    input_schema: t.input_schema,
                })
                .collect();
            Ok(serde_json::json!({ "tools": tools }))
        }
        "tools/call" => {
            #[derive(serde::Deserialize)]
            struct CallParams {
                name: String,
                #[serde(default)]
                arguments: serde_json::Value,
            }
            let params: CallParams = serde_json::from_value(params.clone())
                .map_err(|e| RpcError::new(RpcError::INVALID_PARAMS, e.to_string()))?;
            match call_tool(&params.name, &params.arguments, mgr) {
                Ok(value) => Ok(serde_json::to_value(CallToolResult::text(value)).unwrap()),
                Err(e) => {
                    // Tool execution error → MCP error result (not an RPC error).
                    Ok(serde_json::to_value(CallToolResult::error(e.to_string())).unwrap())
                }
            }
        }
        other => Err(RpcError::new(
            RpcError::METHOD_NOT_FOUND,
            format!("method not found: {other}"),
        )),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    fn sample_config_json() -> serde_json::Value {
        serde_json::json!({
            "units": "reduced",
            "structure": {
                "type": "fcc",
                "lattice_constant": 1.6796,
                "replicas": [2, 2, 2],
                "element": "Ar"
            },
            "model": {
                "type": "lennard_jones",
                "sigma": 1.0,
                "epsilon": 1.0,
                "cutoff": 2.5
            },
            "integrator": { "type": "velocity_verlet", "dt": 0.002 },
            "run": { "steps": 100, "initial_temperature": 1.0, "seed": 42 }
        })
    }

    /// Send a sequence of JSON-RPC lines, return the collected response lines.
    fn run_session(lines: Vec<String>) -> Vec<String> {
        let mgr = RunManager::default();
        let input = lines.join("\n");
        let mut output = Vec::new();
        serve(&mgr, Cursor::new(input), &mut output).unwrap();
        String::from_utf8(output)
            .unwrap()
            .lines()
            .filter(|l| !l.is_empty())
            .map(String::from)
            .collect()
    }

    fn rpc(method: &str, id: i64, params: serde_json::Value) -> String {
        serde_json::json!({ "jsonrpc": "2.0", "id": id, "method": method, "params": params })
            .to_string()
    }

    #[test]
    fn initialize_and_tools_list() {
        let responses = run_session(vec![
            rpc("initialize", 1, serde_json::json!({})),
            rpc("tools/list", 2, serde_json::json!({})),
        ]);
        assert_eq!(responses.len(), 2);
        let init: serde_json::Value = serde_json::from_str(&responses[0]).unwrap();
        assert_eq!(init["result"]["protocol_version"], PROTOCOL_VERSION);
        assert_eq!(init["result"]["server_info"]["name"], "modena-mcp");
        let list: serde_json::Value = serde_json::from_str(&responses[1]).unwrap();
        let names: Vec<&str> = list["result"]["tools"]
            .as_array()
            .unwrap()
            .iter()
            .map(|t| t["name"].as_str().unwrap())
            .collect();
        assert!(names.contains(&"build_system"));
        assert!(names.contains(&"launch_sweep"));
        assert!(names.contains(&"detect_anomaly"));
    }

    #[test]
    fn build_system_via_mcp() {
        let responses = run_session(vec![rpc(
            "tools/call",
            1,
            serde_json::json!({
                "name": "build_system",
                "arguments": { "structure": sample_config_json()["structure"] }
            }),
        )]);
        assert_eq!(responses.len(), 1);
        let resp: serde_json::Value = serde_json::from_str(&responses[0]).unwrap();
        let text = &resp["result"]["content"][0]["text"];
        let payload: serde_json::Value = serde_json::from_str(text.as_str().unwrap()).unwrap();
        assert_eq!(payload["atoms"], 32);
        assert!((payload["density"].as_f64().unwrap() - 0.8442).abs() < 1e-3);
    }

    #[test]
    fn run_probe_via_mcp() {
        let responses = run_session(vec![rpc(
            "tools/call",
            1,
            serde_json::json!({
                "name": "run_probe",
                "arguments": { "config": sample_config_json(), "max_steps": 100 }
            }),
        )]);
        assert_eq!(responses.len(), 1);
        let resp: serde_json::Value = serde_json::from_str(&responses[0]).unwrap();
        // Success omits `is_error` (it is only present when true).
        assert!(
            resp["result"]["is_error"].is_null(),
            "unexpected error result: {resp}"
        );
        let text = resp["result"]["content"][0]["text"].as_str().unwrap();
        let payload: serde_json::Value = serde_json::from_str(text).unwrap();
        assert_eq!(payload["steps"], 100);
    }

    #[test]
    fn unknown_method_returns_rpc_error() {
        let responses = run_session(vec![rpc("nope", 1, serde_json::json!({}))]);
        let resp: serde_json::Value = serde_json::from_str(&responses[0]).unwrap();
        assert_eq!(resp["error"]["code"], RpcError::METHOD_NOT_FOUND);
    }

    #[test]
    fn notifications_get_no_response() {
        // A notification (no id) should produce no output line.
        let mgr = RunManager::default();
        let input = r#"{"jsonrpc":"2.0","method":"notifications/initialized"}"#;
        let mut output = Vec::new();
        serve(&mgr, Cursor::new(input), &mut output).unwrap();
        assert!(
            output.is_empty(),
            "notification produced output: {:?}",
            output
        );
    }

    #[test]
    fn invalid_tool_args_is_error_result() {
        let responses = run_session(vec![rpc(
            "tools/call",
            1,
            serde_json::json!({ "name": "build_system", "arguments": {} }),
        )]);
        let resp: serde_json::Value = serde_json::from_str(&responses[0]).unwrap();
        // Missing required field → reported as an error tool result (is_error true).
        assert_eq!(resp["result"]["is_error"], true);
    }
}

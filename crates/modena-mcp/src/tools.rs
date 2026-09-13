//! Tool descriptors and dispatch — the mapping from MCP `tools/call` to
//! [`RunManager`](crate::RunManager) methods.
//!
//! Each tool declares a JSON-Schema `input_schema` (so the agent's tool-use is
//! type-checked by the client) and a handler that validates input against the
//! existing `InputConfig` types before any run starts. This is the "safe, typed
//! primitives" claim: contrast with blind LAMMPS script generation.

use crate::{McpError, RunManager, SweepField};
use modena_core::{InputConfig, StructureSpec};

/// A tool descriptor: static metadata exposed via `tools/list`.
pub struct Tool {
    pub name: &'static str,
    pub description: &'static str,
    pub input_schema: serde_json::Value,
}

/// The full v1 tool set, in a stable order.
pub fn tool_list() -> Vec<Tool> {
    vec![
        Tool {
            name: "list_potentials",
            description: "List the registered pair potentials and their parameter names.",
            input_schema: schema!({
                "type": "object",
                "properties": {},
                "additionalProperties": false
            }),
        },
        Tool {
            name: "build_system",
            description: "Validate and materialize a structure (FCC lattice or XYZ file) without running. Returns atom count, number density, and box size. Use this to confirm a structure before launching a run.",
            input_schema: schema!({
                "type": "object",
                "properties": {
                    "structure": { "$ref": "#/definitions/structure" }
                },
                "required": ["structure"],
                "additionalProperties": false,
                "definitions": structure_definitions()
            }),
        },
        Tool {
            name: "estimate_run_cost",
            description: "Estimate wall-clock seconds and RSS megabytes for a config from an internal cost model. Helps budget before large runs. Figures are rough estimates, not guarantees.",
            input_schema: schema!({
                "type": "object",
                "properties": {
                    "config": { "$ref": "#/definitions/config" }
                },
                "required": ["config"],
                "additionalProperties": false,
                "definitions": config_definitions()
            }),
        },
        Tool {
            name: "run_probe",
            description: "Run a short synchronous simulation (steps clamped to max_steps, default 500) and return the final thermodynamic state plus the energy-conservation report. Use for quick parameter checks.",
            input_schema: schema!({
                "type": "object",
                "properties": {
                    "config": { "$ref": "#/definitions/config" },
                    "max_steps": { "type": "integer", "minimum": 1, "default": 500 }
                },
                "required": ["config"],
                "additionalProperties": false,
                "definitions": config_definitions()
            }),
        },
        Tool {
            name: "run_async",
            description: "Launch a streaming simulation on a background thread. Returns a run_id immediately. Poll with get_run_status; cancel with cancel_run.",
            input_schema: schema!({
                "type": "object",
                "properties": {
                    "config": { "$ref": "#/definitions/config" }
                },
                "required": ["config"],
                "additionalProperties": false,
                "definitions": config_definitions()
            }),
        },
        Tool {
            name: "get_run_status",
            description: "Snapshot a run's state, current step, latest thermo, and conservation report.",
            input_schema: schema!({
                "type": "object",
                "properties": {
                    "run_id": { "type": "integer", "minimum": 0 }
                },
                "required": ["run_id"],
                "additionalProperties": false
            }),
        },
        Tool {
            name: "get_run_events",
            description: "Return buffered thermodynamic samples for a run with step > since_step.",
            input_schema: schema!({
                "type": "object",
                "properties": {
                    "run_id": { "type": "integer", "minimum": 0 },
                    "since_step": { "type": "integer", "minimum": 0, "default": 0 }
                },
                "required": ["run_id"],
                "additionalProperties": false
            }),
        },
        Tool {
            name: "cancel_run",
            description: "Request a clean stop for a run. It transitions to 'cancelled' within one integration step.",
            input_schema: schema!({
                "type": "object",
                "properties": {
                    "run_id": { "type": "integer", "minimum": 0 }
                },
                "required": ["run_id"],
                "additionalProperties": false
            }),
        },
        Tool {
            name: "launch_sweep",
            description: "Launch one async run per value of a swept parameter, all from a base config. Returns a sweep_id and the per-value run_ids. The pool cap limits the number of concurrent runs.",
            input_schema: schema!({
                "type": "object",
                "properties": {
                    "base_config": { "$ref": "#/definitions/config" },
                    "field": {
                        "type": "string",
                        "enum": ["initial_temperature", "langevin_target_temperature",
                                 "thermostat_target_temperature", "density"]
                    },
                    "values": { "type": "array", "items": { "type": "number" }, "minItems": 1 }
                },
                "required": ["base_config", "field", "values"],
                "additionalProperties": false,
                "definitions": config_definitions()
            }),
        },
        Tool {
            name: "get_sweep_results",
            description: "Aggregate the status of every run in a sweep, each annotated with its swept value.",
            input_schema: schema!({
                "type": "object",
                "properties": {
                    "sweep_id": { "type": "integer", "minimum": 0 }
                },
                "required": ["sweep_id"],
                "additionalProperties": false
            }),
        },
        Tool {
            name: "detect_anomaly",
            description: "Run engine-side anomaly heuristics over a run and/or a sweep: energy drift (verdict 'drifting'), divergence (NaN/Inf), and basic phase-transition detection (largest |d<PE>/d(param)| across a sweep). Returns structured findings.",
            input_schema: schema!({
                "type": "object",
                "properties": {
                    "run_id": { "type": "integer", "minimum": 0 },
                    "sweep_id": { "type": "integer", "minimum": 0 }
                },
                "additionalProperties": false
            }),
        },
    ]
}

/// Dispatch a `tools/call` to the right `RunManager` method. Returns JSON
/// (success) or an [`McpError`] (reported to the agent as an error result).
pub fn call_tool(
    name: &str,
    args: &serde_json::Value,
    mgr: &RunManager,
) -> Result<serde_json::Value, McpError> {
    match name {
        "list_potentials" => Ok(serde_json::to_value(mgr.list_potentials()).unwrap()),
        "build_system" => {
            #[derive(serde::Deserialize)]
            struct In {
                structure: StructureSpec,
            }
            let input: In = parse(args)?;
            let res = mgr.build_system(&input.structure)?;
            Ok(serde_json::to_value(res).unwrap())
        }
        "estimate_run_cost" => {
            #[derive(serde::Deserialize)]
            struct In {
                config: InputConfig,
            }
            let input: In = parse(args)?;
            let res = mgr.estimate_run_cost(&input.config);
            Ok(serde_json::to_value(res).unwrap())
        }
        "run_probe" => {
            #[derive(serde::Deserialize)]
            struct In {
                config: InputConfig,
                #[serde(default = "default_probe_steps")]
                max_steps: usize,
            }
            let input: In = parse(args)?;
            let res = mgr.run_probe(input.config, input.max_steps)?;
            Ok(serde_json::to_value(res).unwrap())
        }
        "run_async" => {
            #[derive(serde::Deserialize)]
            struct In {
                config: InputConfig,
            }
            let input: In = parse(args)?;
            let run_id = mgr.run_async(input.config)?;
            Ok(serde_json::json!({ "run_id": run_id }))
        }
        "get_run_status" => {
            #[derive(serde::Deserialize)]
            struct In {
                run_id: u64,
            }
            let input: In = parse(args)?;
            let res = mgr.get_run_status(input.run_id)?;
            Ok(serde_json::to_value(res).unwrap())
        }
        "get_run_events" => {
            #[derive(serde::Deserialize)]
            struct In {
                run_id: u64,
                #[serde(default)]
                since_step: usize,
            }
            let input: In = parse(args)?;
            let res = mgr.get_run_events(input.run_id, input.since_step)?;
            Ok(serde_json::to_value(res).unwrap())
        }
        "cancel_run" => {
            #[derive(serde::Deserialize)]
            struct In {
                run_id: u64,
            }
            let input: In = parse(args)?;
            let res = mgr.cancel_run(input.run_id)?;
            Ok(serde_json::to_value(res).unwrap())
        }
        "launch_sweep" => {
            #[derive(serde::Deserialize)]
            struct In {
                base_config: InputConfig,
                field: SweepField,
                values: Vec<f64>,
            }
            let input: In = parse(args)?;
            let res = mgr.launch_sweep(input.base_config, input.field, input.values)?;
            Ok(serde_json::to_value(res).unwrap())
        }
        "get_sweep_results" => {
            #[derive(serde::Deserialize)]
            struct In {
                sweep_id: u64,
            }
            let input: In = parse(args)?;
            let res = mgr.get_sweep_results(input.sweep_id)?;
            Ok(serde_json::to_value(res).unwrap())
        }
        "detect_anomaly" => {
            #[derive(serde::Deserialize, Default)]
            struct In {
                run_id: Option<u64>,
                sweep_id: Option<u64>,
            }
            let input: In = parse(args).unwrap_or_default();
            let res = mgr.detect_anomaly(input.run_id, input.sweep_id)?;
            Ok(serde_json::to_value(res).unwrap())
        }
        other => Err(McpError::Other(format!("unknown tool: {other}"))),
    }
}

fn default_probe_steps() -> usize {
    500
}

fn parse<T: serde::de::DeserializeOwned>(args: &serde_json::Value) -> Result<T, McpError> {
    serde_json::from_value(args.clone()).map_err(|e| McpError::InvalidInput(e.to_string()))
}

// ---- JSON-Schema fragments ----
//
// These are deliberately loose (no exhaustive `required` on every nested field)
// so the agent isn't over-constrained; the real validation happens when we
// deserialize into `InputConfig`, which rejects malformed configs with a
// precise error.

fn structure_definitions() -> serde_json::Value {
    json!({
        "structure": {
            "oneOf": [
                {
                    "type": "object",
                    "properties": {
                        "type": { "const": "fcc" },
                        "lattice_constant": { "type": "number" },
                        "replicas": { "type": "array", "items": { "type": "integer" }, "minItems": 3, "maxItems": 3 },
                        "element": { "type": "string" }
                    },
                    "required": ["type", "lattice_constant", "replicas"]
                },
                {
                    "type": "object",
                    "properties": {
                        "type": { "const": "xyz" },
                        "file": { "type": "string" },
                        "box_size": { "type": "array", "items": { "type": "number" }, "minItems": 3, "maxItems": 3 },
                        "element": { "type": "string" }
                    },
                    "required": ["type", "file", "box_size"]
                },
                {
                    "type": "object",
                    "properties": {
                        "type": { "const": "rock_salt" },
                        "lattice_constant": { "type": "number" },
                        "replicas": { "type": "array", "items": { "type": "integer" }, "minItems": 3, "maxItems": 3 },
                        "cation": { "type": "string" },
                        "anion": { "type": "string" },
                        "cation_charge": { "type": "number" },
                        "anion_charge": { "type": "number" }
                    },
                    "required": ["type", "lattice_constant", "replicas", "cation_charge", "anion_charge"]
                }
            ]
        }
    })
}

fn config_definitions() -> serde_json::Value {
    json!({
        "config": {
            "type": "object",
            "properties": {
                "units": { "type": "string", "default": "reduced" },
                "structure": { "$ref": "#/definitions/structure/oneOf/0" },
                "model": {
                    "type": "object",
                    "properties": {
                        "type": { "const": "lennard_jones" },
                        "sigma": { "type": "number" },
                        "epsilon": { "type": "number" },
                        "cutoff": { "type": "number" }
                    },
                    "required": ["type", "sigma", "epsilon", "cutoff"]
                },
                "integrator": {
                    "type": "object",
                    "properties": {
                        "type": {
                            "enum": ["velocity_verlet", "langevin_baoab", "nose_hoover_chain"]
                        },
                        "dt": { "type": "number" },
                        "gamma": { "type": "number" },
                        "target_temperature": { "type": "number" },
                        "seed": { "type": "integer" }
                    },
                    "required": ["type", "dt"]
                },
                "thermostat": {
                    "type": "object",
                    "properties": {
                        "type": { "enum": ["berendsen", "velocity_rescaling", "none"] },
                        "target_temperature": { "type": "number" },
                        "tau": { "type": "number" },
                        "seed": { "type": "integer" }
                    }
                },
                "electrostatics": {
                    "type": "object",
                    "properties": {
                        "kappa": { "type": "number" },
                        "real_cutoff": { "type": "number" },
                        "kspace_cutoff": { "type": "number" }
                    },
                    "required": ["kappa", "real_cutoff", "kspace_cutoff"]
                },
                "run": {
                    "type": "object",
                    "properties": {
                        "steps": { "type": "integer", "minimum": 1 },
                        "initial_temperature": { "type": "number" },
                        "seed": { "type": "integer" },
                        "trajectory_file": { "type": "string" },
                        "trajectory_every": { "type": "integer" },
                        "log_file": { "type": "string" },
                        "log_every": { "type": "integer" }
                    },
                    "required": ["steps"]
                }
            },
            "required": ["structure", "model", "integrator", "run"]
        },
        "structure": structure_definitions().get("structure").cloned().unwrap_or(serde_json::Value::Null)
    })
}

/// `json!` with `additionalProperties: false` baked into a `schema!` macro is
/// verbose; this thin macro just forwards to `serde_json::json!`.
macro_rules! schema {
    ($($tt:tt)*) => { serde_json::json!($($tt)*) };
}
pub(crate) use schema;

/// Re-export `json!` for schema fragments.
macro_rules! json {
    ($($tt:tt)*) => { serde_json::json!($($tt)*) };
}
pub(crate) use json;

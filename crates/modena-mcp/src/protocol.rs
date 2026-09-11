//! Minimal JSON-RPC 2.0 + MCP protocol types for the stdio transport.
//!
//! The MCP "stdio" transport is newline-delimited JSON: each line on stdin is a
//! JSON-RPC request or notification, and each response is a JSON-RPC object on
//! stdout (one per line). This is hand-rolled rather than depending on an
//! external MCP SDK to keep the crate dependency-light and the protocol surface
//! fully under our control — the core stays sync + std.

use serde::{Deserialize, Serialize};

/// MCP protocol version advertised at `initialize`.
pub const PROTOCOL_VERSION: &str = "2024-11-05";

/// A JSON-RPC 2.0 id (number, string, or null). Kept as a raw JSON value so we
/// echo it back verbatim.
pub type Id = serde_json::Value;

/// An incoming JSON-RPC message (request or notification).
#[derive(Debug, Deserialize)]
pub struct Request {
    #[serde(default)]
    pub jsonrpc: Option<String>,
    /// `None`/null => notification (no response expected).
    #[serde(default)]
    pub id: Option<Id>,
    pub method: String,
    #[serde(default)]
    pub params: serde_json::Value,
}

/// A JSON-RPC 2.0 error object.
#[derive(Debug, Serialize)]
pub struct RpcError {
    pub code: i32,
    pub message: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub data: Option<serde_json::Value>,
}

impl RpcError {
    pub const PARSE_ERROR: i32 = -32700;
    pub const INVALID_REQUEST: i32 = -32600;
    pub const METHOD_NOT_FOUND: i32 = -32601;
    pub const INVALID_PARAMS: i32 = -32602;
    pub const INTERNAL_ERROR: i32 = -32603;

    pub fn new(code: i32, message: impl Into<String>) -> Self {
        Self {
            code,
            message: message.into(),
            data: None,
        }
    }
}

/// A JSON-RPC 2.0 response (only for requests that carried an id).
#[derive(Debug, Serialize)]
pub struct Response {
    pub jsonrpc: &'static str,
    pub id: Id,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub result: Option<serde_json::Value>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub error: Option<RpcError>,
}

impl Response {
    pub fn ok(id: Id, result: serde_json::Value) -> Self {
        Self {
            jsonrpc: "2.0",
            id,
            result: Some(result),
            error: None,
        }
    }

    pub fn err(id: Id, error: RpcError) -> Self {
        Self {
            jsonrpc: "2.0",
            id,
            result: None,
            error: Some(error),
        }
    }
}

/// Result of `tools/list`: a tool descriptor with a JSON-Schema input shape.
#[derive(Debug, Serialize, Clone)]
pub struct ToolDescriptor {
    pub name: &'static str,
    pub description: &'static str,
    pub input_schema: serde_json::Value,
}

/// `initialize` result payload.
#[derive(Debug, Serialize)]
pub struct InitializeResult {
    pub protocol_version: &'static str,
    pub capabilities: Capabilities,
    pub server_info: ServerInfo,
}

#[derive(Debug, Serialize)]
pub struct Capabilities {
    pub tools: ToolsCapability,
}

#[derive(Debug, Serialize)]
pub struct ToolsCapability {}

#[derive(Debug, Serialize)]
pub struct ServerInfo {
    pub name: &'static str,
    pub version: &'static str,
}

/// `tools/call` result payload, shaped per the MCP spec.
#[derive(Debug, Serialize)]
pub struct CallToolResult {
    pub content: Vec<Content>,
    #[serde(skip_serializing_if = "std::ops::Not::not")]
    pub is_error: bool,
}

#[derive(Debug, Serialize)]
pub struct Content {
    #[serde(rename = "type")]
    pub kind: &'static str,
    pub text: String,
}

impl CallToolResult {
    /// A successful result carrying JSON text.
    pub fn text(json: serde_json::Value) -> Self {
        Self {
            content: vec![Content {
                kind: "text",
                text: serde_json::to_string_pretty(&json).unwrap_or_else(|_| json.to_string()),
            }],
            is_error: false,
        }
    }

    /// An error result carrying a message.
    pub fn error(message: impl Into<String>) -> Self {
        Self {
            content: vec![Content {
                kind: "text",
                text: message.into(),
            }],
            is_error: true,
        }
    }
}

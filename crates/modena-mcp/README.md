# modena-mcp

A local MCP interface to **MODENA**, a molecular dynamics engine in Rust.
The server uses newline-delimited JSON-RPC over standard input and output and
advertises MCP protocol version `2024-11-05`. Scientific tool payloads use
MODENA's documented configuration fields.

From a [source checkout](https://github.com/ricardoleal20/modena), with Rust 1.90
or newer:

```sh
cargo install --locked --path crates/modena-mcp
```

Configure an MCP client to launch `modena-mcp` as a local stdio process. Use
`tools/list` to discover the supported tools and their input schemas. This
process runs with your filesystem permissions and is not a sandbox; use
operating-system isolation for untrusted jobs. Cost estimates are rough
heuristics, not resource guarantees.

MODENA is in early development. Read the
[validation scope](https://github.com/ricardoleal20/modena/blob/main/docs/VALIDATED_SCOPE.md)
and [security policy](https://github.com/ricardoleal20/modena/blob/main/SECURITY.md).

Licensed under MIT. The complete notice is included in `LICENSE`.

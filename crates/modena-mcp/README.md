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

Live conservation reports contain statistics from every integration step in
reduced energy/time, regardless of log cadence. NVE is `uncalibrated` until
`run.energy_error_thresholds` are selected; NVT/NPT are `not_applicable`.
Anomaly findings retain failures that occur between logged samples. Ordinary
sample means have no burn-in or autocorrelation correction. See the
[B06 diagnostic contract](https://github.com/ricardoleal20/modena/blob/main/docs/B06_ENERGY_DIAGNOSTICS.md).

The model schema and `list_potentials` include optional `ase_emt` alongside
native models. `get_capabilities` declares its real-unit fixed-cell unconstrained
NVE/NVT total-model scope, explicit standard-metal Z/masses, provider-managed
neighbors and one-shot Python startup latency. The local Python runtime needs
`ase>=3.26,<4`; native workflows do not require it. Discover installed pins with
`modena provider-info --python /path/to/venv/bin/python`, then supply the exact
`ase_version`, `model_sha256`, per-atom `atomic_numbers` and
`composition="total"`. Only `pbc=[true,true,true]` and `asap_cutoff=false` are
supported. See the [complete Python example](../../python/README.md#optional-ase-emt).

`run_probe` and `run_async` accept the external config without a cost estimate.
`estimate_run_cost` returns an explicit unsupported-estimate tool error for
`ase_emt`; the native pair-kernel timing/RSS model cannot budget Python/EMT.
Probe metadata and live/final `get_run_status.metadata.provider` carry the
validated core identity, including installed ASE/NumPy versions, model hash,
protocol, capabilities and units. Provider pressure is labeled
`unconstrained_provider_stress`. Cancellation is bounded by the configured
per-evaluation timeout (`timeout_ms`, default 30000) plus normal step completion.
The interface makes no material-validation or throughput claim.

MODENA is in early development. Read the
[validation scope](https://github.com/ricardoleal20/modena/blob/main/docs/VALIDATED_SCOPE.md)
and [security policy](https://github.com/ricardoleal20/modena/blob/main/SECURITY.md).

Licensed under MIT. The complete notice is included in `LICENSE`.

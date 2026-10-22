# Runner handoff risks

CLI check currently constructs Simulation and opens configured outputs; protect existing sentinel files in a real CLI regression. Full physical configuration validation must still occur with writers disabled.

CSVR at exactly zero kinetic energy cannot heat by global rescaling (current apply returns). Reject positive-temperature CSVR with zero-temperature initialization unless a physically correct constrained Maxwell draw is implemented; do not silently keep a frozen ideal gas. Langevin can activate stochastic COM modes, unlike a global rescaler; COM DOF differs.

MC logV criterion independently verified: https://dl_monte.gitlab.io/dl_monte-tutorials-pages/tutorial2-ext1.html. For unconstrained ideal gas p(V) is Gamma(N+1, kT/P), with mean (N+1) kT/P and variance (N+1)(kT/P)^2. State rollback and accepted caches need independent checks, not only correct volume statistics.

Do not label ordinary total-energy drift a conservation fault under NVT/NPT in CLI presentation. Its numerical value can be reported as energy exchanged; MCP NVE classifier gets fixed in Task 4.

Raw programmatic/JSON InputConfig currently bypasses from_str's real-unit conversion: MCP directly calls Simulation::from_config and the builder only validates. Add one idempotent normalization boundary: real config with no unit_factors is converted once, already-normalized config with unit_factors is not converted again. Test real-unit TOML vs direct/JSON configs for identical system/model/dt. Keep physical signature on normalized values. Task4 Python writes original physical TOML and MCP consumes this common boundary; real-unit thermo metadata must state actual output units.

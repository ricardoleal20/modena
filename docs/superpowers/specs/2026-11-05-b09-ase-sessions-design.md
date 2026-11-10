# B09: persistent EMT sessions and MODENA as an ASE calculator

The owner requested B09 followed immediately by B02 after merging B06/B01.
This implements both directions named in B09, each within an explicit subset.

## Persistent external provider

Retain the B01 one-shot path as the compatibility/default transport. Add an
explicit `transport = "persistent"` option to the pinned ASE EMT model, with
matching Python and MCP configuration. Transport is not physical identity and
must not invalidate otherwise compatible checkpoints. The existing provider
identity, physical units and fixed-cell NVE/NVT restrictions remain binding.

One provider instance owns one Python process/calculator session. Discovery and
all evaluations share it. Use a versioned framed protocol with monotonically
identified requests, a validated identity/unit handshake and bounded messages.
The per-request deadline covers writes and reads. Bound stderr and avoid pipe
deadlocks. Drop, timeout, EOF, malformed/stale replies, identity changes and
provider errors kill/reap the process and leave the session terminal; no silent
restart or stale property reuse. Validate snapshots before sending and commit
force/virial additions only after complete validation. Preserve cancellation at
the existing safe setup/step boundaries; an in-flight evaluation is bounded by
its deadline, without claiming immediate preemption.

## MODENA calculator for ASE

Provide optional `modena.ase.ModenaCalculator`, keeping base `import modena`
free of ASE/NumPy requirements. Start with single-species native shifted-force
Lennard-Jones, full periodic orthorhombic cells, no constraints, charges or spin.
Parameters use explicit Angstrom/eV units. A dedicated strict JSON CLI snapshot
command computes energy and forces without integration, velocities or files.
Use the existing core LJ force path and validate radius/cell/shape/finiteness.
The calculator advertises only energy and forces: scalar virial cannot supply
the full ASE stress tensor. Requesting stress must fail explicitly. ASE cache
invalidation must account for positions, numbers, cell/PBC and calculator
parameters, and failures must clear old results. Changing executable/parameters
must not leave an old result cached. Unsupported inputs fail before evaluation.

Snapshot cache comparison is exact, including changes below ASE's default
1e-15 tolerance; no-calculation cache probes preserve ASE's absent-result
semantics. Normalize finite periodic coordinates with remainder arithmetic,
including large positive/negative coordinates. A documented candidate/pair
storage budget must reject dense inputs before neighbor allocation grows past
the budget; atom/grid limits alone do not bound pair storage. Keep normal
unbounded MD neighbor behavior unchanged, and test the bounded path with a
small explicit budget rather than a destructive huge allocation.

## Verification and performance

Use direct installed ASE EMT references for changed positions/cells, assert one
process/calculator startup per persistent session, and exercise restart,
cancellation, deadlines, bad IDs/identity/properties, EOF and terminal failure.
For the reverse calculator use an independent analytical LJ pair/finite
difference oracle, and exercise cache invalidation and absent properties.
Measure repeated cold one-shot versus persistent full MODENA runs on the same
Cu configuration, hardware, versions, steps and outputs. Record all timings,
provenance, medians and matching numerical output. Separate startup and steady
evaluation where useful; do not generalize to DFT/ML solver cost or other engines.
The benchmark validates required physical checkpoint signature and thermo
columns before comparing runs, so shared omissions cannot satisfy parity.

## Global constraints

- Branches begin `ricardo/`; preserve unrelated files and the old B01 evidence.
- Rust floor 1.90; no new mandatory Python/ASE dependency or downloaded model.
- Shell commands begin `rtk`; commits use literal gitmoji plus imperative Action.
- No merge of B09/B02 is authorized by the earlier B06/B01-specific exception.
- Finish and review B09 before implementation of B02; publish reviewable PRs.

ASE implementation contract: https://docs.ase-lib.org/_modules/ase/calculators/calculator.html

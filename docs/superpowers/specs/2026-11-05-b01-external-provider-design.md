# B01: first external force provider

B01 uses ASE EMT and periodic copper as the first reproducible external-provider
reference. It connects a fixed, pinned total model to the ordinary runner,
checkpoints and interfaces within the bounded scope below.

## Scope and decision

Implement one concrete total-energy provider, `ase_emt`, using the official Python `ase.calculators.emt.EMT`. A fresh Python calculation receives each complete snapshot and returns energy, forces and configurational stress. ASE is an optional local dependency; the standard classical build retains no Python runtime requirement. A one-shot process is deliberately the initial parity mechanism. A persistent worker and general ASE calculator selection belong to B09; FAIR-Chem weights and DFT solvers require separate model/task validation.

Alternatives considered: a persistent generic calculator adds lifecycle/protocol scope before a first provider is proven; a snapshot-only library adapter would not exercise the actual TOML runner and restart identity. The selected vertical connects a fixed calculator to the ordinary runner and validates it against direct official calls. It makes no throughput claim.

ASE documents EMT as suitable for demonstrations/tests and supports Cu among its standard metals. The vertical restricts elements to those seven standard metals, not EMT's demonstration H/C/N/O parameterizations. Sources: [official EMT scope](https://docs.ase-lib.org/ase/calculators/emt.html), [official implementation/properties](https://docs.ase-lib.org/_modules/ase/calculators/emt.html), [ASE units](https://docs.ase-lib.org/ase/units.html), consulted 2026-11-05. The implementation invokes the installed dependency; it does not copy EMT equations or parameters.

## Identity and configuration

`ModelSpec::AseEmt(AseEmtSpec)` uses explicit per-atom `atomic_numbers`, distinct from force-field `species` and cosmetic element labels. `System.atomic_numbers` uses zero for unspecified classical identity; external calls require a supported nonzero number per atom. The configured vector is authoritative for the provider. Masses are explicit mechanical inputs, not guessed from Z.

Required fields: `ase_version`, `model_sha256`, `atomic_numbers`, and `composition = "total"`. Optional fields: `python = "python3"`, `pbc = [true,true,true]`, `asap_cutoff = false`, `timeout_ms = 30000`. Reject partial/nonperiodic PBC, non-real units, empty/unsupported identities, invalid pins, charged inputs, topology/electrostatic additions, constraints and pressure coupling for this first vertical. The initial supported run is unconstrained fixed-cell NVE/NVT. Reject residual composition rather than accidentally adding classical terms to the total model.

A `modena provider-info --python PATH` command discovers the protocol, installed ASE/NumPy versions, EMT source SHA-256, capabilities, task, parameters, composition and units. Model parameters are embedded in EMT's source, so that source hash identifies this provider's model data; it is not a learned-weights hash. Version/hash mismatches fail before accepting results. Runtime identity also records NumPy and fixed bridge protocol version, and becomes part of the physical checkpoint signature. Python executable location is transport configuration, not model identity.

## Evaluation and units

The bundled bridge is a core-package Python resource, executed without a shell. Request/response protocol `modena-ase-emt-v1` has explicit `info` and `evaluate` actions. Evaluate transfers atomic numbers, Cartesian positions, the three cell vectors as rows (ASE convention) and PBC per axis; every invocation creates a new Atoms and calculator. No provider result cache survives a call.

Wire properties use Angstrom, eV, eV/Angstrom and eV/Angstrom^3. The core remains reduced; with the external real-unit reference sigma=1 Angstrom and epsilon=1 kcal/mol, convert E and F using the exact-SI eV-to-kcal/mol factor `1.602176634e-19 * 6.02214076e23 / 4184`. For a general library-supplied scale, positions/cell multiply by sigma; reduced forces multiply eV/Angstrom by sigma*conversion/epsilon, and energy multiplies by conversion/epsilon. Obtain scalar virial from the stress derivative, `W = -V trace(stress)`, with the declared tensile-positive ASE convention. Do not invent stress or zero virial when a property is missing. This does not implement B07's stress/constraint-impulse work.

The provider manages its own periodic neighbors; do not advertise an invented cutoff or reuse MODENA's classical half-list as an ASE graph. The runner supplies an unused compatibility list only on the explicitly declared provider-managed path. Known EMT translation invariance may opt into momentum conservation; tests must check it. Report pressure metadata as provider-stress-derived, distinct from a pair virial.

Validate response identity, property presence, exact N-by-3 shape, stress shape/symmetry and finite converted values before adding any force or virial contribution. Failures leave prior accumulator contents intact. Use temporary regular files for stdin/stdout/stderr to avoid pipe deadlocks; enforce timeout and an output-size bound, kill and reap failed children, and preserve a bounded diagnostic error. No shell execution, silent retry, missing-property substitution or stale cache fallback is allowed. This transport is not an OS sandbox. Cancellation is observed before setup, between discovery/evaluation calls and before every integration step. One in-flight provider call may finish under its configured timeout plus ordinary bookkeeping; no hard global wall-time or OS-sandbox guarantee is made.

## Checkpoints and compatibility

Persist atomic numbers separately from force-field types. Older v2 checkpoints lacking that optional field remain usable for unspecified classical Z; external runs require matching Z. Extend only external signatures with PBC, provider identity, model hash and task; preserve native signatures so existing classical checkpoints remain compatible. Continuation must reject changed Z, model pin, version, task or cutoff convention. State-only import remains an intentional new physical protocol, but does not bypass provider validation. Evaluate current identity before restoring cached mechanics.

Keep classic pair construction/composition unchanged in its branch. Expose config fields through Python and MCP; MCP estimates for an external calculator must explicitly be unavailable rather than extrapolated from pair-kernel timings. The default Rust build must compile without ASE, and protocol failure tests must run without it. Optional ASE tests must distinguish absence/skipping from an actual physical reference run.

## Acceptance

1. Direct official EMT and the Rust provider agree on a displaced 4-Cu periodic snapshot for energy, all force components and stress-derived scalar virial, including a skew cell, within fixed absolute tolerances of 1e-10 eV, 1e-10 eV/Angstrom and 1e-9 eV respectively.
2. Central differences at h=1e-5 Angstrom match a selected force within 2e-6 eV/Angstrom; fixed-fractional isotropic strain at h=1e-6 matches W within 2e-6 eV. Common translation preserves energy/forces and force sum within 1e-10 in wire units.
3. Native configuration/regression tests remain passing. A real-unit external TOML executes a short NVE run with B06 full-rate samples and correct provider metadata. Do not equate this short run with ensemble/material validation.
4. Changed numbers, version/hash or cutoff convention reject continuation; valid split continuation reproduces uninterrupted mechanics. Missing properties, nonfinite/incorrect shape output, timeout and child failure fail explicitly without partially changing force accumulators.
5. CLI discovery, Python rendering and MCP config transport expose the same pins and scope. Documentation gives installation, reproducible commands, limitations, executed versions and evidence.

## Accepted preservation exception after whole-branch review

The shared real-pressure conversion is corrected from 69477.8 to the exact-SI
68568.42296625092 atm per kcal/(mol Angstrom³), derived from
4184/(6.02214076e23*1e-30*101325). This explicitly supersedes native pressure
preservation: real pressure output, targets and inverse-pressure compressibility
change; reduced configurations and real fixed-cell integration mechanics do not.
Old real fixed-cell v2 signatures may continue only when replacing the exact old
pressure factor at index5 of seven output factors makes the entire signature
exactly equal. No other factor/provider/Z/physical identity is relaxed. New
checkpoints write the corrected factor. Old real pressure-coupled protocols must
be reinitialized intentionally with state_only. Independent bidirectional SI,
MC/Berendsen normalization and direct ASE cold-Cu CLI tests gate this exception.
Failed Simulation instances are terminal and cannot retry or write checkpoints;
completed progress excludes a failed attempted step. Relative path-like Python
executables resolve lexically from config base_dir without dereferencing symlinks.

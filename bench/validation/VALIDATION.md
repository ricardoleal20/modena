# Scientific validation

Current acceptance uses numerical observations and independent expectations in
[scientific-hardening evidence](../../docs/scientific-hardening/RESULTS.md).
The prebuilt example runs bounded multi-seed CSVR Gamma kinetic moments,
constrained/unconstrained Langevin Maxwell moments, rigid triangle rotational
Gamma moments, force finite differences, NVE timestep convergence, SPME mesh
accuracy against direct Ewald and full-cell strain virials, and stochastic
JSON state continuation. This is regression validation, not converged liquid
thermodynamics or general chemical completeness.

```sh
cargo build --release -p modena-core --example scientific_validation
python3 bench/validate.py --output target/scientific-smoke
# Optional longer sampling (400000 samples per seed; smoke uses 40000):
python3 bench/validate.py --extended --output target/scientific-extended
```

Every numerical row records observed, expected, absolute tolerance and an
independent basis. Correlated sample moments are not formal uncertainty estimates.
The rigid rotation gate uses a conservative OU correlation band. Nonfinite or
out-of-band observations make the runner exit 1. Raw stdout/stderr, binary hash,
commit/dirty status and Python/compiler provenance are retained. CI runs the
bounded smoke and uploads it even on failure.

Earlier Buckingham/Ewald/SHAKE measurements and the original validation narrative
are archived under docs/scientific-hardening/evidence/legacy-docs and original-audit.
Historical 3N temperature, stale reciprocal virials and rebuild-every-step claims
are superseded. A projection removes the constrained fraction of KE in expectation,
not exactly for an individual draw. LAMMPS shifted-force LJ is energy-consistent.
Historical real-unit ionic energies/timings at sigma!=1 used the wrong charge
conversion and cannot validate the corrected Hamiltonian. Fixed-volume molecular
pressure omits unvalidated constraint impulse stress; rigid/anisotropic NPT is
unsupported. Existing decks/harness remain useful historical fixtures, not a
claim of full ensemble parity.

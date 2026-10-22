# Chemical workflow examples

Run from the repository root with the debug or release CLI:

```sh
cargo build -p modena-cli
./target/debug/modena check -i examples/chemical/molecular_nvt.toml
./target/debug/modena run -i examples/chemical/molecular_nvt.toml
./target/debug/modena check -i examples/chemical/continue_nvt.toml
./target/debug/modena run -i examples/chemical/continue_nvt.toml
./target/debug/modena run -i examples/chemical/ionic_spme.toml
./target/debug/modena run -i examples/chemical/mc_npt.toml
./target/debug/modena run -i examples/chemical/water_nvt.toml
```

Continuation requires the molecular run first. Files resolve relative to the
configuration file; outputs go to `examples/chemical/results/`. These short
examples demonstrate accepted workflows and are not equilibrium production
studies or validated material parameter sets. Water is fixed-volume and its
pressure is diagnostic. MC NPT requires unconstrained isotropic geometry and
matching positive canonical controller/barostat temperatures. Berendsen is
available for equilibration only.

Real inputs use Angstrom, ps, amu, elementary charge, kcal/mol, K and atm.
Reduced coordinates use r*=r/sigma, U*=U/epsilon; charge conversion is
q*=q_e sqrt(332.06371/(sigma_Angstrom epsilon_kcal/mol)). The corrected charge
conversion changes previous real electrostatic runs with sigma != 1; prior
results/checkpoints are not interchangeable through normal continuation.
Internal checkpoint mechanics and conservation energies remain reduced.
CSV, trajectory and streaming thermo values use selected output units.
Tables remain explicitly supplied in reduced units.

Typed XYZ requires the exact whitespace-delimited comment token
`Properties=species:S:1:pos:R:3:type:I:1:mass:R:1:charge:R:1`. Each following
row contains `element x y z zero_based_model_type mass charge`. Labels are
presentation; numeric types identify model entries. Unknown model indices,
nonfinite coordinates/charges, nonpositive masses and malformed declarations
fail. Positions/masses/charges convert once from the chosen input units.
The real typed reference mass is 1 amu. Global XYZ mass/charge fields are
allowed only for plain legacy XYZ; plain files use type0/mass1/charge0 unless
those explicit global fields override them.

Graph weights `topology.lj_scales` and `topology.coulomb_scales` correspond to
shortest bond distances 1,2,3 (1-2,1-3,1-4), independent of explicit angles or
dihedrals. Missing arrays are all1 or all0 when `exclude_bonded=true`.
`[[topology.exceptions]]` with `atoms`, `lj`, `coulomb` overrides graph weights;
conflicting duplicate unordered exceptions fail. Legacy angle/dihedral endpoint
exclusions remain when exclude_bonded is true and weighting arrays are absent.
Rigid water always excludes
its three intramolecular pairs; incompatible explicit water exceptions fail. Special Coulomb weights apply to the full
bare pair via Ewald/SPME corrections, including outside the real cutoff.

Python supports the same common Rust validator without dependencies:

```python
from modena import Config, Xyz, LennardJones, NoseHooverChain, Run, check, run
cfg = Config(Xyz("examples/chemical/chain.xyz", [20.0]*3),
             LennardJones(8.0, types=[{"sigma":3.0,"epsilon":0.2},
                                    {"sigma":1.0,"epsilon":0.05}]),
             NoseHooverChain(0.0001,100.0),units="real",
             run=Run(steps=5,initial_temperature=100.0,log_every=1))
# Add the explicit Topology policy for bonded production molecules.
print(check(cfg))
result = run(cfg)
print(result.steps_rows[-1]["temperature"])  # float, K
```

Python transport resolves source/table/restore/output paths from `base_dir`
(default current directory) before writing temporary TOML, preserving intended
restart/trajectory destinations. `run` captures CSV in a temporary log and
returns numeric rows (integer step; floating observables). Both subprocess
operations propagate nonzero exit status with stdout/stderr.

MCP `get_capabilities` reports scientific and execution limitations;
`tools/list` includes all production input sections. CPU/Rayon is the config
backend. GPU/MPI are optional core library backends, with no production input
selection or MCP runtime probe. Status/events/probe report selected thermo
units and constrained-pressure metadata; conservation energy is reduced.
Means/counts accumulate finite post-initial emitted samples in this invocation,
include equilibration, remove no burn-in, and outlive the 4096-event FIFO.
Transition detection uses these actual means and remains a finite-difference
heuristic. MC-coupled target-temperature sweeps reject explicitly; initial
thermal draw sweeps do not change canonical targets. NHC target sweeps are
not implemented.

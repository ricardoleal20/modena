## Task 2: Correct force/energy/pressure and remove measured allocation costs

**Ownership:** force_model.rs (remaining implementation), neighbor.rs, potential/tabulated.rs, ewald.rs, spme.rs, composite.rs, decomp.rs, gpu_comm.rs if contract changes; force/pair/electrostatic tests and core lib exports as needed. Do not edit Task 1 files except integration callsites strictly necessary and report these.

- [ ] Reproduce tabulated energy gradient and direct Ewald isotropic strain virial failures. Joint Hermite interpolation uses force-derived endpoint slopes; shifted-force cutoff remains consistent. Include tabulated in existing derivative properties.
- [ ] Correct reciprocal direct Ewald strain virial and net-charge background virial (W_bg=3 U_bg), with independent finite strain derivative tests. Check SPME background and weighted exclusion corrections.
- [ ] Replace linear exclusion scans with an indexed representation. Add `PairScale { i: usize, j: usize, lj: f64, coulomb: f64 }` and NeighborList `with_pair_scales(cutoff, skin, scales)` plus `pair_scale(i,j)`; retain existing constructors as binary 0/0 exclusion compatibility. Do not drop a pair when only one interaction is zero. Apply independent LJ/Coulomb scales, including reciprocal erf corrections proportional to (1-coulomb scale). Electrostatic constructors `with_pair_scales` match this policy.
- [ ] Implement neighbor_cutoff for pair, electrostatic, composite and backend models. Triclinic neighbor search must fall back to guaranteed correct pair enumeration if bounding-box cell stencil is unsafe; test the audited skew pair.
- [ ] Fix DomainForceModel accumulation contract with an existing nonzero-force fixture. Optional GPU/MPI code must remain compileable; unsupported fractional scaling backends reject it explicitly.
- [ ] Cache SPME FFT plans, grids and scratch per instance using a thread-safe owned workspace; reuse buffers and invalidate influence data by box/parameter key. Eliminate global cross-simulation cache thrashing. Cache pair force scratch if profiling/reproducible microbench supports it, retaining numerical tolerances and deterministic serial path.
- [ ] Run covering force properties, Ewald/SPME references, exclusion policy tests; provide allocation/reuse or timing evidence. Commit and report RED/GREEN, interfaces and restrictions.


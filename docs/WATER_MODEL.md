# Fixed-volume rigid-water model: tested scope and limits

The engine provides SPC/E geometry/charges, intramolecular electrostatic
exclusions, rigid constraints, orthorhombic SPME and checkpoint continuation.
Focused tests cover free-rotor conservation and canonical rotational moments,
geometry/velocity constraints and restart seams. Molecular pressure remains
**diagnostic without validated constraint impulse stress**. Constrained NPT and
full anisotropic coupling are rejected.

The Task4 correction is load-bearing: real charges use
`q*=q*sqrt(332.06371/(sigma*epsilon))`, because `r*=r/sigma` and `U*=U/epsilon`.
The former `sqrt(C*sigma/epsilon)` made physical Coulomb energy wrong by sigma².
Old real electrostatic runs at sigma!=1 therefore have a different Hamiltonian;
physical checkpoint signatures reject normal continuation. Reduced/sigma1 cases
are unchanged. No historical checkpoint migration is promised.

Earlier water timings, energies and unit-factor claims are historical and are
preserved in evidence/legacy-docs/WATER_MODEL.md. In particular the earlier
claim that the ~1.5 kcal/mol/water energy discrepancy could not be a unit bug
was unsupported. Neither structural relaxation nor a correct restart seam
explains away NVE energy drift. No converged bulk SPC/E liquid/LAMMPS parity or
pressure validation is claimed; equilibration, RDF, diffusion and matched static
energy/force checks for the corrected real Hamiltonian remain competition gaps.
See [RESULTS](scientific-hardening/RESULTS.md) for actual numerical coverage.

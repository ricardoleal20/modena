# Real-unit normalization and continuation

Real inputs use Angstrom, ps, kcal/mol, amu, K and elementary charge.
Internal lengths/energies use `r*=r/sigma`, `U*=U/epsilon`; therefore the
Coulomb identity requires `q*=q*sqrt(332.06371/(sigma*epsilon))`.
Then `epsilon*q_i* q_j*/r* = 332.06371*q_i*q_j/r` in kcal/mol.
Real thermo metadata labels output energies kcal/mol, temperature K and
pressure atm; it must not be interpreted as reduced-unit output.

Task4 fixed the earlier factor `sqrt(C*sigma/epsilon)`, whose physical energy
had a sigma² error. At sigma!=1 this changes the real electrostatic Hamiltonian.
Normal restart continuation checks the normalized physical identity including
charge factor and atom charges, and rejects mismatches. Reduced/sigma1 cases
are unchanged. Fresh reinitialization is a separate mode, not evidence that an
old Hamiltonian has been migrated. Historical real ionic benchmarks cannot be
compared as corrected-model energies or strict before/after performance.

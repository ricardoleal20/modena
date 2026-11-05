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


Standard-atmosphere pressure uses the exact-SI factor
`4184/(6.02214076e23*1e-30*101325) = 68568.42296625092` atm per kcal/(mol Å³).
The former 69477.8 scale was about 1.32623% high. This correction also changes
real pressure-target normalization and inverse-pressure compressibility; it
does not change energy/force/time/temperature conversions, reduced inputs or
fixed-cell real integration mechanics. Exact fixed-cell v2 migration accepts
only the sole historical pressure factor in the seven-entry output factors,
with all other signature fields unchanged. Old real pressure-coupled runs need
an explicit new protocol (`restart_mode="state_only"`); unrelated mismatches
remain rejected. See [B01 migration](../B01_EXTERNAL_PROVIDER.md#checkpoints-and-public-api-migration).

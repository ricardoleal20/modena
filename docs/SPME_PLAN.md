# SPME implementation plan (dedicated session)

Status: **ready to implement** (2026-09-27). This is the Tier-B centerpiece
for the computational-biology claim (protein + explicit water + PME). The
direct Ewald is measured ~10.8× behind LAMMPS `ewald` at N=216 even after
Ewald-b; the reciprocal FFT path is the only way a 20k-atom explicit-water
benchmark exists at all.

## The acceptance battery (preconditions, from COMMITTEE_ACTA + MATCH_ENGINES)

1. **Energy**: SPME reciprocal == direct-Ewald reciprocal (the current ±k
   path) to a configured accuracy target (default 5e-5 · |E_rec|) on the
   NaCl validation deck, the displaced config, and a lower-κ config.
2. **Forces**: per-atom FD check `F = −∇U` of the SPME Hamiltonian passes
   (the SPME gradient is analytic **of the mesh functional**, so dU/dr = −F
   holds exactly for that Hamiltonian — Essmann-analytic-gradient argument).
3. **NVE**: drift unchanged in class vs the direct deck (~1e-6 at 216 atoms).
4. **Virial**: the reciprocal **strain-form** virial is computed on the grid
   `W_rec = Σ_k Ṽ_k(1 − k²/2κ²)` and closes the documented pressure defect
   (+0.00632 on the NaCl deck, VALIDATION.md §6).
5. **Benchmark**: 216-atom charged deck ≥ 3× faster than Ewald-b's 107 s; the
   FIRST 10k-atom charged run becomes feasible.
6. LAMMPS `pppm` matched-κ/accuracy protocol for the paper table.

## Stage A — energy-only reciprocal (the independent first artifact)

Purpose: prove the FFT route + fix the factor conventions, validated against
the direct ±k path. ~300 lines + tests. **The α-verification rule**: the only
unknown is the grid normalization constant `α` (the FFT/B-spline convention,
physical values k = 2πk̂/L and κ are fixed definitions). Compute SPME energy
with α=1 and fit α to the direct ±k reciprocal energy on the deck; then verify
α is **identical** on a different-κ config (different k-weighting). A stable α
proves the functional form is right; an unstable α means a shape/factor error.

Form (orthorhombic, order-4 B-spline, reduced units):
- Grid per axis N_a = ceil(L_a·k_target/π)... choose κ-mesh: h ≤ π/k_max with
  k_max = deck `kspace_cutoff` (7.0) → h ≤ 0.449 → N ≥ 34/axis (36³ for the
  15³ deck).
- `M4` weights: b_{-1} = (1−u)³/6, b_0 = (3u³−6u²+4)/6, b_1 = (−3u³+3u²+3u+1)/6,
  b_2 = u³/6, node m = (floor(θ)+k mod N), θ = N·f.
- Spread Q[m] += q_i·b(k1)·b(k2)·b(k3); forward FFT (unscaled) → Q̂.
- `b3(k̂)` = Σ_m b_m·e^{−2πi·k̂m/N} per axis (FT of the spline weights).
- E_rec = α · (2π/V) · Σ_{k̂≠0} e^{−k²/4κ²}/k² · |Q̂(k̂)|² / (|b3(k1)b3(k2)b3(k3)|²),
  physical k = 2πk̂·(1/L) per axis; half-space via Hermitian symmetry.
- Implementation: rustfft (6.1+, `process_strided`); 3-pass 1D FFTs (axis-0
  contiguous; axes 1,2 via chunk transposes into contiguous scratch). Real→
  complex not required for energy (full-complex grid ≤ 36³ ≈ 47k → trivial).

## Stage B — forces (the gradient of the SPME functional)

- ik-differentiation along each axis → 3 inverse FFTs (or FFT + conjugate) →
  back-interpolation with the derivative weights b′ → per-atom force.
- Keeps `dU/dr = −F` for the mesh Hamiltonian; FD gate.

## Stage C — virial + production wiring

- Strain-form on the grid (item 4); `[electrostatics] method = "spme"`;
  composite wiring mirrors Ewald; the real-space erfc/self/background code is
  reused verbatim from `ewald.rs`.

## Notes

- rustfft 6.1+ is a pure-Rust, NEON-capable dependency (add to
  `crates/modena-core/Cargo.toml`).
- The direct Ewald remains the reference; SPME is a second `ForceModel`, like
  LAMMPS keeps `ewald` for its `pppm` test suite.
- Post-paper options it unlocks: GPU charge-spread + FFT (Apple's best-shot),
  MSM for extreme core counts.
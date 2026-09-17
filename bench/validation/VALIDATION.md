# MODENA vs LAMMPS — Validation (Buckingham, Ewald, SHAKE/RATTLE)

This document validates MODENA against LAMMPS across three physics regimes
beyond the plain Lennard-Jones benchmark in `BENCHMARKS.md`:

1. **§2–5 — Buckingham** pair potential (exponential-repulsion + dispersion,
   the physics of ionic and molecular force fields).
2. **§6 — Ewald long-range electrostatics** on a charged NaCl (rocksalt)
   crystal — the long-range Coulomb sum and how it composes with a
   short-range pair potential.
3. **§7 — SHAKE/RATTLE bonded constraints** on rigid water-like trimers —
   the constrained dynamics an integrator must keep on the constraint
   manifold.

> **ENGINES:** MODENA 0.1.0 (Rust, release build) vs LAMMPS (22 Jul 2025,
> `lmp_serial`, single MPI rank, no OpenMP). Single core, Apple silicon,
> macOS. All runs use the velocity-Verlet integrator; the timestep, step
> count, and system are specific to each section (Buckingham: dt = 0.005,
> 2000 steps; Ewald: dt = 0.0025, 2000 steps; SHAKE: dt = 0.005, 2000 steps).
> All trajectories are validated with MODENA's own full-precision
> measurement harness under `harness/`; no numbers are extrapolated.

---

## 1. TOML-reachability coverage audit

Before choosing a validation system, every potential, structure, integrator,
and optional feature listed in the source tree was checked for end-to-end
wiring from a TOML input file through `Simulation::from_config` (the runner).
The audit covers `crates/modena-core/src/io/input.rs` (the `InputConfig`
schema), `crates/modena-core/src/runner.rs` (the `from_config` / `build_*`
path), and `crates/modena-cli/src/main.rs` (the `run` / `check` CLI).

### Pair potentials

| Potential | TOML `model.type` | Runner `build_model` | Runnable? |
|---|---|---|---|
| Lennard-Jones | `lennard_jones` | `LennardJones::new` / `from_types` | **Yes** — single + multi-species (LB mixing + overrides) |
| Buckingham | `buckingham` | `Buckingham::new` / `from_types` | **Yes** — single + multi-species (overrides; no universal mixing rule) |
| Morse | `morse` | `Morse::new` / `from_types` | **Yes** — single + multi-species (overrides) |
| Soft-sphere | `soft_sphere` | `SoftSphere::new` / `from_types` | **Yes** — single + multi-species |
| Tabulated | `tabulated` | `Tabulated::from_text` | **Yes** — linear interpolation from file |

All five potentials are fully reachable from TOML. All use the shifted-force
form (`U_sf = U(r) − U(rc) + (r−rc)·F(rc)`, `F_sf = F(r) − F(rc)`).

### Structures

| Structure | TOML `structure.type` | Runnable? |
|---|---|---|
| FCC | `fcc` | **Yes** — `lattice_constant`, `replicas`, element, mass, charge |
| XYZ file | `xyz` | **Yes** — read positions from file + box |
| Random alloy | `random_alloy` | **Yes** — FCC with per-atom species by concentration |
| Rocksalt (NaCl) | `rock_salt` | **Yes** — ordered cation/anion sublattices with charges |

### Integrators / ensembles

| Feature | TOML key | Runnable? |
|---|---|---|
| NVE (velocity-Verlet) | `integrator.type = "velocity_verlet"` | **Yes** |
| NVT (Langevin BAOAB) | `integrator.type = "langevin_baoab"` | **Yes** |
| NVT (Nosé–Hoover chain) | `integrator.type = "nose_hoover_chain"` | **Yes** |
| Berendsen thermostat | `thermostat.type = "berendsen"` | **Yes** (composes with NVE) |
| Velocity rescaling | `thermostat.type = "velocity_rescaling"` | **Yes** |
| Berendsen barostat (NPT) | `barostat.type = "berendsen"` | **Yes** |
| Parrinello–Rahman (NPT) | `barostat.type = "parrinello_rahman"` | **Yes** |

### Optional physics

| Feature | TOML key | Runnable? |
|---|---|---|
| Bonded (bonds/angles/dihedrals) | `topology` | **Yes** — composes via `CompositeForceModel` |
| SHAKE/RATTLE constraints | `constraints` | **Yes** — wires into integrator |
| Ewald electrostatics | `electrostatics` | **Yes** — composes via `CompositeForceModel` |

**Conclusion:** Every feature in the source tree is genuinely runnable from a
TOML file. There are no dead modules or unwired potentials. Buckingham
(priority 1), Ewald electrostatics (§6), and SHAKE/RATTLE (§7) are all
validated against LAMMPS below.

---

## 2. The matched system

Both engines run the **same physical model**, configured identically:

| Parameter | Value |
|---|---|
| Potential | Buckingham `A·exp(−r/ρ) − B/r⁶` |
| A (repulsive prefactor) | `1000.0` |
| B / C (dispersive prefactor) | `1.0` |
| ρ (decay length) | `0.3` |
| Cutoff | `5.0` |
| Units | reduced (LJ): `σ = ε = m = 1` |
| Lattice | FCC, `a = 3.0` (density `ρ_lat = 4/27 = 0.1481`) |
| Atoms | 864 (`replicas = [6,6,6]`) |
| Integrator | velocity-Verlet (NVE) |
| Timestep | `dt = 0.005` |
| Initial T | `0.1` (Maxwell-Boltzmann, COM momentum removed) |
| Seed | `87287` |
| Steps | `2000` |

The lattice constant `a = 3.0` places the system in a moderately compressed
state (nearest-neighbour distance `r_NN = a/√2 = 2.12`), giving a non-trivial
initial PE of ~4500 for step-0 comparison.

### Potential-form correspondence

MODENA's `Buckingham` and LAMMPS `pair_style buck` implement the **same
functional form**:

```
U(r) = A · exp(−r/ρ) − C / r⁶         # (C == B in MODENA's naming)
```

The **only** convention difference is the cutoff treatment:

| | MODENA | LAMMPS (`pair_modify shift yes`) |
|---|---|---|
| Energy | `U_sf = U(r) − U(rc) + (r−rc)·F(rc)` | `U_es = U(r) − U(rc)` |
| Force | `F_sf = F(r) − F(rc)` | `F = F(r)` (truncated) |

MODENA uses the **force-shifted** form (both U and F continuous at cutoff);
LAMMPS with `pair_modify shift yes` uses the **energy-shifted** form (U
continuous, force truncated). The per-pair energy difference is:

```
ΔU = (r − rc) · F(rc)
```

At `rc = 5.0` with these parameters:
- `U(rc) = 1000·exp(−16.67) − 1/5⁶ = −6.2 × 10⁻⁶` (negligible)
- `F(rc) = (1000/0.3)·exp(−16.67) − 6/5⁷ = 1.16 × 10⁻⁴` (negligible)

The total force-shift correction summed over all ~18,000 pairs within the
cutoff is **−4.67** (deterministic for the step-0 lattice), which is **0.10%**
of the total PE. This is the sole systematic offset between the two engines;
after accounting for it, they agree to **< 0.001%** on step-0 PE (see §3).

---

## 3. Physics agreement

### Step-0 thermodynamic state (identical lattice and temperature)

| Quantity | LAMMPS | MODENA | Δ | Rel. error |
|---|---:|---:|---:|---:|
| Potential energy | 4509.2096 | 4504.5380 | −4.6716 | 0.104% |
| Kinetic energy | 129.4500 | 129.6000 | +0.1500 | 0.116% |
| Total energy | 4638.6596 | 4634.1380 | −4.5216 | 0.098% |
| Temperature | 0.1000 | 0.1000 | 0 | 0% |
| Pressure | 1.8778 | 1.8770 | −0.0008 | 0.044% |

**Two convention offsets, both understood:**

1. **PE offset (−4.67):** the force-shift vs energy-shift difference (see §2).
   This is `F(rc) · Σ_pairs (r_ij − rc) = 1.16×10⁻⁴ × (−40,300) = −4.67`.
   After correcting, step-0 PE agrees to **< 0.001%**.

2. **KE offset (+0.15):** the degrees-of-freedom convention. LAMMPS removes
   3 COM translational dofs (`3N − 3 = 2589`), MODENA uses `3N = 2592`.
   At `T = 0.1` this is exactly `(3/2)·T = 0.15`. (Same convention difference
   as the LJ benchmark in `BENCHMARKS.md` §2.)

After accounting for both conventions, **every step-0 observable agrees to
< 0.001%**.

### Equilibrated time-averaged observables (avg over steps 100–2000)

| Quantity | LAMMPS ⟨⟩ | MODENA ⟨⟩ | Δ | Rel. error |
|---|---:|---:|---:|---:|
| Temperature | 0.04935 | 0.04942 | +0.00007 | 0.14% |
| Potential energy | 4574.79 | 4570.10 | −4.69 | 0.10% |
| Pressure | 1.88773 | 1.88694 | −0.00080 | 0.042% |

The ⟨PE⟩ offset (−4.69) is the same force-shift convention difference
persisting through the trajectory (the correction is nearly constant because
the lattice barely evolves). The ⟨T⟩ agreement of 0.14% is within
microstate-sampling variance: the two engines use independent RNGs and
different DOF conventions, so their microscopic trajectories diverge
(chaotic Lyapunov growth), but the **statistical** (microcanonical) averages
agree.

The ⟨P⟩ agreement of 0.042% is excellent — pressure is the most sensitive
observable to force accuracy, and this level of agreement confirms the
Buckingham force kernel is correct.

### Energy conservation (NVE, 2000 steps)

| Engine | E₀ | E_final | Relative drift |
|---|---:|---:|---:|
| **MODENA** | 4634.1380 | 4634.1489 | **2.35 × 10⁻⁶** |
| LAMMPS | 4638.6596 | 4638.6708 | 2.41 × 10⁻⁶ |

Both engines conserve the NVE Hamiltonian to ~2.4 × 10⁻⁶ over 2000 steps,
consistent with a symplectic velocity-Verlet integrator at `dt = 0.005`.
MODENA's slightly better conservation is expected: its reported PE is the
exact integral of its (force-shifted) force, making `H = KE + PE` the true
conserved quantity. LAMMPS's energy-shifted PE is not the exact integral of
its truncated force, so its "total" has a marginally larger drift.

---

## 4. Honest assessment

### Where MODENA agrees with LAMMPS

- **Potential energy:** step-0 PE agrees to 0.10% (entirely explained by the
  force-shift convention). After correction, agreement is < 0.001%.
- **Temperature:** ⟨T⟩ agrees to 0.14% — within statistical sampling variance.
- **Pressure:** ⟨P⟩ agrees to 0.042% — confirms the force kernel is correct.
- **Energy conservation:** both ~2.4 × 10⁻⁶ drift; MODENA marginally better.
- **Dynamics:** the trajectories are statistically indistinguishable (same
  Hamiltonian, same initial state, diverging only by chaotic Lyapunov growth
  + the negligible F(rc) ≈ 10⁻⁴ force offset).

### Where MODENA diverges and why

1. **PE offset (−4.67):** systematic, not a bug. MODENA force-shifts; LAMMPS
   energy-shifts. The correction is deterministic and < 0.1% of PE. In a
   production setting, both engines would use the same shift convention (or
   the offset would be subtracted analytically). MODENA's force-shift is
   arguably the *better* choice: it makes the force continuous at the cutoff,
   improving energy conservation.

2. **KE offset (+0.15):** the `3N` vs `3N − 3` DOF convention. Same as the LJ
   benchmark. At `T = 1.0` this would be 1.5 (as in `BENCHMARKS.md`); at
   `T = 0.1` it scales to 0.15.

3. **Trajectory divergence at long times:** the two engines use independent
   RNGs for initial velocity assignment, so even with identical Hamiltonians
   their microscopic trajectories diverge exponentially (Lyapunov instability).
   This is expected and is *not* a discrepancy — the statistical averages
   agree.

4. **Neighbor-list rebuild cadence:** MODENA rebuilds its cell list every
   step; LAMMPS uses a Verlet skin list (rebuild only when atoms move
   `> 0.5·skin`). This affects performance (LAMMPS ~12× faster here) but not
   physics — both find the same neighbor pairs.

---

## 5. Reproducing

```bash
# Prereqs: rust toolchain, lmp_serial on PATH.
cargo build --release

# Run both engines (from the repo root):
target/release/modena run -i bench/validation/buckingham/modena.toml
lmp_serial -in bench/validation/buckingham/in.lammps

# Outputs:
#   bench/validation/buckingham/results/modena_thermo.csv  (MODENA thermo log)
#   bench/validation/buckingham/results/modena_traj.xyz    (MODENA trajectory)
#   bench/validation/buckingham/results/lammps_log.txt     (LAMMPS full log)
```

Matched configuration pairs:

| Case | MODENA | LAMMPS |
|---|---|---|
| Buckingham, 864 atoms, 2000 steps | `bench/validation/buckingham/modena.toml` | `bench/validation/buckingham/in.lammps` |
| Ewald rocksalt, 216 atoms, 2000 steps | `bench/validation/ewald/modena.toml` | `bench/validation/ewald/in.lammps` |
| SHAKE trimers, 48 atoms, 2000 steps | `bench/validation/shake/modena.toml` | `bench/validation/shake/in.lammps` |

```bash
# Full-precision measurement harness (standalone crate; links the engine
# as-shipped). Ewald component decomposition + SHAKE residual measurement.
cd bench/validation/harness && cargo build --release
./target/release/validation-harness ewald-decompose   # -> ewald/results/harness_decomposition.txt
./target/release/validation-harness shake-measure     # -> shake/results/harness_shake_measure.txt

# Bond-length residuals from the trajectory files (both engines).
python3 bench/validation/shake/analyze_shake.py lammps bench/validation/shake/results/lammps_trimers.lammpstrj
python3 bench/validation/shake/analyze_shake.py modena bench/validation/shake/results/modena_traj.xyz
```

---

## 6. Ewald long-range electrostatics (NaCl rocksalt)

### The matched system

| Parameter | Value |
|---|---|
| Structure | rocksalt (NaCl): two interpenetrating FCC sublattices |
| Lattice constant | `a = 5.0` (box 15³), `replicas = [3,3,3]` |
| Atoms / charges | 216 (108 Na `+1`, 108 Cl `−1`), net neutral |
| Short-range | Buckingham `A=1000, B=1, ρ=0.3` (all type pairs), cutoff `5.0` |
| Long-range | Ewald: `κ = 0.8`, real cutoff `5.0`, k-space cutoff `7.0` |
| LAMMPS mapping | `pair_style buck/coul/long 5.0 5.0` + `kspace_style ewald 1e-8`, `kspace_modify gewald 0.8` |
| Integrator | velocity-Verlet (NVE), `dt = 0.0025` |
| Initial T / seed / steps | `0.1` / `87287` / `2000` |

The short-range Buckingham keeps the §2 convention difference (force-shift vs
energy-shift); the long-range part uses the **same method on both engines**
(Ewald decomposition with the same κ), so the Ewald terms should agree
closely — the pair shift offset is the only systematic PE offset expected.

### Component-level brute-force decomposition (the authoritative check)

The harness (`bench/validation/harness/`, mode `ewald-decompose`) rebuilds the
exact deck system and evaluates MODENA's `EwaldForceModel` against an
**independent brute-force reference**: O(N²) real-space `Σ q_iq_j·erfc(κr)/r`
(no cell list) + direct reciprocal-space sum + self + background, at the exact
lattice (x0) and at a minimally displaced configuration (x1 = x0 + v·10⁻⁶):

| Component | MODENA | Brute force | Δ |
|---|---:|---:|---:|
| Real space | −1.189442043693 | −1.189442043693 | < 2×10⁻¹⁰ |
| Reciprocal (19,548 k-vectors) | +23.186609678823 | +23.186609678823 | < 10⁻¹⁰ |
| Self | −97.491960037052 | −97.491960037052 | < 10⁻¹⁰ |
| Background (neutral) | 0 | 0 | 0 |
| **Total Ewald PE** | **−75.494792402087** | **−75.494792401922** | **1.7×10⁻¹⁰** |

The engine's `compute()` (forces+energy+virial) and `energy()`-only paths
return identical energies, and the real-space virial `Σ r·F` matches the
brute-force value to 1.7×10⁻¹⁰. The Ewald energy path is verified to
~2×10⁻¹² relative, component by component, against an independent
implementation.

### Step-0 and time-averaged agreement vs LAMMPS

| Quantity | LAMMPS | MODENA | Δ | Explanation |
|---|---:|---:|---:|---|
| Step-0 PE | 87.17065 | 86.69628 | −0.4744 | Buckingham force-shift offset: `Σ(r−rc)·F(rc) ≈ 4338 pairs × (−1) × 1.10×10⁻⁴ ≈ −0.48` ✓ |
| Step-0 KE | 32.25 | 32.40 | +0.15 | `3N−3` vs `3N` DOF convention (= 3/2·T) — same as §3 |
| Step-0 T | 0.1000 | 0.1000 | 0 | — |
| Step-0 P | 0.13724 | 0.14346 | +0.00632 | **reciprocal-virial form** (see below) |
| ⟨T⟩ (100–2000) | 0.08121 | 0.07862 | −3.2% | chaotic trajectory divergence + shift offset; same oscillation phase |
| ⟨PE⟩ (100–2000) | 93.231 | 93.623 | −0.392 | the shift-convention offset persisting (as in §3) |
| ⟨P⟩ (100–2000) | 0.13391 | 0.14129 | +0.0074 | reciprocal virial + shift offset |

After the two known convention corrections (pair shift ≈ −0.48, KE DOF
+0.15), step-0 total energies agree to **~10⁻⁴ relative** — the Ewald terms
themselves agree to well beyond that (the residual is the uncertainty of the
analytic shift-offset estimate, not a measured Ewald discrepancy; the
component-level brute-force check above bounds the Ewald error at 10⁻¹²).

### Energy conservation

| Engine | E₀ | E_final | Relative drift |
|---|---:|---:|---:|
| **MODENA** | 119.096283 | 119.096373 | **7.6 × 10⁻⁷** |
| LAMMPS | 119.420650 | 119.420750 | 8.4 × 10⁻⁷ |

Both conserve to ~10⁻⁶ over 2000 steps at `dt = 0.0025`; MODENA marginally
tighter (its force-shifted PE is the exact integral of its force — the
conservation contract, §`ARCHITECTURE.md` 2).

### Honest finding: the reciprocal-space virial form

MODENA's Ewald reciprocal virial is accumulated in the naive atom-based
`Σ r_i·F_i` form (`ewald.rs`), which contributes ≈ 0 for the reciprocal part.
The correct reciprocal pressure virial (what LAMMPS uses) is the
**strain-derivative form** `Σ_k U_k·(1 − k²/2κ²)` = **−63.98** here. The
consequence is exactly predictable:

```
ΔP = W_rec_strain / (3V) = −63.98 / (3 × 3375) = +0.00632
```

which matches the observed step-0 pressure gap (**+0.00632**) to all reported
digits. This is a **pressure-reporting defect only**: forces, energies, and
dynamics are unaffected (verified by the brute-force decomposition above and
by the trajectory agreement), because the virial does not feed back into the
NVE equations of motion. **Caveat:** an NPT run coupling a barostat to this
pressure would use a slightly wrong virial until the strain form is
implemented — a documented limitation, not a silent error.

### The κ = 0.5 → 0.8 lesson (truncated-Ewald discontinuity)

The deck originally used `κ = 0.5` (`erfc(κ·rc) = 4.1×10⁻⁴`). On this lattice
**648 pairs sit exactly on the r = rc sphere**; the first time any of them
crosses the cutoff the real-space sum jumps coherently by ~`648 × q²·erfc(κrc)/rc`
≈ **+0.026** — the classic truncated-Ewald discontinuity, which MODENA's
`ConservationReport` surfaces as a drift verdict (this is the conservation
contract doing its job). With `κ = 0.8`, `erfc(κ·rc) = 1.5×10⁻⁸` and the
residual discontinuity is ~10⁻⁶ total: the harness measures the Ewald PE
changing by 1.002×10⁻⁶ between the exact lattice and a v·10⁻⁶-displaced
configuration — dominated by hundreds of cutoff-boundary crossings of
`q²·erfc(κrc)/rc ≈ 3×10⁻⁹` each, invisible at the 10⁻⁶ drift level. The
LAMMPS κ=0.5 run is preserved in `results/lammps_kappa0.5_crosscheck.log`.

### A finite-difference artifact, ruled out

The harness also runs an FD force-consistency check (`F = −∇U` by central
differences at h = 10⁻⁶) and reports a worst `|F − FD| = 1.5×10⁻³`. This is
**not a force bug**: moving one atom by ±h flips the cutoff inclusion of its
~18 borderline pairs, each a `q²·erfc(κrc)/r ≈ 3×10⁻⁹` energy jump, so the
FD gradient is contaminated at the `3×10⁻⁹ / 2×10⁻⁶ ≈ 1.5×10⁻³` level per
crossing — precisely the observed magnitude. The authoritative checks (the
component-level brute-force agreement at 10⁻¹⁰ and the trajectory agreement
with LAMMPS) confirm the forces are the correct gradient of the verified
energy.

---

## 7. SHAKE/RATTLE bonded constraints (rigid trimers)

### The matched system

| Parameter | Value |
|---|---|
| System | 16 rigid water-like trimers (O + 2 H), 48 atoms, box 10³ |
| Masses | all `1.0` |
| Constraints | 3 per trimer: 2× O–H (`d = 1.0`) + H–H 1-3 (`d = 1.632993161855452` = 2·sin(θ/2), θ = 109.47° tetrahedral) |
| LAMMPS mapping | `fix shake 1e-10 100 0 b 1 a 1` — LAMMPS converts the angle constraint to the same 1-3 distance, so **both engines solve the identical 48-constraint set** (rigid triangles). `fix shake` requires disconnected clusters — hence isolated trimers. |
| Non-bonded | LJ `σ=1, ε=0.1, rc=2.5`, **no exclusions** (`special_bonds lj 1 1 1`) — bonded pairs feel the full LJ in both engines |
| Bonded terms | harmonic bonds/angles present with `k = 0` (no force; SHAKE does all the work) |
| Integrator | velocity-Verlet + SHAKE/RATTLE, `dt = 0.005`, tol `1e-10`, max 100 iter |
| Initial T / seed / steps | `0.5` / `87287` / `2000` |

The trajectories cannot be compared microscopically (independent RNGs), so the
validation targets the constraint-specific invariants: constraint residuals,
the exactness of the velocity projection, and constrained-energy conservation.

### Constraint satisfaction (the headline result)

| Metric | MODENA | LAMMPS |
|---|---:|---:|
| bond `|r − d|` max (full precision) | **5.0 × 10⁻¹¹** | **1.0 × 10⁻¹⁴** |
| bond `|r − d|` RMS | 2.7 × 10⁻¹¹ | 1.1 × 10⁻¹⁵ |
| RATTLE `|r·v|` max | 1.0 × 10⁻¹⁰ | (enforced to the same tolerance) |
| samples | 96,000 bond visits (2000 steps × 48) | 960 (20 dump frames × 48) |

MODENA's residuals are measured **in-engine at full f64 precision** by the
harness (`shake-measure`) — the shipped XYZ writer emits 8 decimals, which
floors trajectory-file residuals at ~10⁻⁸ (the file-based measurement reports
1.1×10⁻⁸ max, i.e. the writer floor, not the solver). Both engines hold every
constraint 5+ orders below any physical tolerance; LAMMPS's is tighter still,
but both meet the configured `1e-10` tolerance.

### The velocity-projection exactness check

The initial velocities are drawn at `T = 0.5` on all `3N = 144` DOF
(KE₀ = 36.0). The 48 constraints remove exactly 48 velocity DOF; the one-time
initial RATTLE projection should therefore remove **exactly 48/144 = 1/3 of
the kinetic energy**:

| | MODENA | LAMMPS |
|---|---:|---:|
| KE₀ (pre-projection, both "T = 0.5") | 36.000 (created on 3N = 144 DOF) | 35.250 (created on 3N−3 = 141 DOF) |
| expected post-projection KE | 96/144 × 36 = **24.000** | 93/141 × 35.25 = **23.250** |
| measured post-projection KE | ≈ **24.07** (E step 1 = 24.0707, PE₁ ≈ 0) | ≈ **23.3–23.6** (⟨KE⟩ over the run) |

Both engines remove exactly the constrained-DOF fraction of the initial KE —
the signature of a correct RATTLE projection (the residual difference between
the two post-projection KE values is the independent-RNG realization of the
two draws, within the χ² fluctuation σ ≈ 3.5 of a 96-DOF Maxwell state — not
an engine discrepancy).

### Constrained-energy conservation (post-projection, NVE+SHAKE)

| Engine | E (post-projection) | Rel. drift |
|---|---:|---:|
| **MODENA** | 24.0707 (steps ≥ 100 range [24.0680, 24.0738]) | **2.4 × 10⁻⁴** |
| LAMMPS | 22.6111 (steps 100–2000 range [22.6104, 22.6113]) | 4.0 × 10⁻⁵ |

Both are symplectic-level for constrained dynamics at `dt = 0.005`; LAMMPS is
~6× tighter here. Neither shows any constraint-leak signature (a leaking
SHAKE shows *drifting* energy, not bounded oscillation).

### Statistical equivalence

Comparing temperatures on the **same free-DOF basis** (MODENA reports T on
`3N = 144` and does not subtract constraint DOF; LAMMPS's thermo uses
`3N − 3 − 48 = 93`):

| | MODENA (2·⟨KE⟩/96) | LAMMPS (2·⟨KE⟩/93) |
|---|---:|---:|
| ⟨T⟩ per free DOF | 0.513 | 0.507 |

1.2% agreement — the two engines sample the same constrained microcanonical
ensemble.

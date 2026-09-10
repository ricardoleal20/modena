# MODENA — Novelty & Features Proposal for Journal Publication

**Status:** Proposal for review (no implementation yet). SOTA assessment is from
the author's knowledge of the field (web search was unavailable during drafting);
items marked **[VERIFY]** should be confirmed against live literature before a
paper is written.

**Date:** 2026-09-10

---

## 0. The question

MODENA today is a correct, well-tested, LAMMPS-class MD core in Rust: shifted-force
LJ, cell-list neighbors, velocity-Verlet + BAOAB-Langevin, Berendsen/Bussi
thermostats, NVE/NVT, rayon parallelism, and a `ForceModel` trait seam for ML.
Benchmarked vs LAMMPS it is ~7× slower (no Verlet skin yet, roadmap M7) but uses
~3× less RAM and conserves NVE energy ~6× better (drift 4.7e-6 vs 2.7e-5).

The question this document answers: **what genuinely novel, research-worthy
features distinguish MODENA from LAMMPS / GROMACS / HOOMD-blue / OpenMM and
justify a journal paper?**

The honest starting point: *"rewritten in Rust"* is an engineering contribution,
not a scientific one, and will not clear the novelty bar of a research journal on
its own. The novelty must be in *what the engine can do that others cannot*, not
the language it is written in.

---

## 1. SOTA landscape (what already exists)

| Capability | SOTA tools / papers | Saturation |
|---|---|---|
| **Differentiable MD** | **JAX-MD** (Schoenholz & Cubuk, 2021) — end-to-end differentiable, hardware-accelerated; **TorchMD / TorchMD-Net** | **Saturated.** JAX-MD owns this space. Rust has no reverse-mode autodiff comparable to JAX, so a Rust differentiable MD is both crowded and technically weak. **Not a viable novelty vector.** |
| **Equivariant ML force fields** | **NequIP** (Batzner 2022), **MACE** (Batatia 2022), **Allegro** (Musaelian 2023, *Nature Comm.*), **MACE-MP-0** foundation models | **Crowded at the method level.** The science is done. A Rust-native ML potential would be an *engineering* contribution, not a *scientific* one. |
| **Active-learning MD** | **FLARE** (Vandermause 2022) — on-the-fly GP force fields with uncertainty-driven training; **Active Learning for MLPs** surveys | **Moderately crowded.** The methodology exists; a safe, fast Rust core with a `ForceModel` hot-swap seam is a real but incremental contribution. |
| **Agent / LLM-driven science** | **Coscientist** (Boiko 2023, *Nature*) — GPT-4 autonomous synthesis; **ChemCrow** (Bran 2024, *Nature Mach. Intell.*) — LLM chemistry agent; **ChemAgent**, **ProtAgent** | **Emerging & under-served for *simulation*.** These agents drive *experiments/chemistry*, not MD engines. No MD engine ships a native programmatic/agent control interface — agents shell out to LAMMPS/GROMACS as black boxes. **This is the open gap.** |
| **GPU MD** | **HOOMD-blue**, **OpenMM**, LAMMPS GPU | Saturated. MODENA should add GPU support eventually but it is not a novelty. |
| **Browser/WASM MD** | **NGL**, **Mol\*** (viewers, not engines); scattered WebMD demos | **Sparse.** Full production-grade MD running in-browser via WASM is rare; low-memory laptop-scale MD is a defensible *engineering* angle, not primary science. |
| **Energy conservation / shifted-force** | **Toxvaerd & Dyre** (2011) — shifted-force LJ, energy conservation analysis | Known methodology. "Better drift" is a *consequence* of a consistent Hamiltonian, not a new idea — but *formalizing* it as a construction principle + property-based testing is a legitimate methodology contribution. |
| **Rust in HPC/science** | Growing interest, few *published* Rust MD frameworks | Novelty as an artifact, not as science. |

### Key conclusion

The single most defensible, timely, under-served novelty vector is **agent-native
MD**: an MD engine that exposes a native, structured control interface so that
LLM agents can autonomously set up, parameterize, run, and analyze simulations —
not by shelling out to a black-box binary, but by calling typed simulation
primitives. This sits squarely in the 2024–2026 "AI-driven science" wave
(Coscientist, ChemCrow) but is *unoccupied* on the simulation-engine side.

---

## 2. Candidate novelty features (ranked)

### 🥇 Feature A — Agent-native simulation interface (MCP tool layer)
**The primary novelty.** Expose MODENA's simulation primitives as a structured
tool interface (Model Context Protocol, or an equivalent typed RPC) so an LLM
agent can:
- construct a system from a natural-language / structured spec ("FCC argon,
  864 atoms, ρ=0.8442"),
- run short NVE/NVT probes and read back thermodynamics,
- launch parameter sweeps (temperature ladder, density scan) as parallel streams,
- detect anomalies (phase transition, energy drift, melting) and adapt the sweep,
- checkpoint, branch, and compare simulations.

**Why novel:** No existing MD engine ships a native agent-control surface.
Coscientist/ChemCrow drive *chemistry experiments*; for MD they would have to
generate LAMMPS input scripts blindly. MODENA turns the engine itself into an
agent-callable tool with typed, validated primitives.

**Why MODENA-specific:** MODENA's low RAM (~3× less than LAMMPS) and fast startup
make *many parallel agent-driven simulation streams* economically feasible on a
single workstation — exactly the regime agent sweeps need. LAMMPS' 12 MB
baseline × 50 concurrent streams is wasteful; MODENA's 3–6 MB makes it tractable.

**Paper thesis enabled:** *"Agent-native molecular dynamics: LLM-orchestrated
simulation discovery over a safe, lightweight engine."*

**Effort:** Medium. Needs an MCP/RPC server crate (`modena-mcp`) wrapping the
existing `Simulation` API + a schema for simulation primitives. No core physics
change. The `runner.rs` API is already close to what an agent would call.

**Risk:** Novelty hinges on a compelling *demonstration* (an agent autonomously
discovering, e.g., a phase transition or an optimal annealing schedule). The
engineering is straightforward; the paper lives or dies on the demo.

---

### 🥈 Feature B — Consistent-Hamiltonian energy conservation, formalized + property-tested
Turn MODENA's existing 6×-better-drift result into a *methodology contribution*:
- Make **"energy is the exact integral of force"** a *construction principle*
  enforced by the `PairPotential` trait (the trait returns `(U, f_over_r)` from a
  single analytic kernel — already true; document and test the invariant).
- Add **property-based tests** (`proptest` / `quickcheck`) that verify, for random
  configurations: (1) finite-difference force ≈ analytic force, (2) finite-difference
  energy gradient ≈ analytic force, (3) symplectic one-step energy drift bounded
  by O(dt⁴), (4) no NaN/Inf under random inputs.
- Publish the LAMMPS-vs-MODENA energy-drift benchmark as a reproducible figure.

**Why novel:** Energy-conservation *results* are common; making *constructive
force/energy consistency + property-based verification* the advertised guarantee
of an engine is not. It reframes a known property (shifted-force → better drift)
as an engineered, tested contract.

**Why MODENA-specific:** Builds directly on the existing shifted-force kernel and
the 6× drift result. Pure upside, low risk.

**Paper thesis enabled:** Supports Feature A as a "trust" pillar — agents can rely
on the engine's conservation guarantee when interpreting drift as a signal.

**Effort:** Low–Medium. Mostly tests + a trait doc/invariant + benchmark figure.

---

### 🥉 Feature C — Active-learning ML force fields on the `ForceModel` seam
Implement the `ForceModel` trait with an ML backend (start with a lightweight
equivariant or Behler-Parrinello-style network, or bind to an existing model via
ONNX/Candle) plus an **active-learning loop**: run with the cheap classical
potential, query an uncertainty estimator, and retrain the ML surrogate on-the-fly
for regions of configuration space that matter.

**Why (moderately) novel:** The methodology (FLARE) exists; the contribution is a
*safe, single-binary, fast-startup Rust core* that hot-swaps classical ↔ ML
potentials behind one trait — attractive for the agent-driven sweeps of Feature A
(where you want a cheap exploratory potential and an expensive accurate one).

**Effort:** High. ML backend + uncertainty + training loop is a large body of work.
**Recommendation:** implement only the `ForceModel`-trait ML stub + a classical/ML
hot-swap demo for v1; defer real active-learning to a follow-up. Do not let this
block the paper.

---

### Feature D — WASM / browser-deployable MD (enabler, not headline)
Compile `modena-core` to `wasm32-unknown-unknown` and run interactive MD in the
browser, streaming trajectory to a Mol*/NGL viewer. Supports Feature A's "agent
shows the user a live simulation" story.

**Effort:** Medium. Rust→WASM is straightforward; the work is the viewer streaming
+ UI. **Secondary.**

---

### Feature E — Verlet skin neighbor lists (roadmap M7, not novelty)
Closes the 7× perf gap vs LAMMPS. **Engineering, not novelty** — mention as
ongoing work, do not sell as a contribution.

---

## 3. Recommended paper thesis

**Title (working):** *"MODENA: An Agent-Native, Energy-Conserving Molecular
Dynamics Framework for Autonomous Simulation Discovery"*

**One-paragraph pitch:**

> We present MODENA, a molecular dynamics engine whose distinguishing design goal
> is *agent nativeness*: the simulation itself is exposed as a set of typed,
> validated primitives callable by large-language-model agents, enabling
> autonomous setup, execution, and adaptive analysis of simulation campaigns
> without blind shell-out to a legacy binary. Two properties make this tractable
> and trustworthy: (1) a low memory footprint (~3× below LAMMPS) that makes many
> concurrent agent-driven streams economical on a single workstation, and (2) a
> *constructively energy-conserving* force/energy kernel — the reported potential
> energy is the exact integral of the reported force — giving ~6× lower NVE drift
> than LAMMPS and a drift signal an agent can interpret as a correctness check.
> We validate physics agreement with LAMMPS to <1e-4 and demonstrate an LLM agent
> autonomously discovering [a phase transition / an optimal annealing schedule]
> over MODENA.

**Novelty claims (in priority order):**
1. First MD engine with a native agent/MCP control interface — agents call typed
   simulation primitives, not scripts. **[Primary]**
2. Constructive force/energy consistency as an engine contract, verified by
   property-based testing; ~6× lower NVE drift than LAMMPS. **[Methodology]**
3. Low-RAM core enabling dense parallel agent-driven simulation streams. **[Enabler]**
4. (Optional) Classical↔ML potential hot-swap behind a single trait for cheap
   exploration + accurate confirmation. **[Stretch]**

**Honest non-claims:** we do *not* claim novelty in differentiable MD (JAX-MD),
equivariant ML potentials (NequIP/MACE/Allegro), or GPU MD (HOOMD/OpenMM).

---

## 4. Suggested build order (for after approval)

1. **Feature B** (low risk, low effort) — property tests + conservation contract +
   benchmark figure. De-risks the "trust" pillar immediately.
2. **Feature A** (primary novelty) — `modena-mcp` crate wrapping `Simulation`;
   define the primitive schema; build one compelling agent demo. This is the paper.
3. **Feature D** (secondary) — WASM build + live viewer, in service of the demo.
4. **Feature C** (stretch) — `ForceModel` ML stub + hot-swap demo only; defer real
   active learning.
5. **Feature E** (M7) — Verlet skin, shipped alongside as a performance note.

---

## 5. Items to verify with live literature before writing [VERIFY]

- Confirm no MD engine has shipped a native agent/MCP/tool-calling interface
  since training cutoff (search: "MCP molecular dynamics", "LLM agent MD engine",
  "agent-driven simulation 2025 2026").
- Confirm Coscientist/ChemCrow do not already orchestrate MD engines via a typed
  interface (they may have evolved).
- Confirm the FLARE / active-learning MD landscape hasn't produced a Rust-native
  competitor.
- Confirm the exact JAX-MD / TorchMD positioning to frame the non-claim correctly.
- Check whether a "consistent Hamiltonian / property-tested symplectic MD" angle
  has been published as an engine contribution (likely not, but verify).
- Verify MODENA's name isn't already taken by a published MD code.

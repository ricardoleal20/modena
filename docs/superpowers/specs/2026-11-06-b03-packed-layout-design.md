# B03 packed pair layout experiment

The user authorized all remaining backlog cases. B03 compares internal SoA,
AoSoA and spatially clustered views with B02's interleaved pair path. Public
System storage, logical atom identity, topology, restart and serial reference
remain the reference. Experimental code initially lives only in modena-bench.

The first deliverable is a repeatable full-step experiment supporting an
accept/reject decision, not an asserted speedup. Compare baseline PairForceModel
with SoA and width-four AoSoA views. Add a spatial permutation and width-four
cluster mask view derived from the exact half neighbor list. Sorted lanes carry
logical IDs; padding never creates an interaction. Cluster masks must reproduce
the original scaled/excluded logical pair set exactly, once per half pair.
Packing reads current coordinates on every force evaluation. Cluster list
rebuild follows neighbor rebuild epochs; masks cannot become stale at fixed N.

Conservative first domain: orthorhombic periodic single-species shifted-force
LJ, sigma=epsilon=mass=1, cutoff2.5, skin0.3. Unsupported cells/potentials must
reject explicitly in experiment; generic/scaled/excluded logical-pair fixtures
still exercise mapping and masks. Do not change production System or infer
universal SIMD acceleration from a layout. Scalar lane evaluation is allowed
but must be labeled; no unsafe intrinsics are necessary.

Numerical gates fixed before running: energy/virial absolute1e-8 plus relative
1e-12, each force component absolute1e-10 plus relative1e-12 against independent
scalar half-pair calculation. Identity/permutation and mask sets exact. Check
cutoff boundaries, periodic crossings, empty/padded clusters, coordinate and box
changes, reorder, resize, exclusions/scales and additive preexisting forces.
Timestep convergence must retain B06's existing observations; changed parallel
sum order need not produce bit-identical trajectories, but the production
one-thread/deterministic path must remain unchanged.

Compare full velocity-Verlet steps, including position packing, force scatter,
neighbor rebuild and spatial sorting/mask construction where needed. Separate
setup, no-rebuild and forced-rebuild timings so hidden conversion costs cannot
be excluded. Frozen release binaries with source/compiler/input/output hashes,
alternating AB/BA five paired repetitions plus separate warmups, N864/2048/5000/
10000, 1/4/8 threads and actual default count. Record wall/CPU/RSS/load, all raw
rows and retained bytes; show actual pair counts and rebuild counts. Select
durations via untimed pilot before freezing the campaign. No unrelated heavy
team jobs while timing and do not stop unrelated user processes.

Production promotion requires at least5% median full-step improvement in a
declared case, favorable direction in at least4/5 matched pairs, reproduced in
a separate five-pair campaign, passing numerical/invariant gates and an explicit
memory/regression assessment. Failure to meet this bar means keep production
unchanged and publish the reproducible rejected experiment. Any accepted
candidate must integrate behind explicit experimental selection with unchanged
default and reject unsupported inputs; selection/API details get a focused
spec amendment before implementation. No engine/chemistry parity claim.

Rust floor1.90; no mandatory dependency; branches ricardo/; literal gitmoji +
imperative Action commits. Agent shell commands start rtk, product subprocesses
do not depend on RTK. Preserve earlier evidence, no force-push or extra merges.

# Focused force handoff

Review risks for Task 2: Ewald energy() currently calls exclusion_correction with an empty forces slice although correction writes into it. Cover energy-only special pairs so this cannot panic. K-vector cache key currently bounding edges only; same-size triclinic basis changes require invalidation. Tabulated lower-range extrapolation must retain U/F consistency; flat U with nonzero F is another derivative defect. Keep empty pair-scale fast path cheap; indexed lookup on every ordinary pair can erase LJ gains.

Special-pair long-range conventions should be documented and tested even when a special pair lies outside the real cutoff. Compare the implemented total energy correction to subtraction of (1-s)*q_i*q_j/r, or clearly account for the chosen finite erfc truncation. Do not infer chemical correctness from only distances inside rc.

CI runs warnings-as-errors, fmt and clippy default/GPU/MPI. Modified source and tests should not leave ignored Result warnings or inconsistent docs.

Ewald cache concurrency: ensure_kvec_cache currently unlocks before compute reacquires the mutex. Another simulation sharing one ForceModel can replace the box-keyed vectors between ensure/use. Return an immutable Arc snapshot for the requested full basis or keep ensure+use in one critical section; test concurrent different boxes against separate-model references. Per-instance SPME workspace must likewise not let cached box influence data cross simulations.

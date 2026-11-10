# B05 frozen SPME campaign results

The single frozen campaign on 2026-11-10 produced **no selected configurations,
no confirmed target IDs, and no demonstrated equal-accuracy speedup**. Two of
ten family/thread cases completed the measurement tooling; eight stopped at a
strict timing guard or projected work admission. B05 therefore has partial
measured coverage, not full campaign acceptance.

The [protocol](B05-SPME-PROTOCOL.md) was unchanged. Execution used commit
`9a25e2914beb0d477a8f296aec15299aa61e0700`, Apple M3 Pro, macOS Darwin 25.6.0,
Rust/Cargo 1.97.1, default features and the release example. The launcher ran
once, sequentially, from 10:52:04 to 11:28:37 UTC, with the six-level
`--reference-extended` preset and threads `1 4`. Independent inspection confirmed
all ten reports completed with successful report children, monitors and wrappers.
The historical launcher exit 0 alone did not enforce report completion; it also
did not mean the scientific targets passed.

These are historical measurements at `9a25e29`, documented at `a17cf2f`. The
subsequent final review fixes enforce and retain the cubic replication comparison,
aggregate preparation/report completion into `outcome.json` with a nonzero launcher
exit on operational failure, and correct new reports' elapsed provenance. The
registered execution source identity has changed. Frozen receipts and the retained
binary remain immutable evidence for the measured source only; any execution of
the fixed source requires fresh preparation and registration. No new campaign or
stronger scientific acceptance is claimed by these engineering fixes.

| Family | Threads | Tooling outcome | Target `1e-3` | Target `1e-4` |
| --- | ---: | --- | --- | --- |
| cubic64 | 1 | Failed: timing guard exit 87 | Pending | Pending |
| ortho64 | 1 | FullScopeMeasured | TargetNotReached | TargetNotReached |
| special3 | 1 | Failed: timing guard exit 87 | Pending | Pending |
| cubic216 | 1 | Failed: projected work BudgetExceeded | Pending | Pending |
| ortho216 | 1 | Failed: projected work BudgetExceeded | Pending | Pending |
| cubic64 | 4 | Failed: timing guard exit 87 | Pending | Pending |
| ortho64 | 4 | FullScopeMeasured | TargetNotReached | TargetNotReached |
| special3 | 4 | Failed: timing guard exit 87 | Pending | Pending |
| cubic216 | 4 | Failed: projected work BudgetExceeded | Pending | Pending |
| ortho216 | 4 | Failed: projected work BudgetExceeded | Pending | Pending |

All ten physical verifications passed. Special3 retained 1,509 passing controls
at each thread count; cubic64/216 retained their own converged crystal and
Madelung sign/normalization checks. The historical cubic replication comparison
was neither enforced nor separately retained by the CLI. Its retained per-pair
energies, `-1.5886950860301647` and `-1.588695086030165` at both thread counts,
differ by `2.220446049250313e-16`, below `1e-9/2.2`; the later review confirmed
this from the immutable individual checks. Completed selection references used six
levels, two cross-kappa checks and all seven translations. The independent
oracle was all-image real-space/full-positive-and-negative-k Ewald.

At each thread count, the 64-ion cases retained 3,556 selection error rows and
special3 retained 7,112, including baseline126. Cubic64 had one feasible
`1e-3` candidate (104); special3 had nine. Their `alloc_zeroed` callbacks during
guarded primary timing prevented complete costs, frozen choices and held-out
target confirmation. These operational failures do not establish target
unattainability. Neither family had a feasible `1e-4` candidate in selection,
but its interrupted decision state remains Pending. Both 216-ion cases were
rejected before candidate evaluation by projected complete-family work;
their scientific targets were not evaluated. No roster, cap or grid was reduced.

Ortho64 had no feasible candidate for either target. Its smallest grid-wide
worst relative force error was `0.0010978665855974885` (ID104), already above
`1e-3`. Baseline126 retained 28 selection and 21 held-out rows separately;
its relative RMS range was `0.765564–1.686782` on selection and
`1.046773–1.571038` on held-out snapshots. It was not target-accurate.
Each ortho64 baseline completed ordinary cold/rebuild/five-warmup/5×20-call
timing, 42 functional controls, and a separate 106-call, 12-field phase profile
on the first unshifted Selection snapshot. Ordinary median costs were
23,783.35 ns/call (one thread) and 23,958.35 ns/call (four threads).
These representative baseline diagnostics establish no target speedup.

All actual Rust children had separately recorded monitor/time-wrapper exits 0;
the Rust sweep exits were 87, 4 or 2 as described above. Numeric macOS RSS was
present for every child; the largest was 85,393,408 bytes. Requested Rust heap
stayed below 256 MiB and RSS below 512 MiB. No active-family timeout occurred;
the largest authoritative elapsed value was 375.680 s against 1,800 s.
Projected work failures remain failures despite these unused time/memory budgets.

The historical reports retain their elapsed-report limitation: ortho64 canonical
`provenance.family_elapsed_ns` stops at confirmation (311.970/311.375 s).
Final internal receipts report 326.451/325.703 s; the authoritative parent sums
of verify/sweep/confirm/profile are 332.070/331.268 s at one/four threads.
Use the conservative maximum of parent active sum and final internal elapsed;
the historical canonical field alone omits later profiling and parent overhead.
Newly generated reports from the fixed source use that conservative maximum and
retain both components in diagnostics; no frozen report has been rewritten.

Raw records are retained outside Git in
[`task8-frozen-20261110-v1`](../../target/b05-evidence/task8-frozen-20261110-v1/).
The [freeze manifest](../../target/b05-evidence/task8-audit-20261110/frozen-evidence-manifest.json)
binds 752 files/469,147,052 bytes, including
[source/binary/environment/test receipts](../../target/b05-evidence/task8-preflight-20261110/).
The [read-only audit](../../target/b05-evidence/task8-audit-20261110/raw-audit.json)
records every case's metric envelopes, raw translation ranges, reference deltas,
partial ordinary stages, phase fields, work, heap, RSS and remaining budget.

| Identity | SHA256 |
| --- | --- |
| Frozen evidence manifest | `0af7f838e2a991d51a64d9816432c90dd1dd94627875a39cd300218dbf27aee3` |
| Launcher registration | `cc2283587f0f8df35ce040f617ffd23b13a68d4bac2b6ba518bb2d106f1aa222` |
| Retained release executable | `e0e64e2fdaf33980852de19f186607c7a114a14ee6a108412ecb5eaf6dd9a59f` |
| Frozen protocol specification | `25c6e3a8af6c9520330295a4d276a20cdfd3fc6adb7b69a15489ed22a78a1f3f` |
| Derived raw audit | `7e01744b65896d6c070ccc3882881947a839d0069c712b234ca177d6aca629e8` |

Final pre-measurement checks passed: 513 default-workspace tests, seven ignored
harnesses/doc examples, 13 launcher tests, and the release build. The workspace
used `RAYON_NUM_THREADS=1` with default Rust harness concurrency; launcher tests
used the real debug report helper. Optional ASE/MACE, GPU and multirank MPI
coverage is not established by this run. Existing explicit parameters,
analytic order-4 differentiation and defaults are unchanged. Results cover the
frozen reduced-unit orthorhombic snapshots only: no universal error bound,
material/trajectory guarantee, or full stress tensor certification follows.

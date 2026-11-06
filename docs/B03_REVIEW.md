# B03 review closure and delivery

B03 is delivered as a reproducible rejected scalar-layout experiment in
[PR #57](https://github.com/ricardoleal20/modena/pull/57), open against B02's
`ricardo/b02-pair-scratch`. Assignee: ricardoleal20; labels: enhancement, rust,
python. B09 #55 and B02 #56 remain open. No additional merge was performed.

## Review and validation

All three task gates approved. Task 2 closed its three original validator/test
findings and the subsequent CI prerequisite finding in two focused rounds.
Task 3 independently verified all 115 new evidence files and the all-case
receipts. The final whole-branch reviewer approved all 16 changed files across
`604039e1ecde99b9594acfe9d86cb90bf1fae9ab` through
`63e5c466dc2b185609578a199532e894e0a22388`, with no Critical or Important finding.
Later delivery metadata changes only documentation, preserving that reviewed
implementation. This is publication approval, without merge authorization.

Fresh checks pass: 47 Python benchmark tests, seven Rust integration tests,
four example tests, scoped formatting/Clippy with `-D warnings`, and Rust1.90
checks. Production core remains tree
`39b46c7496c6820219056fe38a821378a32633e8`. No existing production kernel,
storage, threshold or default was changed.

Both protocols retain 1,344 invocations and 96 rejected comparisons; 32 pilot
invocations remain separate. Full saved-state primary validation, snapshot-free
limitations, original measured ELF5001bb6 versus fresh contract ELF822b000,
source epochs and exact commands are recorded in
[B03 verification](B03_VERIFICATION.md) and the [all-case guide](B03_PACKED_LAYOUT.md).
This rejects these scalar implementations in the stated synthetic domain;
unimplemented SIMD/parallel layouts and chemical accuracy remain untested.

## Deferred minor triage

The stdout/snapshot subset issue is closed by strict typed shared-value checks,
documentation and the mandatory actual Rust two-mode regression.
Unbuffered post-timer snapshot output remains a nonblocking Minor. The final
reviewer explicitly requires no correction before this PR. Before another
collection, a buffered writer should propagate explicit flush errors, run
focused correctness checks and capture fresh executable provenance. Historical
source, ELF and measurements must remain unchanged. The syscall explanation
is an inference; no syscall trace was measured.

## Rulings I made

| Decision and reason | Cost if wrong |
| --- | --- |
| Test benchmark-only packed views before production promotion; the backlog requires measured benefit. | Experimental work may be discarded and production promotion delayed. |
| Reuse the isolated B02 checkout and stack B03 on open B02; preserve predecessor work without an unauthorized merge. | PR base synchronization after predecessor integration. The managed checkout remains available for feedback. |
| Compare scalar packed views with equivalent interleaved-serial execution and keep production-auto separate. | These data cannot decide unimplemented parallel/SIMD benefit; another matched implementation/campaign may be needed. |
| Add a complete separate snapshot-free confirmation matrix because retained primary JSON affects cache/RSS. | About8–12 extra minutes and analysis complexity; the secondary lacks full vectors/state SHA and cannot replace the primary. |
| Defer BufWriter correction to a future collection; the current writer is outside VV and its Minor is nonblocking. | Future collection retains avoidable IO overhead until a focused writer/provenance change. |

## Actual remote CI

Observed Community run
[42205811533](https://github.com/ricardoleal20/modena/actions/runs/42205811533),
[job127852023197](https://github.com/ricardoleal20/modena/actions/runs/42205811533/job/127852023197),
at reviewed source63e5c46 failed with an empty steps list. Its annotation states:
“The job was not started because your account is locked due to a billing issue.”
The reported title-check failure therefore is not evidence of an executed title
validation. The main CI workflow filters PR base `main`; this stacked feature-base
PR did not trigger that workflow. Local validation is not remote CI success.
No protections or check statuses were changed.

## Workflow preservation and remaining work

The controller copied and independently SHA256-verified this plan's 46 ignored
workflow files /1,059,282 bytes, including future-case preparation, into
`target/b03-evidence/review-archive-pr57`. Manifest SHA256:
`b138b99f3246a55d51f0c5c9a211140160efbbcd18b11686931e1e6bcdc50427`.
Its separate receipt records the actual verification. The managed worktree and
all original campaign artifacts remain available; only this plan's disposable
SDD directory is eligible for cleanup. Those local files are not in the PR;
a fresh checkout must reproduce the documented commands.

The eight-case goal remains active. B04 is independently active in its own
worktree; B05/B07/B08/B10/B11/B12 still require their own implementation,
scientific/runtime evidence, review and PRs. B10's selected material is FCC
aluminium. See the [execution record](CORE_RESEARCH_EXECUTION.md).

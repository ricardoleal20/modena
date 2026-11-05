# B02 review closure and delivery

[PR #56](https://github.com/ricardoleal20/modena/pull/56) publishes
`ricardo/b02-pair-scratch` against `ricardo/b09-ase-sessions`, the open
[B09 PR #55](https://github.com/ricardoleal20/modena/pull/55). Both PRs remain
open; no B09/B02 merge or force-push was performed. The reviewed B02 head is
`b99db28cd501e54af2ff228aa51e9c3fb0943fa4`; subsequent shipping changes are
documentation only. After a B09 squash, B02 needs base synchronization before
merging into main.

Separate read-only task reviews approved the implementation, measured
acceptance and fresh validation record. The whole-branch review found no
Critical/Important defect or design fault, and one Minor benchmark-validation
gap. Its single grouped fix wave and scoped rereview closed that gap with no
new breakage or residual findings. These are AI-assisted engineering reviews,
not human peer review or material/ensemble validation.

| Closed finding | Correction and evidence |
| --- | --- |
| Product benchmark provenance invoked RTK | Direct system-tool argv; subprocess regression with RTK absent. RTK remains only the agent shell convention. |
| Numeric equality was insufficient for a bit-exact one-thread gate | Rust bit comparisons and actual checkpoint/thermo hashes, with signed-zero regression. |
| Scratch-order comment was inaccurate | Documentation corrected without changing kernel bytes or timing samples. |
| Benchmark accepted matching malformed saved continuation state | Actual saved identity buffers, typed DOF and explicit null controller states are checked independently of the declared signature. Complete v2 fixture passes the actual Rust reader; shared bad-mass and missing-field regressions fail. |

The final fix is `940c95b`, followed by documentation `b99db28`. Its eight
focused driver tests and all **37 benchmark-tool tests** passed. The stronger
validator rechecked all **380** retained subprocess rows and **48** exact
one-thread checkpoint/thermo pairs; all **1,943** original artifact hashes
remain unchanged. All 48 production-source/manifest hashes match the prior
whole-branch review base `8d20fed`; no kernel change, rebuild or retiming was
needed for this tooling fix.

The earlier fresh broad campaign remains attached to source `360defb`:
**357 Rust passes**, zero failures, one existing ignored doc example;
**35 actual-ASE Python passes**, five worker tests and **67 scientific
observations**. Formatting, default/GPU all-target Clippy and Rust 1.90
default/GPU/standalone-harness checks passed. That campaign's 34 tooling tests
are historical; 37 is the final tooling count. GPU checks are compile/lint
only; no GPU/MPI runtime or other-engine campaign is claimed.

[Verification](B02_VERIFICATION.md) preserves commands, epochs and hashes.
[Performance guide](B02_PAIR_SCRATCH.md) retains every timing case, including
the rejected dense variant and unfavorable rows. In the declared N5000/4T
case, median full-CLI time fell from 1.355683 to 1.248346 seconds, about 7.9%;
all five pairs and a separate three-pair pilot favored the candidate. This
does not establish uniform large/default-pool speedup. O(TN) scratch remains
retained between calls, and process RSS can increase even when Rust-requested
peak bytes fall. Local raw evidence is ignored under `target/b02-evidence/`
and is unavailable in a fresh clone.

## Decisions and costs

| Decision | Reason | Cost if the choice needs revision |
| --- | --- | --- |
| Stack B02 on open B09 | Keep both requested features usable without an unauthorized merge | PR dependency and base synchronization after B09 squash; no force-push assumed. |
| Retain the parallel-block variant for the repeated N5000/4T benefit | Meets the declared measured-case acceptance; the initial dense variant failed | Other workloads/load regimes may be neutral or slower; retained memory and RSS costs remain. New speed claims require measurement. |

The four B09 decisions and their costs remain in
[B09's closure record](B09_REVIEW.md).

## Remote verification

At initial publication, Community run
[42090260274](https://github.com/ricardoleal20/modena/actions/runs/42090260274),
job `127466127520`, executed **no steps**. Its annotation says the account
billing lock prevented the job from starting. This is not evidence that the
title checker or tests ran and failed.

The main CI workflow filters pull requests to base `main`; it does not run
for this stacked feature-branch base. Neither local validation nor an absent
main-CI check establishes remote CI success. No protection settings or check
statuses were changed. The PR is assigned to `ricardoleal20` and labeled
`enhancement`, `rust` and `python`.

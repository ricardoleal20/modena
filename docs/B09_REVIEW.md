# B09 review closure and delivery

[PR #55](https://github.com/ricardoleal20/modena/pull/55) publishes
`ricardo/b09-ase-sessions` against main. It is open; no B09 merge was performed.
The functional fix epoch is `f14d6efd4a2a6e84e658a2f60a39805cf783ef55`;
`3e7278a` adds its measured verification record. Subsequent shipping edits
are documentation only.

Separate read-only task reviews approved the persistent provider, reverse
calculator and benchmark. The whole-branch review identified three Important
and two Minor findings; its scoped rereview closed all five with no new
Critical/Important defect. These are AI-assisted engineering reviews, not human
peer review or material/ensemble validation.

| Closed final finding | Correction |
| --- | --- |
| ASE's default tolerance reused a singular changed snapshot | Exact state comparison; threshold-crossing failure clears cached results. |
| Large finite periodic translations lost their remainder | Remainder-based snapshot normalization, tested for both signs. |
| Dense snapshots could allocate billions of neighbor entries | Snapshot-only streamed builder limits retained pairs to 1,000,000; force pass consumes the bounded list, including zero skin. Normal MD remains unrestricted. |
| Benchmark accepted identical missing physical data | Required physical signature and eight thermo columns are validated before parity. |
| Empty no-calculation ASE probe raised instead of returning None | Supported cache probes preserve ASE's absent-result semantics. |

Earlier persistent-session findings (future replies, pre-limit materialization,
idle fatal cleanup, retained event senders) and the fixture synchronization race
were also corrected and reviewed before the final wave.

Fresh final validation passed 355 Rust tests with zero failures and one
intentional ignored doc example, Python 35 actual-ASE tests (26 plus nine skips
without ASE), five worker tests, 29 tooling tests and 67 scientific observations.
Formatting, default/GPU Clippy and Rust 1.90 default/GPU/required-harness checks
passed. GPU checks are compile/lint only; MPI/GPU runtime campaigns were not run.
[Verification](B09_VERIFICATION.md) preserves all source epochs, failed initial
checks, fixes, exact commands, hashes, timings and untracked local artifact paths.

## Decisions and costs

| Decision | Reason | Cost if the choice needs revision |
| --- | --- | --- |
| Keep `one_shot` as default; persistence is explicit | Preserve B01 configuration behavior and a controlled comparison | Users must opt in for startup savings; a later default change needs migration. |
| Start the inverse adapter with single-species native LJ E/F | Existing scalar virial cannot supply ASE tensor stress | Other models, cells and stress require a separate extension. |
| Synchronize the inherited standalone-harness lock | Required Rust 1.90 CI check lacked B01's tempfile dependency | Extra lock-only review surface; nine packages were added without changing existing pins. |
| Add a documented snapshot pair-storage budget | Atom/grid caps admitted about 80 GB of dense pair entries | Some dense inputs now reject; the bounded serial search has its own cost. |

## Remote verification

At publication, Actions run `42073080899`, job `127405936991`, executed no steps.
Its annotation says the account billing lock prevented the job from starting.
Local verification does not establish remote CI success. No protection settings
or successful check statuses were changed. The next requested work is B02.

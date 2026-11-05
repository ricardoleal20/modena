# B01 review closure and decisions

B01 implements the first concrete external total-energy provider: pinned ASE
EMT, with explicit atomic numbers, full 3D periodicity, real-unit boundaries,
fixed-cell unconstrained NVE/NVT, and a fresh calculator process per snapshot.
This record closes local engineering review; it does not certify material
accuracy, ensemble convergence, throughput or parity with another engine.

Logical base: `0af8ec9c72b2b64b42bf41ef5b26d5b638c5a2f4`, the separately reviewed
B06 PR head. Source reviewed through `05b56429e169adf7db3723d59784e1dfb1b75e8c`.
The implementation, runner/restart, interface and documentation tasks each
received a separate read-only review. A subsequent whole-branch review found
four issues. Corrections `65ec427` and `119a221`, then narrow closure `05b5642`,
were checked by scoped reviews against those findings.

## Closed findings

| Finding | Implemented correction and evidence |
| --- | --- |
| Real pressure used a bar-like factor labeled atm | Shared SI factor corrected to `68568.42296625092`. Independent conversion, MC/Berendsen normalization and actual cold-Cu CLI versus direct ASE/SI regressions pass. Absolute cold-Cu pressure error is `2.61e-9 atm`. |
| Failed integration could overwrite a good checkpoint or run again | Failed instances are terminal; writes reject them before touching a file. Completed progress advances only on a successful step. The prior good checkpoint remains byte-for-byte intact and resumable. |
| Cancellation could launch additional provider calls during setup | Typed cancellation and checks before/between setup evaluations and before steps. Synchronized subprocess latches verify invocation counts; clean zero-step cancellation remains reusable. |
| Relative interpreter names differed across interfaces | Shared lexical base-directory resolution preserves executable symlinks and bare PATH lookup. A dotted PATH executable wins over a same-named input-directory decoy; non-UTF-8 resolved paths fail explicitly. |

All four findings and the related UTF-8 issue were closed at `05b5642`, with no
new fix breakage reported. The earlier async-test wording concern is closed:
its evidence is lifecycle/final metadata, not a guaranteed nonterminal snapshot.
The known upstream ASE 3.29/NumPy 2.5 deprecation warning remains visible and
nonblocking. These were AI-assisted read-only reviews, not independent human
peer review or a general scientific validation campaign.

## Verification boundaries

The complete numerical fix epoch passed 330 reported Rust functions, one
existing ignored doc example and six explicitly skipped optional references.
Separate actual-dependency runs passed 24 core, 53 CLI/MCP, 25 wrapper and
three worker tests, plus 22 tooling tests and 67 scientific observations.
Default/GPU Clippy, Rust 1.90 checks and formatting passed.

The last interpreter-only change subsequently passed its two helper regressions,
32 actual core/CLI/MCP checks, four Python external checks and compile/lint/format
checks. Its PATH/decoy test failed before the correction and passed afterward.
Unchanged numerical campaigns were not rerun or relabeled as that later epoch.
[Verification](B01_VERIFICATION.md) records exact commands, hashes, local
artifacts and optional skips; [the guide](B01_EXTERNAL_PROVIDER.md) records
supported inputs, installation and migration.

## Decisions made during implementation

1. **Introduce `System.atomic_numbers` with the provider unit.** This lets the
   real ForceModel consume chemical identity before runner integration, without
   a test-only interface. If that boundary changes, the public field and its
   consumers require coordinated revision.
2. **Use ASE EMT / periodic Cu as the first reference.** It is reproducible
   without downloading model weights. A different requested model requires an
   adapter/task revision; the identity/unit/parity infrastructure remains useful.
3. **Preserve existing alloy mass defaults.** The native non-optional field
   cannot distinguish omitted mass from explicit default after parsing. Optional
   FCC/XYZ/rocksalt mass fields require explicit values for external runs. A
   stricter alloy declaration can be added later without changing classical
   defaults; no mass is inferred from atomic numbers.
4. **Preserve the by-value final event payload.** A localized large-enum lint
   allowance retains its public shape while carrying owned provider identity.
   Boxing it instead would require an explicit API change and review.
5. **Correct the existing real-pressure conversion.** This is an explicit
   exception to blanket native-output preservation: reduced mechanics and real
   fixed-cell dynamics stay unchanged, while real pressure output and pressure
   coupling normalization change. Only exact pressure-factor-only fixed-cell
   v2 signatures may migrate. Old pressure-coupled runs require an explicit new
   protocol; unrelated identity differences reject. Tightening migration to
   reject all legacy checkpoints would sacrifice safe fixed-cell convenience.
6. **Finish a narrowly scoped residual path fix.** The first fix review exposed
   a remaining wrong-executable selection case. Rather than park a confirmed
   Important defect at a workflow iteration cap, implementation added the small
   PATH/decoy and UTF-8 correction and one scoped closure review. Its cost was
   one additional focused patch/review, without reopening closed physics work.

GitHub shipping is separate from local review closure. B06 was subsequently
merged as PR #53 with the owner's explicit administrator exception; B01 was
rebased without content changes onto that main. See
[shipping provenance](B01_VERIFICATION.md#shipping-provenance). Local checks
and the administrator exception do not establish remote CI success.

# B04 implementation and verification record

B04 numerical work is starting in the isolated `ricardo/b04-periodic-neighbors`
branch based on reviewed B03 source `0c54db355b830c49e687e97e2e9aa8502bbaa5d1`.
The core tree is `39b46c7496c6820219056fe38a821378a32633e8`, identical to reviewed
B03 prototype66af0e4. No core/CLI files changed between those revisions.
B03 primary/secondary timing and direct audit finished before this work; its
CPU lease was explicitly released. B03 reports/review remain separate closure
work. Before B04 publication, synchronize its unpublished base with final
reviewed B03 without force-push or an unauthorized merge.

Controller preflight confirmed native managed worktree creation, branch prefix
`ricardo/`, source-tree equivalence and the actual Vec3 dot/wrap operation
orders used by the proof. No B04 implementation tests have run yet.

| Gate | Status | Required evidence |
| --- | --- | --- |
| G1 | pending | Fixed-cube/resource tests and official ASE exact tuples |
| G2 | pending | Actual initialized official MACE graph/E/F parity |
| G3 | pending | Scalar/stencil pairs, E/F/W, skin/bounded regressions |
| G4 | pending | Paired complete-step benefit and fallback overhead |

Source epochs, commands, numerical tolerances, package/model hashes and measured
results will be appended when actually executed. Prepared source contracts do
not count as acceptance evidence.

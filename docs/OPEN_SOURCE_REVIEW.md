# Open source preparation review — historical snapshot

Snapshot date: 2026-10-23. Baseline: `a3e475c` on `main`, MODENA 0.1.1.
This record captures local preparation evidence gathered before the changes
were committed or pushed. At that point no package, release, or preparation PR
had been published, and GitHub had no open PR when checked. Pre-existing local
audit material was left untouched. This is a historical preparation snapshot,
not a statement of the current branch or publication status. For hosted
verification after push, consult the associated pull request's Checks tab.

## Repository and distribution

The existing MIT license was retained. Exact copies were added to the three
distributable Rust crates and the Python package. Cargo package manifests
included their README and license notices. Python wheel and source distribution
declared SPDX MIT, included the complete notice, and provided project metadata.
No new runtime dependency was added.

Contribution and commit rules, three issue forms, a PR template, CODEOWNERS,
conduct and security policies, governance, support, and an editor configuration
were prepared. CI used read-only workflow permissions, timeouts, locked Cargo
commands, and concurrency cancellation. New checks were prepared for PR titles
and Python distributions; Dependabot configuration covered all four dependency
locations.

[DESIGN.md](../DESIGN.md) separates a proposed visual system from the static
prototype and records color contrast, typography, units, scientific
communication, accessibility, and motion guidance. Logo concepts were explored
and inspected locally during preparation; that artwork and its previews were
deferred and excluded from the contribution.

## Code review and fixes

Targeted review covered core setup and validation, structure allocation, runner,
XYZ/restart/log I/O, CLI entry points, MCP protocol and run management, Python
configuration and subprocess execution, package manifests, and CI.

The following reproducible MCP defects were fixed and independently reviewed:

| Finding | Result |
| --- | --- |
| Snake-case protocol response keys | Standard `protocolVersion`, `serverInfo`, `inputSchema`, and `isError` wire keys. |
| Overflowing generated structure counts or nonfinite boxes | Recoverable input errors before lattice allocation; subsequent stdio requests still succeed. |
| Concurrent admission counted only running jobs | Queued and running slots reserved atomically under one registry lock. |
| Invalid or oversized sweeps started inaccessible partial work | Every configuration is validated and the complete batch reserves capacity before launching. |
| Cost estimates overflowed or hid unreadable/empty XYZ input | Checked counts and geometry, explicit errors through `Result`, and finite budgets. |

Regression tests were demonstrated failing before each fix and passing after
it. Exact protocol keys were checked through the server's stdio session. New
admission tests exercised simultaneous callers, queued capacity, sweep rejection,
terminal-slot reuse, invalid boxes, and estimate errors.

A redacted scan of 629 tracked UTF-8 files found no recognizable private-key,
AWS, GitHub, or OpenAI credential candidates. That scan excluded the user's
untracked audit folder. No vendored source or copied third-party headers were
identified by targeted scans. These observations are not a complete provenance,
legal, dependency-advisory, or Git-history audit.

## Local verification recorded at the snapshot

| Check | Result |
| --- | --- |
| `cargo test --locked --workspace` | 284 passed, one ignored documentation example, 35 suites. |
| `cargo fmt --all -- --check` | Passed. |
| `cargo clippy --locked --workspace --all-targets -- -D warnings` | Passed. |
| `cargo build --locked --release --workspace` | Passed. |
| Rust 1.90 workspace/all-target check, default and `modena-core/gpu` | Both passed with locked, offline resolution. |
| `python3 -m unittest discover -s python/modena/tests -v` | 18 passed, with the rebuilt release CLI available. |
| `python3 -m unittest discover -s bench/tests -v` | 10 passed outside the macOS resource-timer sandbox restriction. |
| CLI quickstart `check` and temporary source installation | Passed. |
| Python wheel and source distribution | Built; exact MIT bytes and metadata checked; isolated installed-wheel import/config rendering passed. |
| Cargo package file lists | LICENSE and README present in core, CLI, and MCP packages. |
| `actionlint` 1.7.12 | Both workflows passed; official download checksum verified. |
| PR title matcher | Six accepted and seven rejected examples. |
| Structured files | Seven YAML and six TOML files parsed; six deferred local SVG files were also inspected. |
| Local documentation links | 60 relative links resolved. |
| `git diff --check` | Passed. |

The full Rust suite ran with stable Rust/Cargo 1.97.1 on macOS arm64. Python
regressions ran with Python 3.14.5; distribution installation also used an
isolated Python 3.12 environment. Temporary build artifacts and installations
were placed outside the repository. No physical-model implementation changed.

## GitHub configuration

Read back and confirmed after configuration:

- Public repository, MIT, `main` default branch, issues enabled.
- Private vulnerability reporting enabled.
- Dependabot alerts and security fixes enabled.
- Secret scanning and secret push protection enabled.
- `main` protection requires an up-to-date PR, resolved conversations, and eight
  real CI checks from successful runs; force pushes and branch deletion blocked.
  Required approval count is zero for this single-maintainer project; owner
  administrator bypass is retained.

At snapshot time, the existing successful hosted checks belonged to the
baseline; the changed workflows had not been pushed or run on GitHub. This
snapshot therefore establishes local evidence only. The associated pull
request's Checks tab is the source for subsequent hosted results. Issue forms
and version-update configuration take effect after integration into the default
branch. Add community/distribution check names to branch protection only after
confirming their actual names in successful hosted runs.

## Scope limits

This preparation covered contribution, distribution, design, and release processes.
It does not establish complete scientific correctness, production suitability,
dependency license provenance, or registry publication. GPU compilation was
checked locally; GPU hardware validation, multirank MPI, and Linux execution
were not rerun in this local preparation. The CI configuration retained those paths.
See [VALIDATED_SCOPE.md](VALIDATED_SCOPE.md) for scientific evidence and
[MAINTAINING.md](MAINTAINING.md) for future releases.

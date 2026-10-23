# Maintaining MODENA

The maintainer and decision process are documented in
[GOVERNANCE.md](../GOVERNANCE.md). This guide covers routine repository upkeep
and release preparation. It does not publish packages or promise release dates.

## Review and merge

Keep changes focused and follow the [commit and PR title rules](../COMMIT_RULES.md).
The community workflow checks new pull request titles. Prefer squash merging
with the reviewed PR title; individual work-in-progress commits need not be
rewritten. Require a passing CI run for the exact commit being merged. Optional
GPU and MPI checks must report what ran and what skipped.

For numerical changes, review the independent reference, units, tolerances,
sampling limits, and affected validation claims. Update
[VALIDATED_SCOPE.md](VALIDATED_SCOPE.md) when the evidence or limitations change.
Keep historical results intact and attach fresh provenance to new measurements.

Do not combine a release with unrelated refactoring. Do not publish a version
whose documented installation path or package contents have not been checked.

## Repository settings

The public repository is <https://github.com/ricardoleal20/modena>, with `main`
as the default branch. The source uses MIT. Issue
templates and contribution files become active on GitHub after they reach the
default branch.

The following security settings were enabled and checked on 2026-10-23:

- Private vulnerability reporting, linked from [SECURITY.md](../SECURITY.md).
- Dependabot vulnerability alerts and automated security fixes.
- Secret scanning and secret push protection.

[Dependabot configuration](../.github/dependabot.yml) adds weekly version-update
checks for the root Cargo workspace, the standalone validation harness, Python,
and GitHub Actions after it reaches `main`. Review dependency updates with the
same numerical and compatibility checks as other changes. Compatible Cargo
updates are grouped per workspace, and Actions updates share one weekly PR.
Major Cargo updates remain separate. The Python policy avoids raising a minimum
that already admits the new version unless a compatibility or security reason
requires it. The numeric `dtolnay/rust-toolchain` ref is excluded from automatic
version bumps because it selects the declared Rust compiler floor; change that
floor deliberately in both Cargo metadata and CI.

Keep both Cargo lockfiles coherent when changing shared workspace dependencies.
The compiler-floor CI job builds the standalone validation harness as well as
the root workspace. Python distribution jobs verify Python 3.11 and 3.12,
including installation and import of the optional MCP tooling.

Branch protection was also enabled for `main`: pull requests, an up-to-date
branch, resolved conversations, and the eight existing CI checks are required;
force pushes and branch deletion are blocked. Check names were taken from
successful GitHub runs. The owner retains administrator bypass, and the rule
does not require a second person's approval while the project has one maintainer.
The required checks are format, Clippy, Linux tests, macOS tests, Rust 1.90,
MPI, GPU, and benchmark compilation. Before adding community and Python
distribution checks to the rule, confirm their actual names in successful
hosted runs. Add review requirements when another maintainer can satisfy
them. Invented check names can prevent every merge.

## Release checklist

1. Review open regressions, the latest validation scope, dependency alerts, and
   the diff since the previous release. Record accepted limitations.
2. Update the workspace version in `Cargo.toml`, the Python version in
   `python/pyproject.toml` and `python/modena/__init__.py`, `CITATION.cff`, and
   `CHANGELOG.md`. Internal Cargo dependency versions and lockfiles must agree.
   Use the actual release date; do not change historical release records.
3. Confirm that `LICENSE`, each distributable Rust crate's `LICENSE`, and
   `python/LICENSE` are byte-identical. Preserve third-party notices if future
   contributions include third-party material.
4. Run the baseline checks from the repository root:

   ```sh
   cargo fmt --all -- --check
   cargo clippy --locked --workspace --all-targets -- -D warnings
   cargo test --locked --workspace
   cargo build --locked --release --workspace
   python3 -m unittest discover -s python/modena/tests -v
   python3 -m unittest discover -s bench/tests -v
   cargo run --locked --release --bin modena -- check -i examples/argon_nve.toml
   ```

5. Run the compiler-floor, GPU, MPI, and numerical smoke checks defined in
   [CI](../.github/workflows/ci.yml). Confirm two real MPI ranks; a single-rank
   skip does not establish MPI behavior. A successful GPU compile or adapter
   skip is not evidence of hardware accuracy or speed.
6. Build the Python sdist and wheel from `python/` with `python3 -m build` in an
   isolated environment containing the build frontend. Inspect both archives
   for `LICENSE` and metadata, install the wheel into a fresh environment,
   and import `modena` from outside the checkout. Verify a CLI-backed example
   separately: the Python wheel does not contain the Rust executable.
7. Use `cargo package --list -p modena-core`, and the equivalent commands for
   `modena-cli` and `modena-mcp`, to inspect the package README and MIT notice.
   Package from a clean commit for a release. If publishing to crates.io,
   validate and publish `modena-core` before packages that depend on it.
   Registry publication is a separate maintainer action; source availability
   does not imply that a package name or version is available on a registry.
8. Record toolchain, platform, exact commit, commands, results, skips, and
   validation limitations in a versioned release note under `docs/releases/`.
   Confirm the GitHub CI run for that commit, then tag and publish the release
   deliberately. Mark early development releases as prereleases.

GitHub releases and registry uploads are separate steps. A successful local
check is not evidence that a package was uploaded or that hosted CI passed.

## Identity changes

Follow [DESIGN.md](../DESIGN.md). The visual system is proposed, and logo artwork
is deferred. Any future identity contribution should include editable sources,
accessible descriptions, and a preview covering light, dark, monochrome, and
small-size use. New images, fonts, datasets, and copied source need a known
redistribution license and any required notices before they are added.

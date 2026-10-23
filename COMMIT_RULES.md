# Commit and pull request titles

Use Conventional Commits for new commits and pull request titles:

```text
type(scope): short imperative description
```

The scope is optional. Use a lowercase type and a clear description of the
change. A body can explain the problem, approach, numerical evidence, and
tradeoffs. Link related issues in the body or footer.

| Type | Use |
| --- | --- |
| `feat` | New behavior or capability |
| `fix` | A bug or numerical correctness fix |
| `docs` | Documentation only |
| `test` | Tests or validation coverage |
| `perf` | A measured performance improvement |
| `refactor` | Restructuring without intended behavior changes |
| `build` | Dependencies or build configuration |
| `ci` | Continuous integration |
| `chore` | Repository maintenance |
| `revert` | Reverting an earlier change |

Scopes should identify the affected area, such as `core`, `cli`, `python`,
`mcp`, `gpu`, `mpi`, or `bench`.

```text
fix(core): preserve force cache after rejected volume moves
docs: explain numerical validation limits
test(python): cover CLI discovery with a custom target directory
```

Mark an incompatible API, input-format, or checkpoint change with `!` after
the type or scope, and explain the impact and migration in a `BREAKING CHANGE:`
footer:

```text
feat(cli)!: require an explicit unit system

BREAKING CHANGE: inputs must declare their unit system; update existing TOML files.
```

The example illustrates the format, not current MODENA behavior. Describe
changes to physical interpretation explicitly, even when input syntax stays
compatible. A correction to units or a Hamiltonian can invalidate earlier
scientific comparisons.

Keep commits focused and avoid mixing unrelated formatting with behavior changes.
These rules apply going forward. Do not rewrite historical commits to adopt
them. A maintainer may adjust the final title when merging; this policy does
not itself configure automated releases or commit enforcement.

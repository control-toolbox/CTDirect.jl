# CTDirect.jl — Agent Navigation Guide

Quick-reference for any agent working on this repository.

---

## Repository Layout

CTDirect is a single flat module — no submodule split, no `ext/`.

```text
src/           # DOCP construction (DOCP_*.jl), collocation.jl, direct_shooting.jl, ode/ schemes
test/          # Test suite: flat files (not test/suite/) — aqua.jl, problems/, gpu tests
docs/          # Documentation site (DocumenterVitepress)
.conventions/  # ct-conventions, cloned at the pinned tag (git-ignored)
```

---

## Conventions

The coding conventions live in
[ct-conventions](https://github.com/control-toolbox/ct-conventions), pinned below to a
SemVer release tag (`vX.Y.Z`). Read them from a local clone in `.conventions/`
(git-ignored, not a submodule). If it is missing, or `git -C .conventions describe --tags`
is not the pinned tag, (re)create it:

```bash
rm -rf .conventions
git -c advice.detachedHead=false clone --quiet --depth 1 --branch v0.7.0 https://github.com/control-toolbox/ct-conventions.git .conventions
```

**Pinned version: ct-conventions v0.7.0.** Cite rules by identifier (`MOD3`, `TST1`), and
state the version you applied in plans and reports.

Read in this order:

| File | Content |
| --- | --- |
| [`rules/general.md`](.conventions/rules/general.md) | Rules that hold everywhere — **always read** |
| [`README.md`](.conventions/README.md) | Index of the rule files, one per area (`MOD`, `TYP`, `DIM`, `EXC`, `PERF`, `DOC`, `TST`, `SITE`, `CHG`, `REL`) |
| [`method.md`](.conventions/method.md) | How to work: tests, coverage, docs build, output capture |
| [`plan-template.md`](.conventions/plan-template.md) | Template for a substantial piece of work |

How agents behave on these repositories is described in
[ct-agents](https://github.com/control-toolbox/ct-agents). The CI workflows and their
trigger labels are described in
[CTActions `docs/WORKFLOWS.md`](https://github.com/control-toolbox/CTActions/blob/main/docs/WORKFLOWS.md).

---

## Key Conventions

- **English everywhere** — code, docs, commits, issues, PRs (`GEN1`).
- **No top-level exports** — use `Package.symbol` everywhere (`MOD8`).
- **Qualified imports** — `using Pkg: Pkg`, never bare `using Pkg`; `import` is never used (`MOD3`, `MOD4`). The current source still uses `import` in places: the rule is the target for new and edited code, not a description of existing lines.
- **Fake types at module top-level** — never inside test functions (`TST3`).
- **Structured errors** — seven typed exceptions; pick by the choice rule (`EXC1`).
- **Type stability enforced** — hot paths `@inferred`-clean, verified with JET (`TYP8`, `PERF1`–`PERF5`).
- **1-D is a scalar** — a one-dimensional state/control/variable is a `Number` (`DIM1`).
- **Plan before code** — a local agent proposes a plan ([`plan-template.md`](.conventions/plan-template.md)) and waits for confirmation before working on code. The process for remote agents is still to be defined.
- **Docstrings last** — written once the implementation is stable (`DOC1`–`DOC8`).
- **The human merges** — an agent never merges a pull request, except a release PR in `autonomous` mode, and then only with `release.jl merge` (see below).

---

## Recurring tasks

- **Release.** When asked, in any language, to prepare a release, follow
  [`recipes/release.md`](https://github.com/control-toolbox/ct-agents/blob/main/recipes/release.md)
  of ct-agents (rules `REL1`–`REL7`). Read it with
  `gh api repos/control-toolbox/ct-agents/contents/recipes/release.md -H "Accept: application/vnd.github.raw"`.
  Stable or beta is the human's choice: if the request does not say, ask; never infer it.

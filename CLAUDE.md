# CTDirect.jl — Agent Navigation Guide

Quick-reference for any agent working on this repository.

---

## Repository Layout

CTDirect is a single flat module — no submodule split, no `ext/`.

```text
src/        # DOCP construction (DOCP_*.jl), collocation.jl, direct_shooting.jl, ode/ schemes
test/       # Test suite: flat files (not test/suite/) — aqua.jl, problems/, gpu tests
docs/       # Documentation site (DocumenterVitepress)
```

---

## Developer Resources

Design philosophy, operational rules, plan templates, and CI/CD conventions live in the
[control-toolbox Handbook](https://github.com/control-toolbox/Handbook):

| Topic | Link |
| --- | --- |
| Code philosophy (modules, types/traits, exceptions, docstrings, testing, docs) | [`PHILOSOPHY.md`](https://github.com/control-toolbox/Handbook/blob/main/PHILOSOPHY.md) |
| Operational rules (tests, coverage, docs, git) | [`RULES.md`](https://github.com/control-toolbox/Handbook/blob/main/RULES.md) |
| Plan template | [`PLAN.md`](https://github.com/control-toolbox/Handbook/blob/main/PLAN.md) |
| CI/CD workflows (centralized reusable workflows, label-gated triggers) | [`WORKFLOWS.md`](https://github.com/control-toolbox/Handbook/blob/main/WORKFLOWS.md) |

---

## Key Conventions

- **No top-level exports** — use `Package.symbol` everywhere.
- **Qualified imports** — `using Pkg: Pkg`, never bare `using Pkg`; `import` is never
  used. (The current source still uses `import` in places — treat the Handbook rule as
  the target for new/edited code, not a description of existing lines.)
- **Fake types at module top-level** — never inside test functions.
- **Structured errors** — seven typed exceptions under `CTException`; pick by the
  IncorrectArgument / PreconditionError / NotImplemented rule.
- **Type stability enforced** — hot paths must be `@inferred`-clean, verified with JET;
  setup-path dispatch is fine.
- **1-D is a scalar** — a one-dimensional state/control/variable is a `Number`, never a
  length-1 vector.
- **Plans before code** — write a plan and confirm with the user before touching files.
- **Docstrings last** — written only after all implementation steps are stable.
- **Never commit or push without explicit user approval.**

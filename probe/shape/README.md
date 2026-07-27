# Shape probe — "1-D = scalar" boundary audit

A **diagnostic** (not a test suite) that measures, against the *current, unmodified*
CTDirect, what shapes the transcription actually hands to the user OCP functions stored
in a `CTModels.Model`, and what the cost of each candidate fix would be.

It backs the roadmap item
[#610 §7 / discussion #609 §7](https://github.com/control-toolbox/CTDirect.jl/discussions/609)
and the Handbook rule
[`philosophy/dimension-and-shape.md`](https://github.com/control-toolbox/Handbook/blob/main/philosophy/dimension-and-shape.md):

> A one-dimensional quantity is a `scalar`, never a length-1 vector.

Nothing here asserts pass/fail. Every experiment is wrapped in `try/catch`, so a single
failure never stops the run; each outcome is collected into capability matrices printed
at the bottom of the log.

## Running it

```console
julia --project=probe/shape probe/shape/probe_shape.jl
```

The header activates this folder's environment and `Pkg.develop`s the checked-out
CTDirect, so the probe always measures the working-tree source.

## What it probes

- **A. Boundary shapes today** — a *type-recording* OCP (dynamics, Lagrange, Mayer, path
  and boundary constraints each record the types they were handed, and compute with `[1]`
  so they stay correct for either shape) is driven through `__objective` /
  `__constraints!` for every collocation scheme, once with `n=m=q=1` and once with
  `n=m=q=2`. The cell is `arg1/arg2/arg3` — `x/u/v` for dynamics, Lagrange and path;
  `x0/xf/v` for Mayer and boundary.
- **B. Scalar-written user functions** — an OCP whose dynamics is written in pure scalar
  style (`r[1] = -x + u`, no `[1]` on the *inputs*). This is exactly the code the
  convention exists to make legal.
- **C. Option A audit** — coercing *inside the getters* would hand a scalar to every
  internal consumer too, not only to the user functions. This block measures whether each
  such consumer survives a scalar: the free-time accessors, the solution-building writes,
  the midpoint/trapeze state algebra, and the IRK stage state.
- **D. Solution side** — solve a 1-D problem and inspect what `state(sol)(t)`,
  `control(sol)(t)`, `costate(sol)(t)` return.
- **E. `only` as the coercion** — idempotence on a scalar, behaviour on a `SubArray`
  view, transparency under `ForwardDiff`, behaviour on an *empty* collection, and the
  allocation count.

## Measured (2026-07-27, branch `refactor/docp-dispatch`)

| Block | Result |
| --- | --- |
| A | **Uniformly non-compliant**: all 7 schemes × 5 user functions pass `Vec(1)` for every 1-D quantity. n-D is correct (`Vec(n)`). |
| B | **Fails on every scheme** (`MethodError`): scalar-written user functions are currently illegal. |
| C | Free-time accessors take **both** shapes (`CTModels` has a `ctNumber` method). Broken by a scalar: `build_OCP_solution`'s `v[:] = …` (`ArgumentError`) and trapeze's `x_next += h*(w_i+w_{i+1})` (`MethodError`). Unreachable by any getter: the IRK stage state, which is a **work buffer**. |
| D | **Already compliant** — `CTModels` wraps 1-D trajectories with `only` (`Solutions/interpolation_helpers.jl`), so `state`/`control`/`costate` already return `Number`. |
| E | `only` is idempotent on a scalar, 0-allocation on a length-1 view, AD-transparent, and **throws on an empty collection** — so the dimension-driven map must be `dim == 1 ? only : identity` (never `dim <= 1`). |

The full report built on these numbers lives in `.reports/` (local, git-ignored).

<!-- markdownlint-disable MD024 -->
# Changelog

All notable changes to CTDirect will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.1.3-beta] - 2026-08-28

Aligns every error site in `src/` on the CTBase exception types and adds the missing
`Collocation` / `DirectShooting` docstrings. Dependency-wise it is identical to
1.1.1-beta.

### 🔧 Changed

- **Every error site in `src/` now raises a typed `CTBase.Exceptions` subtype**
  ([#627](https://github.com/control-toolbox/CTDirect.jl/issues/627)), the same alignment
  CTParser did in 0.9.0. Six untyped `error(...)` / `ArgumentError` throws in
  `src/DOCP_data.jl` and `src/ode/common.jl` become `IncorrectArgument` (non-strictly-increasing
  `time_grid`, unknown discretization scheme, unknown getter value) and `NotImplemented`
  (the `DOCP_Jacobian_pattern` / `DOCP_Hessian_pattern` interface stubs), each with
  `got` / `expected` / `suggestion` / `context` fields so callers can `catch` precisely and
  the user gets the structured `Reason / Context / Hint` block. The unknown-scheme
  `expected` is derived from the real dispatch branches rather than a hand-written list
  (the string that drifted in [#624](https://github.com/control-toolbox/CTDirect.jl/issues/624)),
  and the `suggestion` points at the smaller `:exa` set, which CTParser validates
  separately. See the non-breaking note in [BREAKING.md](BREAKING.md).

### 🧪 Testing

- The `:variable`-scheme regression test (`test/ci/test_discretization.jl`) is tightened
  to `@test_throws Exceptions.IncorrectArgument`; new `@test_throws` for a non-increasing
  `time_grid`; `Exceptions` imported in `test/test_common.jl`.
- The `:truck_trailer` `:manual`-backend ADNLP test is disabled, same
  [ADNLPModels.jl#383](https://github.com/JuliaSmoothOptimizers/ADNLPModels.jl/issues/383)
  root cause as `:moonlander` / `:quadrotor` in 1.1.2-beta. Re-enable tracked by
  [#626](https://github.com/control-toolbox/CTDirect.jl/issues/626).

### 📚 Documentation

- **`Collocation` and `DirectShooting` now carry docstrings**
  ([#623](https://github.com/control-toolbox/CTDirect.jl/issues/623)). Both discretizer
  structs were undocumented, so OptimalControl.jl's API reference produced a "no docs
  found" warning when transcluding `CTDirect.Collocation`. Each now states what the
  transcription is, when to prefer it over the other, its supported modeler backends, and
  how to select it from the explicit-mode API.

### ✅ Compatibility

- **No breaking changes.** The error *types* change from `ArgumentError` / `ErrorException`
  to `CTException` subtypes, but the throws fire only on invalid input, so no working code
  can depend on them. See the non-breaking note in [BREAKING.md](BREAKING.md).

## [1.1.2-beta] - 2026-08-27

Drops a discretization scheme that was advertised but never functional, and adds the
`CHANGELOG.md` / `BREAKING.md` pair the Handbook requires. Dependency-wise it is
identical to 1.1.1-beta.

### 🐛 Bug Fixes

- **The `:variable` discretization scheme is removed** ([#624](https://github.com/control-toolbox/CTDirect.jl/issues/624)).
  It was listed in the `:scheme` option description (so `describe(:collocation)` offered
  it) and had a live dispatch branch, but its implementation `src/ode/variable.jl` is not
  compiled — `#include("ode/variable.jl")` in `src/CTDirect.jl` is commented out — so
  `solve(...; scheme=:variable)` failed with a raw
  `UndefVarError(:VariableStepODE, …, CTDirect)` from CTDirect's internals. The
  advertisement and the dead branch are gone; `scheme=:variable` now raises the normal
  `"Unknown discretization method"` error like any other unknown symbol. The WIP file and
  its commented `include` stay in place. See the non-breaking note in
  [BREAKING.md](BREAKING.md).
- **The unknown-scheme error message matches the real dispatch** (`src/DOCP_data.jl`). It
  listed `:gauss_legendre_2_stagewise` / `:gauss_legendre_3_stagewise`, which no branch
  accepts, and omitted `:gauss_legendre_2_constant_control` /
  `:gauss_legendre_3_constant_control`, which are accepted.

### 🧪 Testing

- Regression test for the removed `:variable` scheme (`test/ci/test_discretization.jl`).
- The `:moonlander` and `:quadrotor` `:manual`-backend ADNLP tests are disabled pending
  [JuliaSmoothOptimizers/ADNLPModels.jl#383](https://github.com/JuliaSmoothOptimizers/ADNLPModels.jl/issues/383)
  — `SparseReverseADHessian` records its tape on uninitialised memory, and moonlander's
  unbounded `theta` inside `cos`/`sin` turns a garbage `-Inf` seed into a `DomainError`
  at model-build time. Deterministic on the Julia 1.12 runners, flaky on 1.10. Re-enable
  tracked by [#626](https://github.com/control-toolbox/CTDirect.jl/issues/626).

### 📚 Documentation

- Added `CHANGELOG.md` and `BREAKING.md` (Handbook two-file convention), with a baseline
  entry at 1.0.12.
- `docs/src/index.md`: the default collocation scheme is `:midpoint`, not `:trapeze`.

### ✅ Compatibility

- **No breaking changes.** The `:variable` scheme was never functional, so its removal
  affects no working code. See the non-breaking note in [BREAKING.md](BREAKING.md).

## [1.1.1-beta] - 2026-08-26

Adds the "1-D is a scalar" call-boundary contract and realigns every compat bound with
the stable ecosystem (CTBase 0.29, CTModels 0.18, CTSolvers 0.5) and CTParser 0.9.0-beta.

### 💥 Breaking Changes

#### 1-D state / control / variable reach user functions as scalars

A per-DOCP `DOCPshape` is precomputed and applied immediately before every call into a
user OCP function (dynamics, Lagrange cost, Mayer cost, path and boundary constraints),
across all seven schemes. A state, control or optimisation variable declared with
dimension 1 is now passed as a `Number`, driven by the declared dimension rather than the
runtime type — matching the convention already used by CTModels and CTFlows. The in-place
derivative buffer handed to the dynamics stays a length-1 vector.

**Migration**: problems that index 1-D quantities with `[1]` or call `only(...)` are
unaffected (`x[1] == x` and `length(x) == 1` for a scalar). Code that branches on
`x isa AbstractVector`, or that calls vector-only methods on a 1-D quantity, must accept a
`Number` too — `x isa Union{Number,AbstractVector}`, or normalise with `only(x)`.

```julia
# Before — 1-D state always a length-1 vector
function dynamics!(dx, t, x, u, v)
    dx[1] = x[1] * u[1]
end

# After — 1-D state is a scalar; the derivative buffer is still a vector
function dynamics!(dx, t, x, u, v)
    dx[1] = x * u
end
```

See [BREAKING.md](BREAKING.md).

### 🔧 CI

- The `CI`, `Documentation` and `Breakage` workflows are gated behind PR labels.
- `GPU.yml` (with its weekly cron) is folded into `CI.yml` as a GPU job running the
  dedicated `test/test_gpu.jl` suite; that job then moves off the retired `kkt` runner
  onto `occidata`.
- `test/test_gpu.jl` loads `CUDSS` explicitly — MadNLPGPU 0.10 no longer pulls it in
  transitively, and it is needed for MadNLPGPU's CUDA sparse-KKT extension to arm.

### 📚 Documentation

- `CLAUDE.md` / `AGENTS.md` repository navigation guide.

### 📦 Dependencies

| | from | to |
| --- | --- | --- |
| `CTBase` | `0.28` | `0.29` |
| `CTModels` | `0.15` | `0.18` |
| `CTParser` | `0.8` | `0.9` (0.9.0-beta) |
| `CTSolvers` | `0.4` | `0.5` |
| `ExaModels` | `0.11` | `0.12` |

CTParser 0.9's own breaking changes (the `:exa` emission now requires ExaModels ≥ 0.12;
its `String` errors become `CTException` subtypes) do not reach CTDirect: the ExaModels
builder is a closure from `CTModels.get_build_examodel`, never called directly, and no
test catches a string error from CTParser.

### ✅ Compatibility

- **Breaking**: the 1-D scalar call-boundary contract — see the entry above and
  [BREAKING.md](BREAKING.md). Problems that only index with `[1]` are unaffected. The
  dependency floors are raised; see the non-breaking note in [BREAKING.md](BREAKING.md).

## [1.1.0-beta] - 2026-07-28

Moves the discretization architecture onto the shared CTSolvers / CTBase machinery. This
is the entry point of the 1.1 line and its largest set of breaking changes.

### 💥 Breaking Changes

#### The discretization API follows the CTSolvers dispatch contract

`AbstractDiscretizer` is now owned by CTSolvers; CTDirect's local abstract type and its
`discretize` wrappers are removed, and `Strategies` / `Options` / `Core` / `Exceptions`
are imported from CTBase. The `(discretizer)(ocp)` closure factory is replaced by named
`CTSolvers.discretize` / `CTSolvers.build_model` / `CTSolvers.build_solution` methods
that dispatch on `(DiscretizedModel{<:Discretizer}, modeler)`.

**Migration**:

```julia
# Before
discretizer = CTDirect.Collocation()
docp = CTDirect.discretize(ocp, discretizer)
nlp  = docp.adnlp_model_builder(init; backend=:manual)

# After
discretizer = CTDirect.Collocation()
docp = CTSolvers.discretize(ocp, discretizer)
nlp  = CTSolvers.nlp_model(docp, init, CTSolvers.ADNLP(; backend=:manual))
```

See [BREAKING.md](BREAKING.md).

#### `DOCPCache` is immutable; `build_solution` consumes a `BuiltModel`

`DOCPCache` holds only the DOCP. The ExaModels getter, produced with the ExaModels
object, now rides an `ExaBuildCache` carried by the `CTSolvers.BuiltModel` that
`build_model` returns (`NoCache` for the ADNLP and direct-shooting paths). Nothing
mutates the shared `DiscretizedModel` cache any more, and `build_solution` dispatches on
the `BuiltModel`.

See [BREAKING.md](BREAKING.md).

#### Concrete discretizers must define `Strategies.parameter`

The CTBase 0.28 parameter contract requires
`Strategies.parameter(::Type{<:Collocation}) = nothing` and the same for
`DirectShooting`. Downstream discretizer types must add the method.

See [BREAKING.md](BREAKING.md).

### 📦 Dependencies

| | from | to |
| --- | --- | --- |
| `CTBase` | `0.18` | `0.28` (promoted to `[deps]`) |
| `CTModels` | `0.10` | `0.15` |
| `CTParser` | `[extras]` | `[deps]` |
| `ExaModels` | `0.9` | `0.11` |
| `CUDA` (test) | `5` | `5, 6` |
| `MadNLP` (test) | `0.9` | `0.9, 0.10` |
| `MadNLPGPU` (test) | `0.8` | `0.8, 0.10` |

### ✅ Compatibility

- **Breaking**: see the three entries above and [BREAKING.md](BREAKING.md). CTDirect is a
  low-level package; the affected surface is CTSolvers, OptimalControl.jl and any custom
  discretizer, not `@def` problem definitions.

## [1.0.12] - 2026-05-08 — baseline

This is the reference version. No changelog was maintained before this point; use
`git log` for earlier history. Breaking changes from this version onward are tracked in
[BREAKING.md](BREAKING.md).

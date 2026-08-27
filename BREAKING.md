<!-- markdownlint-disable MD024 -->
# Breaking Changes

Breaking changes in CTDirect releases, and how to migrate. Tracked from the 1.0.12
baseline onward; see [CHANGELOG.md](CHANGELOG.md) for the full record.

## [1.1.2-beta] - 2026-08-27

No breaking changes in this release.

## Non-breaking note (1.1.2-beta)

- **The `:variable` discretization scheme is removed.** It was advertised in the
  `:scheme` option description (and therefore in `describe(:collocation)`) and had a live
  dispatch branch, but its implementation was never compiled — every
  `solve(...; scheme=:variable)` raised `UndefVarError(:VariableStepODE, …, CTDirect)`.
  **No breaking change**: no working code can have depended on it. `scheme=:variable` now
  raises the standard `"Unknown discretization method"` error. The `src/ode/variable.jl`
  WIP file and its commented `#include` are kept for a future variable-step
  implementation.

## [1.1.1-beta] - 2026-08-26

### 1-D state / control / variable are scalars at the user-function boundary

CTDirect coerces every quantity crossing into a user OCP function (dynamics, Lagrange
cost, Mayer cost, path and boundary constraints) to its declared shape, for all seven
discretization schemes. A state, control or optimisation variable of dimension 1 is now
passed as a `Number` rather than a length-1 `Vector`, driven by the declared dimension
and not the runtime type. This matches the "1-D is a scalar" convention already used by
CTModels and CTFlows. The in-place derivative buffer given to the dynamics stays a
length-1 vector.

**Who is affected**: OCP function definitions that branch on `x isa AbstractVector` for a
1-D quantity, call vector-only methods on it, or otherwise assume it is always a
container. Definitions that only index with `[1]`, iterate, or take `length` are
unaffected — `x[1] == x`, `length(x) == 1` and `only(x) == x` all hold for a scalar.

**Migration**:

```julia
# Before — 1-D state is a length-1 vector
function dynamics!(dx, t, x, u, v)
    dx[1] = x[1] * u[1]
end

# After — 1-D state and control are scalars; dx stays a vector
function dynamics!(dx, t, x, u, v)
    dx[1] = x * u
end
```

```julia
# Before
val = x isa AbstractVector ? x[1] : x

# After
val = only(x)                       # normalises scalar and length-1 vector alike
```

## Non-breaking note (1.1.1-beta)

- **Dependency floors raised**: `CTBase 0.29`, `CTModels 0.18`, `CTSolvers 0.5`,
  `CTParser 0.9` (resolves `0.9.0-beta` from the ct-registry), `ExaModels 0.12`. **No
  breaking change**: no CTDirect source change was required, and CTParser 0.9's own
  breaking changes do not reach CTDirect — it never calls the ExaModels builder API
  directly (the closure comes from `CTModels.get_build_examodel`) and no test catches a
  string error from CTParser. Callers pinning older bounds must move them up.

## [1.1.0-beta] - 2026-07-28

### The discretization API follows the CTSolvers dispatch contract

`AbstractDiscretizer` moved to CTSolvers. CTDirect's local abstract type and its
`discretize` wrappers were removed, and the `(discretizer)(ocp)` closure-factory pattern
was replaced by named methods — `CTSolvers.discretize`, `CTSolvers.build_model`,
`CTSolvers.build_solution` — dispatching on `(DiscretizedModel{<:Discretizer}, modeler)`.
`Strategies`, `Options`, `Core` and `Exceptions` are imported from CTBase.

**Who is affected**: code that calls `CTDirect.discretize` or builds an NLP through the
`docp.*_model_builder` fields — i.e. CTSolvers, OptimalControl.jl and any custom
discretizer. `@def` problem definitions are unaffected.

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

### `DOCPCache` is immutable; `build_solution` consumes a `BuiltModel`

`DOCPCache` now holds only the DOCP. The ExaModels getter rides an `ExaBuildCache`
carried by the `CTSolvers.BuiltModel` returned by `build_model` (`NoCache` for the ADNLP
and direct-shooting paths); the shared `DiscretizedModel` cache is no longer mutated
between `build_model` and `build_solution`.

**Who is affected**: code constructing `DOCPCache` directly or reading the getter off the
discretized-model cache.

**Migration**:

```julia
# Before — getter mutated into the shared cache
built = build_model(dm, modeler)
getter = dm.cache.exa_getter

# After — getter travels with the BuiltModel
built  = CTSolvers.build_model(dm, modeler)
getter = built.cache.exa_getter        # ExaBuildCache; NoCache for ADNLP
```

### Concrete discretizers must define `Strategies.parameter`

The CTBase 0.28 parameter contract requires every discretizer type to answer
`Strategies.parameter`. CTDirect adds it for `Collocation` and `DirectShooting`.

**Who is affected**: packages defining their own `<: CTSolvers.AbstractDiscretizer`.

**Migration**:

```julia
# Add for each custom discretizer type
CTBase.Strategies.parameter(::Type{<:MyDiscretizer}) = nothing
```

## [1.0.12] - 2026-05-08 — baseline

Reference version. Breaking changes are tracked from here onward; use `git log` for
earlier history.

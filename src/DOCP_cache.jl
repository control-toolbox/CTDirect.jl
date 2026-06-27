# DOCP cache
#
# Backend cache attached to a `CTSolvers.DiscretizedModel`.

"""
$(TYPEDEF)

Mutable cache attached to a `CTSolvers.DiscretizedModel`.

Created by `CTSolvers.discretize`, then partially mutated by
`CTSolvers.build_model` for the Exa backend (which sets `exa_getter`). Making
the `build_model` → `build_solution` coupling explicit through this cache
replaces the previously captured closure variable.

# Fields
- `docp::D`: The internal discretized OCP structure (bounds, dimensions, ...).
- `exa_getter`: Getter produced by the ExaModels constructor; `nothing` before an
  Exa model is built. Typed `Any` because the cache is created (with `nothing`) by
  `discretize`, then reassigned to a backend-specific getter by `build_model`; a
  concrete field parameter would freeze it to `Nothing` and forbid the reassignment.
"""
mutable struct DOCPCache{D<:DOCP} <: Core.AbstractCache
    docp::D
    exa_getter::Any
end

DOCPCache(docp::DOCP) = DOCPCache(docp, nothing)

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
- `exa_getter::E`: Getter produced by the ExaModels constructor; `nothing`
  before an Exa model is built.
"""
mutable struct DOCPCache{D<:DOCP,E} <: Core.AbstractCache
    docp::D
    exa_getter::E
end

DOCPCache(docp::DOCP) = DOCPCache(docp, nothing)

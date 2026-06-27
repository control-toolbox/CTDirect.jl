# DOCP cache
#
# Backend cache attached to a `CTSolvers.DiscretizedModel`.

"""
$(TYPEDEF)

Immutable discretize-time cache attached to a `CTSolvers.DiscretizedModel`.

Created by `CTSolvers.discretize`, it holds the precomputed [`DOCP`](@ref) reused
across model/solution builds. It is never mutated: any build-time auxiliary (e.g.
the ExaModels getter) is carried by the [`ExaBuildCache`](@ref) inside the
`CTSolvers.BuiltModel` returned by `CTSolvers.build_model`.

# Fields
- `docp::D`: The internal discretized OCP structure (bounds, dimensions, ...).
"""
struct DOCPCache{D<:DOCP} <: Core.AbstractCache
    docp::D
end

"""
$(TYPEDEF)

Immutable build-time cache produced by `CTSolvers.build_model` for the Exa backend.

Carries the getter produced together with the `ExaModel` so that
`CTSolvers.build_solution` can reconstruct the OCP solution. Stored in the
`cache` field of the `CTSolvers.BuiltModel`; never mutated.

# Fields
- `exa_getter::G`: Getter produced by the ExaModels constructor.
"""
struct ExaBuildCache{G} <: Core.AbstractCache
    exa_getter::G
end

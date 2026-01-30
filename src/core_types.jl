# ---------------------------------------------------------------------------
# Discretization schemes: see disc/
# ---------------------------------------------------------------------------
"""
$(TYPEDEF)

Abstract type representing a discretization strategy for an optimal
control problem.  

Concrete subtypes of `Discretization` define specific schemes for
transforming a continuous-time problem into a discrete-time
representation suitable for numerical solution.

# Example

```julia-repl
julia> struct MyDiscretization <: Discretization end
MyDiscretization
```
"""
abstract type Discretization end

# ---------------------------------------------------------------------------
# Abstract discretizer type
# ---------------------------------------------------------------------------
# must be a subtype of CTModels.AbstractOCPTool to get benefit of options handling
abstract type AbstractOptimalControlDiscretizer <: CTModels.AbstractOCPTool end

# ---------------------------------------------------------------------------
# Collocation discretizer
# ---------------------------------------------------------------------------
mutable struct Collocation{D<:CTDirect.Discretization} <: AbstractOptimalControlDiscretizer
    #{D<:CTDirect.Discretization, O<:CTModels.Model} 
    
    # required to be able to use default CTModels.AbstractOCPTool getters
    options_values
    options_sources

    disc::D
    docp
    exa_getter
    #+++ put here contents from CTDirect.DOCP and adjust functions ?
end

# useful for OptimalControl. Should be a field of AbstractOptimalControlDiscretizer with default getter, for consistency.
CTModels.get_symbol(::Type{<:Collocation}) = :collocation

# default options (this replaces default.jl)
__grid_size()::Int = 250
__scheme()::Symbol = :midpoint
__grid()::Union{Int,AbstractVector} = __grid_size()

# options specs: for each option, we define the type, default value, and description.
function CTModels._option_specs(::Type{<:Collocation})
    return (
        grid=CTModels.OptionSpec(;
            type=Union{Int,AbstractVector},
            default=__grid(),
            description="Collocation grid (Int = number of time steps, Vector = explicit time grid).",
        ),
        scheme=CTModels.OptionSpec(;
            type=Symbol,
            default=__scheme(),
            description="Time integration scheme used by the collocation discretizer.",
        ),
    )
end

# constructor: kwargs contains the options values
function Collocation(; kwargs...)

    # parsing options from CTModels
    values, sources = CTModels._build_ocp_tool_options(Collocation; kwargs..., strict_keys=true)
    #scheme = values.scheme

    # disc
    disc, _, _ = CTDirect.Midpoint(1,1,1,1,1,1)
    # docp
    docp = nothing
    # exa getter
    exa_getter = nothing

    return Collocation{typeof(disc)}(values, sources, disc, docp, exa_getter)
end

# ---------------------------------------------------------------------------
# Discretizers registration
# ---------------------------------------------------------------------------
# useful for OptimalControl.
const REGISTERED_DISCRETIZERS = (Collocation,) # TODO: add DirectShooting
registered_discretizer_types() = REGISTERED_DISCRETIZERS
discretizer_symbols() = Tuple(CTModels.get_symbol(T) for T in REGISTERED_DISCRETIZERS)
function _discretizer_type_from_symbol(sym::Symbol)
    for T in REGISTERED_DISCRETIZERS
        if CTModels.get_symbol(T) === sym
            return T
        end
    end
    msg = "Unknown discretizer symbol $(sym). Supported discretizers: $(discretizer_symbols())."
    throw(CTBase.IncorrectArgument(msg))
end
function build_discretizer_from_symbol(sym::Symbol; kwargs...)
    T = _discretizer_type_from_symbol(sym)
    return T(; kwargs...)
end
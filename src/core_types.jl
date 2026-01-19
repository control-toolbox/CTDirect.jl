abstract type AbstractIntegratorScheme end
struct Midpoint <: AbstractIntegratorScheme end
struct Trapezoidal <: AbstractIntegratorScheme end
const Trapeze = Trapezoidal

abstract type AbstractOptimalControlDiscretizer <: CTModels.AbstractOCPTool end

struct Collocation{T<:AbstractIntegratorScheme} <: AbstractOptimalControlDiscretizer
    options_values
    options_sources
end

__grid_size()::Int = 250
__scheme()::AbstractIntegratorScheme = Midpoint()

__grid()::Union{Int,AbstractVector} = __grid_size()

function _option_specs(::Type{<:Collocation})
    return (
        grid=OptionSpec(;
            type=Union{Int,AbstractVector},
            default=__grid(),
            description="Collocation grid (Int = number of time steps, Vector = explicit time grid).",
        ),
        lagrange_to_mayer=OptionSpec(;
            type=Bool,
            default=false,
            description="Whether to transform the Lagrange integral cost into an equivalent Mayer terminal cost.",
        ),
        scheme=OptionSpec(;
            type=AbstractIntegratorScheme,
            default=__scheme(),
            description="Time integration scheme used by the collocation discretizer.",
        ),
    )
end

function Collocation(; kwargs...)
    values, sources = _build_ocp_tool_options(Collocation; kwargs..., strict_keys=true)
    scheme = values.scheme
    return Collocation{typeof(scheme)}(values, sources)
end

get_symbol(::Type{<:Collocation}) = :collocation

const REGISTERED_DISCRETIZERS = (Collocation,)

registered_discretizer_types() = REGISTERED_DISCRETIZERS

discretizer_symbols() = Tuple(get_symbol(T) for T in REGISTERED_DISCRETIZERS)

function _discretizer_type_from_symbol(sym::Symbol)
    for T in REGISTERED_DISCRETIZERS
        if get_symbol(T) === sym
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
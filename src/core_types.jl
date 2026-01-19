# ---------------------------------------------------------------------------
# Discretization schemes
# ---------------------------------------------------------------------------
# to merge with structs in `disc/`
abstract type AbstractIntegratorScheme end
struct MidpointScheme <: AbstractIntegratorScheme end
struct TrapezoidalScheme <: AbstractIntegratorScheme end
const TrapezeScheme = TrapezoidalScheme

# ---------------------------------------------------------------------------
# Abstract discretizer type
# ---------------------------------------------------------------------------
# must be a subtype of CTModels.AbstractOCPTool to get benefit of options handling
abstract type AbstractOptimalControlDiscretizer <: CTModels.AbstractOCPTool end

# ---------------------------------------------------------------------------
# Collocation discretizer
# ---------------------------------------------------------------------------
# must have the fields `options_values` and `options_sources` to be compatible with CTModels.AbstractOCPTool getters
struct Collocation{T<:AbstractIntegratorScheme} <: AbstractOptimalControlDiscretizer
    options_values
    options_sources
end

# useful for OptimalControl. Should be a field of AbstractOptimalControlDiscretizer with default getter, for consistency.
CTModels.get_symbol(::Type{<:Collocation}) = :collocation

# default options (related to default.jl)
__grid_size()::Int = 250
__scheme()::AbstractIntegratorScheme = MidpointScheme()
__grid()::Union{Int,AbstractVector} = __grid_size()

# options specs: for each option, we define the type, default value, and description.
function CTModels._option_specs(::Type{<:Collocation})
    return (
        grid=CTModels.OptionSpec(;
            type=Union{Int,AbstractVector},
            default=__grid(),
            description="Collocation grid (Int = number of time steps, Vector = explicit time grid).",
        ),
        lagrange_to_mayer=CTModels.OptionSpec(;
            type=Bool,
            default=false,
            description="Whether to transform the Lagrange integral cost into an equivalent Mayer terminal cost.",
        ),
        scheme=CTModels.OptionSpec(;
            type=AbstractIntegratorScheme, # maybe we should use a Symbol instead?
            default=__scheme(),
            description="Time integration scheme used by the collocation discretizer.",
        ),
    )
end

# constructor: kwargs contains the options values
function Collocation(; kwargs...)
    # parsing options from CTModels
    values, sources = CTModels._build_ocp_tool_options(Collocation; kwargs..., strict_keys=true)
    scheme = values.scheme
    return Collocation{typeof(scheme)}(values, sources)
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
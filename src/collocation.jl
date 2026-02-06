# ---------------------------------------------------------------------------
# Implementation of Collocation discretizer
# ---------------------------------------------------------------------------

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


function discretize(
    ocp::AbstractOptimalControlProblem, 
    discretizer::AbstractOptimalControlDiscretizer
)
    return discretizer(ocp)
end

__discretizer()::AbstractOptimalControlDiscretizer = Collocation()

function discretize(
    ocp::AbstractOptimalControlProblem;
    discretizer::AbstractOptimalControlDiscretizer=__discretizer(),
)
    return discretize(ocp, discretizer)
end

# ---------------------------------------------------------------------------
# Collocation discretizer
# ---------------------------------------------------------------------------
struct Collocation <: AbstractOptimalControlDiscretizer
    # required to be able to use default CTModels.AbstractOCPTool getters
    options_values
    options_sources
end

# useful for OptimalControl. Should be a field of AbstractOptimalControlDiscretizer with default getter, for consistency.
CTModels.get_symbol(::Type{<:Collocation}) = :collocation

# default options
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

    return Collocation(values, sources)
end

# ---------------------------------------------------------------------------
# Discretizers registration
# ---------------------------------------------------------------------------
# useful for OptimalControl.
const REGISTERED_DISCRETIZERS = (Collocation,)
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


# default options for modelers backend
# +++ recheck kwargs passing / default with Olivier
__adnlp_backend() = :optimized
__exa_backend() = nothing

# ==========================================================================================
# Scheme symbol mapping
# ==========================================================================================
function get_scheme(discretizer::Collocation)
    return CTModels.get_option_value(discretizer, :scheme)
end

# ==========================================================================================
# Grid options mapping
# unified grid option: Int => grid_size, Vector => explicit time_grid
# ==========================================================================================
function grid_options(discretizer::Collocation)

    grid = CTModels.get_option_value(discretizer, :grid)
    if grid isa Int
        grid_size = grid
        time_grid = nothing
    else
        grid_size = length(grid)
        time_grid = grid
    end

    return grid_size, time_grid
end

# ==========================================================================================
# Build core DOCP structure with discretization information (ADNLP)
# ==========================================================================================
function get_docp(discretizer::Collocation, ocp::AbstractOptimalControlProblem)
    
    # recover discretization scheme and options
    scheme = get_scheme(discretizer)
    grid_size, time_grid = grid_options(discretizer)

    # initialize DOCP
    docp = DOCP(ocp; grid_size=grid_size, time_grid=time_grid, scheme=scheme)

    # set bounds in DOCP
    variables_bounds!(docp)
    constraints_bounds!(docp)

    return docp
end

# ==========================================================================================
# Build initial guess for discretized problem
# ==========================================================================================
function get_docp_initial_guess(modeler::Symbol, docp,
        initial_guess::Union{CTModels.AbstractOptimalControlInitialGuess,Nothing},
        )

        ocp = ocp_model(docp)

        # build functional initial guess
        functional_init = CTModels.build_initial_guess(ocp, initial_guess)

        # build discretized initial guess
        x0 = DOCP_initial_guess(docp, functional_init)
   
        if modeler == :adnlp
            return x0
        elseif modeler == :exa
            # reshape initial guess for ExaModel variables layout
            # - do not broadcast, apparently fails on GPU arrays
            # - unused final control in examodel / euler, hence the different x0 sizes
            n = CTModels.state_dimension(ocp)
            m = CTModels.control_dimension(ocp)
            q = CTModels.variable_dimension(ocp)
            N = docp.time.steps
            # N + 1 states, N controls
            state = hcat([x0[(1 + i * (n + m)):(1 + i * (n + m) + n - 1)] 
            for i in 0:N]...)
            control = hcat([x0[(n + 1 + i * (n + m)):(n + 1 + i * (n + m) + m - 1)] 
            for i in 0:(N - 1)]...,)
            # +++? todo: pass indeed to grid_size only for euler(_b), trapeze and midpoint
            control = [control control[:, end]] 
            variable = x0[(end - q + 1):end]
            
            return (variable, state, control)
        else
            error("unknown modeler in get_docp_initial_guess: ", modeler)
        end
    end


# ==========================================================================================
# Build discretizer API (return sets of model/solution builders)
# ==========================================================================================
function (discretizer::Collocation)(ocp::AbstractOptimalControlProblem)

    # common parts for builders
    docp = get_docp(discretizer, ocp)
    exa_getter = nothing # will be set in build_exa_model

    # ==========================================================================================
    # The needed builders for the construction of the final DiscretizedOptimalControlProblem
    # ==========================================================================================
    # +++ recheck kwargs passing / default with Olivier
    function build_adnlp_model(
        initial_guess::CTModels.AbstractOptimalControlInitialGuess;
        adnlp_backend=__adnlp_backend(),
        show_time=false,
        kwargs...
    )::ADNLPModels.ADNLPModel

        # functions for objective and constraints
        f = x -> CTDirect.DOCP_objective(x, docp)
        c! = (c, x) -> CTDirect.DOCP_constraints!(c, x, docp)

        # build initial guess
        init = get_docp_initial_guess(:adnlp, docp, initial_guess)

        # unused backends (option excluded_backend = [:jprod_backend, :jtprod_backend, :hprod_backend, :ghjvprod_backend] does not seem to work)
        unused_backends = (
        hprod_backend=ADNLPModels.EmptyADbackend,
        jtprod_backend=ADNLPModels.EmptyADbackend,
        jprod_backend=ADNLPModels.EmptyADbackend,
        ghjvprod_backend=ADNLPModels.EmptyADbackend,
        )

        # set adnlp backends
        if adnlp_backend == :manual

            # build sparsity patterns for Jacobian and Hessian
            J_backend = ADNLPModels.SparseADJacobian(
                docp.dim_NLP_variables, f,
                docp.dim_NLP_constraints, c!,
                CTDirect.DOCP_Jacobian_pattern(docp),
            )
            H_backend = ADNLPModels.SparseReverseADHessian(
                docp.dim_NLP_variables, f,
                docp.dim_NLP_constraints, c!,
                CTDirect.DOCP_Hessian_pattern(docp),
            )
            backend_options = (
                gradient_backend=ADNLPModels.ReverseDiffADGradient,
                jacobian_backend=J_backend,
                hessian_backend=H_backend,
            )
        else
            # use backend preset
            backend_options = (backend=adnlp_backend,)
        end

        # build NLP
        nlp = ADNLPModel!(
            f,
            init,
            docp.bounds.var_l,
            docp.bounds.var_u,
            c!,
            docp.bounds.con_l,
            docp.bounds.con_u;
            minimize=(!docp.flags.max),
            backend_options...,
            unused_backends...,
            show_time=show_time,
        )

        return nlp
    end

    # Solution builder for ADNLPModels
    function build_adnlp_solution(nlp_solution::SolverCore.AbstractExecutionStats)
        
        # retrieve data from NLP solver
        minimize = !docp.flags.max
        objective, iterations, constraints_violation, message, status, successful = CTModels.extract_solver_infos(nlp_solution, minimize)

        # retrieve time grid
        T = get_time_grid(nlp_solution.solution, docp)

        # build OCP solution from NLP solution
        sol = CTDirect.build_OCP_solution(docp, nlp_solution, T, 
        objective, iterations, constraints_violation, message, status, successful)
        
        return sol
    end

    # NLP builder for ExaModels
    # +++ recheck kwargs passing / default with Olivier
    function build_exa_model(
        ::Type{BaseType}, 
        initial_guess::CTModels.AbstractOptimalControlInitialGuess; 
        exa_backend=CTDirect.__exa_backend(),
        kwargs...
    )::ExaModels.ExaModel where {BaseType<:AbstractFloat}

        # recover discretization scheme and options
        # since exa part does not reuse the docp struct
        scheme = get_scheme(discretizer)
        grid_size, time_grid = grid_options(discretizer)

        # build initial guess
        init = get_docp_initial_guess(:exa, docp, initial_guess)

        # build Exa model and getters
        # +++ later try to call Exa constructor here if possible, reusing existing functions...
        build_exa = CTModels.get_build_examodel(ocp)
        nlp, exa_getter = build_exa(;
            grid_size=grid_size,
            backend=exa_backend,
            scheme=scheme,
            init=init,
        )

        return nlp
    end

    # Solution builder for ExaModels
    function build_exa_solution(nlp_solution::SolverCore.AbstractExecutionStats)

        # NB exa_getter is set during build_exa_model call !
        if isnothing(exa_getter)
            error("build_exa_solution: exa_getter is nothing")
        end

        # retrieve data from NLP solver
        minimize = !docp.flags.max
        objective, iterations, constraints_violation, message, status, successful = CTModels.extract_solver_infos(nlp_solution, minimize)
  
        # retrieve time grid
        T = get_time_grid_exa(nlp_solution, docp, exa_getter)

        # build OCP solution from NLP solution
        sol = CTDirect.build_OCP_solution(docp, nlp_solution, T,
        objective, iterations, constraints_violation, message, status, successful; 
        exa_getter=exa_getter)
        
        return sol
    end

    #NB. it would be better to return builders as model/solution pairs since they are linked
    return CTModels.DiscretizedOptimalControlProblem(
        ocp,
        CTModels.ADNLPModelBuilder(build_adnlp_model),
        CTModels.ExaModelBuilder(build_exa_model),
        CTModels.ADNLPSolutionBuilder(build_adnlp_solution),
        CTModels.ExaSolutionBuilder(build_exa_solution),
    )
end

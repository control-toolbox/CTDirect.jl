# ---------------------------------------------------------------------------
# Implementation of Direct shooting discretizer
# ---------------------------------------------------------------------------
struct DirectShooting <: AbstractDiscretizer
    options::Strategies.StrategyOptions
end

# useful for OptimalControl
Strategies.id(::Type{<:DirectShooting}) = :direct_shooting

# default options
__direct_shooting_grid_size()::Int = 100
__direct_shooting_scheme()::Symbol = :midpoint # later use variable step ode solver
__direct_shooting_scheme_steps() = 10

function Strategies.metadata(::Type{<:Collocation})
    return Strategies.StrategyMetadata(
        Options.OptionDefinition(
        name = :grid_size,
        type = Int,
        default = __direct_shooting_grid_size(),
        description = "Number of time steps for the direct shooting grid",
        ),
        Options.OptionDefinition(
        name = :scheme,
        type = Symbol,
        default = __direct_shooting_scheme(),
        description = "Time integration scheme (e.g., :midpoint, :trapeze)",
        ),
    )
end

# constructor: kwargs contains the options values
function DirectShooting(; mode::Symbol = :strict, kwargs...)
    opts = Strategies.build_strategy_options(DirectShooting; mode = mode, kwargs...)
    return DirectShooting(opts)
end

Strategies.options(c::DirectShooting) = c.options


# ==========================================================================================
# Build discretizer API (return sets of model/solution builders)
# ==========================================================================================
function (discretizer::DirectShooting)(ocp::AbstractModel)

    # common parts for builders
    docp = get_docp(discretizer, ocp) ?
    exa_getter = nothing # will be set in build_exa_model

    # ==========================================================================================
    # The needed builders for the construction of the final DiscretizedModel
    # ==========================================================================================
    
    # NLP builder for ADNLPModels
    function build_adnlp_model(
        initial_guess::CTModels.AbstractInitialGuess;
        backend,
        kwargs...
    )::ADNLPModels.ADNLPModel

        # functions for objective and constraints
        f = x -> CTDirect.DirectShooting_objective(x, docp)
        c! = (c, x) -> CTDirect.DirectShooting_constraints!(c, x, docp)

        # build initial guess
        init = get_docp_initial_guess(:adnlp, docp, initial_guess)

        # unused backends (option excluded_backend = [:jprod_backend, :jtprod_backend, :hprod_backend, :ghjvprod_backend] does not seem to work)
        unused_backends = (
        hprod_backend=ADNLPModels.EmptyADbackend,
        jtprod_backend=ADNLPModels.EmptyADbackend,
        jprod_backend=ADNLPModels.EmptyADbackend,
        ghjvprod_backend=ADNLPModels.EmptyADbackend,
        )


        # use backend preset
        backend_options = (backend=backend,)

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
            kwargs...,
        )

        return nlp
    end

    # Solution builder for ADNLPModels
    function build_adnlp_solution(nlp_solution::SolverCore.AbstractExecutionStats)
        
        # retrieve data from NLP solver
        minimize = !docp.flags.max
        objective, iterations, constraints_violation, message, status, successful = CTSolvers.extract_solver_infos(nlp_solution, minimize)

        # retrieve time grid
        T = get_time_grid(nlp_solution.solution, docp)

        # build OCP solution from NLP solution
        sol = CTDirect.build_OCP_solution(docp, nlp_solution, T, 
        objective, iterations, constraints_violation, message, status, successful)
        
        return sol
    end

    # NLP builder for ExaModels
    function build_exa_model(
        ::Type{BaseType}, 
        initial_guess::CTModels.AbstractInitialGuess; 
        backend
    )::ExaModels.ExaModel where {BaseType<:AbstractFloat}
    end

    # Solution builder for ExaModels
    function build_exa_solution(nlp_solution::SolverCore.AbstractExecutionStats)
    end

    #NB. it would be better to return builders as model/solution pairs since they are linked
    return CTSolvers.DiscretizedModel(
        ocp,
        CTSolvers.ADNLPModelBuilder(build_adnlp_model),
        CTSolvers.ExaModelBuilder(build_exa_model),
        CTSolvers.ADNLPSolutionBuilder(build_adnlp_solution),
        CTSolvers.ExaSolutionBuilder(build_exa_solution),
    )
end

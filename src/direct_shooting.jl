# ---------------------------------------------------------------------------
# Direct shooting discretizer
#
# Implements the CTSolvers contract (discretize / build_model / build_solution)
# for the DirectShooting transcription, for the ADNLP backend only. The Exa
# backend is not implemented: the CTSolvers NotImplemented stub applies.
# ---------------------------------------------------------------------------
import SparseConnectivityTracer.TracerLocalSparsityDetector

struct DirectShooting <: CTSolvers.AbstractDiscretizer
    options::Strategies.StrategyOptions
end

# useful for OptimalControl
Strategies.id(::Type{<:DirectShooting}) = :direct_shooting

# default options
__direct_shooting_grid_size()::Int = 250
__direct_shooting_control_steps()::Int = 1 # ie number of controls per time step
__direct_shooting_scheme()::Symbol = :midpoint # later use variable step ode solver

function Strategies.metadata(::Type{<:DirectShooting})
    return Strategies.StrategyMetadata(
        Options.OptionDefinition(
        name = :grid_size,
        type = Int,
        default = __direct_shooting_grid_size(),
        description = "Number of time steps for the direct shooting grid",
        ),
        Options.OptionDefinition(
        name = :control_steps,
        type = Int,
        default = __direct_shooting_control_steps(),
        description = "Number of controls per time step for the direct shooting",
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
# CTSolvers contract: discretize
# ==========================================================================================
"""
$(TYPEDSIGNATURES)

Discretize an OCP with the DirectShooting strategy into a `CTSolvers.DiscretizedModel`
holding a [`DOCPCache`](@ref) with the precomputed DOCP.
"""
function CTSolvers.discretize(ocp::AbstractModel, discretizer::DirectShooting)
    docp = get_docp(discretizer, ocp)
    return CTSolvers.DiscretizedModel(ocp, discretizer, DOCPCache(docp))
end

# ==========================================================================================
# CTSolvers contract: ADNLP backend
# ==========================================================================================
"""
$(TYPEDSIGNATURES)

Build a `CTSolvers.BuiltModel` wrapping an `ADNLPModel` for a DirectShooting-discretized
problem. The ADNLP backend needs no build-time auxiliary, hence `CTSolvers.NoCache`.
"""
function CTSolvers.build_model(
    dm::CTSolvers.DiscretizedModel{<:Any,<:DirectShooting},
    initial_guess::CTModels.AbstractInitialGuess,
    modeler::CTSolvers.Modelers.ADNLP,
)
    docp = dm.cache.docp
    ocp = dm.ocp

    # modeler options; backend is consumed here, the rest is forwarded to ADNLPModel!
    options = Strategies.options_dict(modeler)
    backend = pop!(options, :backend)

    # functions for objective and constraints
    f = x -> __objective(x, docp)
    c! = (c, x) -> __constraints!(c, x, docp)

    # build initial guess
    functional_init = CTModels.build_initial_guess(ocp, initial_guess)
    x0 = __initial_guess(docp, functional_init)

    # unused backends (option excluded_backend = [:jprod_backend, :jtprod_backend, :hprod_backend, :ghjvprod_backend] does not seem to work)
    unused_backends = (
    hprod_backend=ADNLPModels.EmptyADbackend,
    jtprod_backend=ADNLPModels.EmptyADbackend,
    jprod_backend=ADNLPModels.EmptyADbackend,
    ghjvprod_backend=ADNLPModels.EmptyADbackend,
    )

    # use backend preset
    # NB. problems with variable step, even with :generic (also, dense...)
    backend_options = (backend=backend,)

    #= error
    backend_options = (
    jacobian_backend = ADNLPModels.SparseADJacobian(
        docp.dim_NLP_variables, f,
        docp.dim_NLP_constraints, c!,
        detector=TracerLocalSparsityDetector()),
    hessian_backend = ADNLPModels.SparseADHessian(
        docp.dim_NLP_variables, f,
        docp.dim_NLP_constraints, c!,
        detector=TracerLocalSparsityDetector()),
    )=#

    # build NLP
    nlp = ADNLPModels.ADNLPModel!(
        f,
        x0,
        docp.bounds.var_l,
        docp.bounds.var_u,
        c!,
        docp.bounds.con_l,
        docp.bounds.con_u;
        minimize=(!docp.flags.max),
        backend_options...,
        unused_backends...,
        options...,
    )

    return CTSolvers.BuiltModel(dm, nlp, CTSolvers.NoCache())
end

"""
$(TYPEDSIGNATURES)

Build an OCP solution from an ADNLP solver result for a DirectShooting-discretized problem.
"""
function CTSolvers.build_solution(
    built::CTSolvers.BuiltModel{<:CTSolvers.DiscretizedModel{<:Any,<:DirectShooting}},
    nlp_solution::SolverCore.AbstractExecutionStats,
    ::CTSolvers.Modelers.ADNLP,
)
    docp = built.problem.cache.docp

    # retrieve data from NLP solver
    objective, iterations, constraints_violation, message, status, successful = CTSolvers.extract_solver_infos(nlp_solution)

    # retrieve time grid +++ put in build_OCP_solution
    T = get_time_grid(nlp_solution.solution, docp)

    # build OCP solution from NLP solution
    sol = build_OCP_solution(docp, nlp_solution, T,
    objective, iterations, constraints_violation, message, status, successful)

    return sol
end

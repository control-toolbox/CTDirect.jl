# ---------------------------------------------------------------------------
# Collocation discretizer
#
# Implements the CTSolvers contract (discretize / build_model / build_solution)
# for the Collocation transcription, for the ADNLP and Exa backends.
# ---------------------------------------------------------------------------
struct Collocation <: CTSolvers.AbstractDiscretizer
    options::Strategies.StrategyOptions
end

# useful for OptimalControl
Strategies.id(::Type{<:Collocation}) = :collocation

# default options
__collocation_grid_size()::Int = 250
__collocation_scheme()::Symbol = :midpoint
__collocation_time_grid() = nothing

function Strategies.metadata(::Type{<:Collocation})
    return Strategies.StrategyMetadata(
        Options.OptionDefinition(
        name = :grid_size,
        type = Int,
        default = __collocation_grid_size(),
        description = "Number of time steps for the collocation grid",
        ),
        Options.OptionDefinition(
        name = :scheme,
        type = Symbol,
        default = __collocation_scheme(),
        aliases=(:disc_method,),
        description = "Time integration scheme (:trapeze, :midpoint, :euler (or :euler_explicit, :euler_forward), :euler_implicit (or :euler_backward), :gauss_legendre_2, :gauss_legendre_3, :variable)",
        ),
        Options.OptionDefinition(
        name = :time_grid,
        type = Union{Nothing,AbstractVector},
        default = __collocation_time_grid(),
        description = "Explicit time grid (possibly non uniform) for the collocation",
        ),
    )
end

# constructor: kwargs contains the options values
function Collocation(; mode::Symbol = :strict, kwargs...)
    opts = Strategies.build_strategy_options(Collocation; mode = mode, kwargs...)
    return Collocation(opts)
end

Strategies.options(c::Collocation) = c.options

# ==========================================================================================
# CTSolvers contract: discretize
# ==========================================================================================
"""
$(TYPEDSIGNATURES)

Discretize an OCP with the Collocation strategy into a `CTSolvers.DiscretizedModel`
holding a [`DOCPCache`](@ref) with the precomputed DOCP.
"""
function CTSolvers.discretize(ocp::AbstractModel, discretizer::Collocation)
    docp = get_docp(discretizer, ocp)
    return CTSolvers.DiscretizedModel(ocp, discretizer, DOCPCache(docp))
end

# ==========================================================================================
# CTSolvers contract: ADNLP backend
# ==========================================================================================
"""
$(TYPEDSIGNATURES)

Build an `ADNLPModel` for a Collocation-discretized problem.
"""
function CTSolvers.build_model(
    dm::CTSolvers.DiscretizedModel{<:Any,<:Collocation},
    initial_guess::CTModels.AbstractInitialGuess,
    modeler::CTSolvers.Modelers.ADNLP,
)::ADNLPModels.ADNLPModel

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

    # set adnlp backends
    if backend == :manual

        # build sparsity patterns for Jacobian and Hessian
        J_backend = ADNLPModels.SparseADJacobian(
            docp.dim_NLP_variables, f,
            docp.dim_NLP_constraints, c!,
            DOCP_Jacobian_pattern(docp),
        )
        H_backend = ADNLPModels.SparseReverseADHessian(
            docp.dim_NLP_variables, f,
            docp.dim_NLP_constraints, c!,
            DOCP_Hessian_pattern(docp),
        )
        backend_options = (
            gradient_backend=ADNLPModels.ReverseDiffADGradient,
            jacobian_backend=J_backend,
            hessian_backend=H_backend,
        )
    else
        # use backend preset
        backend_options = (backend=backend,)
    end

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

    return nlp
end

"""
$(TYPEDSIGNATURES)

Build an OCP solution from an ADNLP solver result for a Collocation-discretized problem.
"""
function CTSolvers.build_solution(
    dm::CTSolvers.DiscretizedModel{<:Any,<:Collocation},
    nlp_solution::SolverCore.AbstractExecutionStats,
    ::CTSolvers.Modelers.ADNLP,
)
    docp = dm.cache.docp

    # retrieve data from NLP solver
    objective, iterations, constraints_violation, message, status, successful = CTSolvers.extract_solver_infos(nlp_solution)

    # retrieve time grid  +++ put in build_OCP_solution
    T = get_time_grid(nlp_solution.solution, docp)

    # build OCP solution from NLP solution
    sol = build_OCP_solution(docp, nlp_solution, T,
    objective, iterations, constraints_violation, message, status, successful)

    return sol
end

# ==========================================================================================
# CTSolvers contract: Exa backend
# ==========================================================================================
"""
$(TYPEDSIGNATURES)

Build an `ExaModel` for a Collocation-discretized problem.

Mutates `dm.cache.exa_getter` so that `build_solution` can reuse it.
"""
function CTSolvers.build_model(
    dm::CTSolvers.DiscretizedModel{<:Any,<:Collocation},
    initial_guess::CTModels.AbstractInitialGuess,
    modeler::CTSolvers.Modelers.Exa,
)::ExaModels.ExaModel

    docp = dm.cache.docp
    ocp = dm.ocp

    BaseType = modeler[:base_type]
    backend = modeler[:backend]

    # recover discretization scheme and size
    scheme = Strategies.options(dm.discretizer)[:scheme]
    grid_size = Strategies.options(dm.discretizer)[:grid_size]

    # build initial guess
    functional_init = CTModels.build_initial_guess(ocp, initial_guess)
    x0 = __initial_guess(docp, functional_init)
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
    if m > 0
        control = hcat([x0[(n + 1 + i * (n + m)):(n + 1 + i * (n + m) + m - 1)]
        for i in 0:(N - 1)]...,)
        # see with JB: pass indeed to grid_size only for euler(_b), trapeze and midpoint
        control = [control control[:, end]]
    else
        # zero control dimension: create empty matrix with correct type
        control = similar(x0, 0, N + 1)
    end
    variable = x0[(end - q + 1):end]
    init = (variable, state, control)

    # build Exa model and getters
    # see with JB. later try to call Exa constructor here if possible, reusing existing functions...
    build_exa = CTModels.get_build_examodel(ocp)
    nlp, exa_getter = build_exa(;
        grid_size=grid_size,
        backend=backend,
        scheme=scheme,
        init=init,
        base_type=BaseType,
    )

    # store getter in the cache for build_solution -- explicit, no captured closure
    dm.cache.exa_getter = exa_getter

    return nlp
end

"""
$(TYPEDSIGNATURES)

Build an OCP solution from an Exa solver result for a Collocation-discretized problem.

# Precondition
`build_model` must have been called first (it populates `dm.cache.exa_getter`).
"""
function CTSolvers.build_solution(
    dm::CTSolvers.DiscretizedModel{<:Any,<:Collocation},
    nlp_solution::SolverCore.AbstractExecutionStats,
    ::CTSolvers.Modelers.Exa,
)
    docp = dm.cache.docp
    exa_getter = dm.cache.exa_getter

    isnothing(exa_getter) && throw(
        CTBase.Exceptions.IncorrectArgument(
            "build_solution (Exa): exa_getter is not set";
            got="dm.cache.exa_getter === nothing",
            expected="a getter populated by build_model",
            suggestion="Call build_model with an Exa modeler before build_solution",
            context="CTDirect.build_solution for Collocation/Exa",
        ),
    )

    # retrieve data from NLP solver
    objective, iterations, constraints_violation, message, status, successful = CTSolvers.extract_solver_infos(nlp_solution)

    # retrieve time grid
    T = get_time_grid_exa(nlp_solution, docp, exa_getter)

    # build OCP solution from NLP solution
    sol = build_OCP_solution(docp, nlp_solution, T,
    objective, iterations, constraints_violation, message, status, successful;
    exa_getter=exa_getter)

    return sol
end

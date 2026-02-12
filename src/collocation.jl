# ---------------------------------------------------------------------------
# Implementation of Collocation discretizer
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Collocation discretizer
# ---------------------------------------------------------------------------
struct Collocation <: AbstractOptimalControlDiscretizer
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
        description = "Time integration scheme (e.g., :midpoint, :trapeze)",
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

# default options for modelers backend (now in CTSolvers ?)
#__adnlp_backend() = :optimized
#__exa_backend() = nothing


# ==========================================================================================
# Build core DOCP structure with discretization information (ADNLP)
# ==========================================================================================
function get_docp(discretizer::Collocation, ocp::AbstractOptimalControlProblem)
    
    # recover discretization scheme and options
    scheme = Strategies.options(discretizer)[:scheme]
    grid_size = Strategies.options(discretizer)[:grid_size]
    time_grid = Strategies.options(discretizer)[:time_grid]

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
            # see with JB: pass indeed to grid_size only for euler(_b), trapeze and midpoint
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
    
    # NLP builder for ADNLPModels
    function build_adnlp_model(
        initial_guess::CTModels.AbstractOptimalControlInitialGuess;
        backend,
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
        if backend == :manual

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
            backend_options = (backend=backend,)
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
        initial_guess::CTModels.AbstractOptimalControlInitialGuess; 
        backend
    )::ExaModels.ExaModel where {BaseType<:AbstractFloat}

        # recover discretization scheme and size
        scheme = Strategies.options(discretizer)[:scheme]
        grid_size = Strategies.options(discretizer)[:grid_size]

        # build initial guess
        init = get_docp_initial_guess(:exa, docp, initial_guess)

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
        objective, iterations, constraints_violation, message, status, successful = CTSolvers.extract_solver_infos(nlp_solution, minimize)
  
        # retrieve time grid
        T = get_time_grid_exa(nlp_solution, docp, exa_getter)

        # build OCP solution from NLP solution
        sol = CTDirect.build_OCP_solution(docp, nlp_solution, T,
        objective, iterations, constraints_violation, message, status, successful; 
        exa_getter=exa_getter)
        
        return sol
    end

    #NB. it would be better to return builders as model/solution pairs since they are linked
    return CTSolvers.DiscretizedOptimalControlProblem(
        ocp,
        CTSolvers.ADNLPModelBuilder(build_adnlp_model),
        CTSolvers.ExaModelBuilder(build_exa_model),
        CTSolvers.ADNLPSolutionBuilder(build_adnlp_solution),
        CTSolvers.ExaSolutionBuilder(build_exa_solution),
    )
end

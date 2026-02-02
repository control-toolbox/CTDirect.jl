# ---------------------------------------------------------------------------
# Implementation of Collocation discretizer
# ---------------------------------------------------------------------------

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
    docp = DOCP(
        ocp;
        grid_size=grid_size,
        time_grid=time_grid,
        scheme=scheme,
    )

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

        ocp = docp.ocp

        # set initial guess data
        if (initial_guess === nothing)
            init = nothing
        else
            init = 
            (
                state=CTModels.state(initial_guess),
                control=CTModels.control(initial_guess),
                variable=CTModels.variable(initial_guess),
            )
        end

        # build functional initial guess
        functional_init = CTModels.build_initial_guess(ocp, init)

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

    # construct common data for builders
    discretizer.docp = get_docp(discretizer, ocp)

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

        # recover docp
        docp = discretizer.docp
        
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
                docp.dim_NLP_variables,
                f,
                docp.dim_NLP_constraints,
                c!,
                CTDirect.DOCP_Jacobian_pattern(docp),
            )
            H_backend = ADNLPModels.SparseReverseADHessian(
                docp.dim_NLP_variables,
                f,
                docp.dim_NLP_constraints,
                c!,
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
        
        docp = discretizer.docp

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

        # build initial guess (ADNLP format)
        docp = discretizer.docp
        init = get_docp_initial_guess(:exa, docp, initial_guess)

        # build Exa model and getters
        # +++ later try to call Exa constructor here if possible, reusing existing functions...
        build_exa = CTModels.get_build_examodel(ocp)
        nlp, discretizer.exa_getter = build_exa(;
            grid_size=grid_size,
            backend=exa_backend,
            scheme=scheme,
            init=init,
        )

        return nlp
    end

    # Solution builder for ExaModels
    function build_exa_solution(nlp_solution::SolverCore.AbstractExecutionStats)

        docp = discretizer.docp

        # retrieve data from NLP solver
        minimize = !docp.flags.max
        objective, iterations, constraints_violation, message, status, successful = CTModels.extract_solver_infos(nlp_solution, minimize)
  
        # retrieve time grid
        T = get_time_grid_exa(nlp_solution, docp, discretizer.exa_getter)

        # build OCP solution from NLP solution
        sol = CTDirect.build_OCP_solution(docp, nlp_solution, T,
        objective, iterations, constraints_violation, message, status, successful; 
        exa_getter=discretizer.exa_getter)
        
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

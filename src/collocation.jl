# ---------------------------------------------------------------------------
# Implementation of Collocation discretizer
# ---------------------------------------------------------------------------
function (discretizer::Collocation)(ocp::AbstractOptimalControlProblem)

    # ==========================================================================================
    # Scheme symbol mapping
    # ==========================================================================================
    SchemeSymbol = Dict(MidpointScheme => :midpoint, TrapezoidalScheme => :trapeze, TrapezeScheme => :trapeze)
    function scheme_symbol(discretizer::Collocation)
        scheme = CTModels.get_option_value(discretizer, :scheme)
        return SchemeSymbol[typeof(scheme)]
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
    function get_docp()
        
        # recover discretization scheme and options
        disc_method = scheme_symbol(discretizer)
        grid_size, time_grid = grid_options(discretizer)

        # initialize DOCP
        docp = DOCP(
            ocp;
            grid_size=grid_size,
            time_grid=time_grid,
            disc_method=disc_method,
        )

        # set bounds in DOCP
        variables_bounds!(docp)
        constraints_bounds!(docp)

        return docp
    end

    # ==========================================================================================
    # Build initial guess for discretized problem
    # ==========================================================================================
    function get_x0(
        initial_guess::Union{CTModels.AbstractOptimalControlInitialGuess,Nothing},
        docp
        )

        # build initial guess data
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
   
        #+++ add adjustment for exa case here, with input flag for modeler ?

        return x0
    end

    #+++get_x0_exa() see exa model builder below, call get_x0 and return init triplet

    # construct common data for builders
    discretizer.docp = get_docp()

    # ==========================================================================================
    # The needed builders for the construction of the final DiscretizedOptimalControlProblem
    # ==========================================================================================
    function build_adnlp_model(
        initial_guess::CTModels.AbstractOptimalControlInitialGuess;
        adnlp_backend=CTDirect.__adnlp_backend(),
        show_time=false,
        kwargs...
    )::ADNLPModels.ADNLPModel

        # build docp (to be renamed later as disc_core ?)
        docp = discretizer.docp
        
        # functions for objective and constraints
        f = x -> CTDirect.DOCP_objective(x, docp)
        c! = (c, x) -> CTDirect.DOCP_constraints!(c, x, docp)

        # build initial guess
        x0 = get_x0(initial_guess, docp)

        # set adnlp backends
        backend_options = (backend=adnlp_backend,)

        # build NLP
        nlp = ADNLPModel!(
            f,
            x0,
            docp.bounds.var_l,
            docp.bounds.var_u,
            c!,
            docp.bounds.con_l,
            docp.bounds.con_u;
            minimize=(!docp.flags.max),
            backend_options...,
            #unused_backends...,
            show_time=show_time,
        )

        return nlp
    end

    # Solution builder for ADNLPModels
    function build_adnlp_solution(nlp_solution::SolverCore.AbstractExecutionStats)
        
        # build docp (to be renamed later as disc_core ?)
        docp = discretizer.docp

        #retrieve data from NLP solver
        #objective, iterations, constraints_violation, message, status, successful = CTModels.extract_solver_infos(nlp_solution)
        objective, iterations, constraints_violation, message, status, successful = CTDirect.SolverInfos(nlp_solution)

        # retrieve time grid
        T = get_time_grid(nlp_solution.solution, docp)

        # build OCP solution from NLP solution
        sol = CTDirect.build_OCP_solution(docp, nlp_solution, objective, iterations, constraints_violation, message, status, successful, T)
        
        return sol
    end

    # NLP builder for ExaModels
    function build_exa_model(
        ::Type{BaseType}, 
        initial_guess::CTModels.AbstractOptimalControlInitialGuess; 
        exa_backend=CTDirect.__exa_backend(),
        kwargs...
    )::ExaModels.ExaModel where {BaseType<:AbstractFloat}

        # recover discretization scheme and options
        # since exa part does not reuse the docp struct
        disc_method = scheme_symbol(discretizer)
        grid_size, time_grid = grid_options(discretizer)

        # build initial guess (ADNLP format)
        docp = discretizer.docp
        x0 = get_x0(initial_guess, docp)

        #+++ use aux function here  get_x0_exa()
        # reshape initial guess for ExaModel variables layout
        # - do not broadcast, apparently fails on GPU arrays
        # - unused final control in examodel / euler, hence the different x0 sizes
        n = CTModels.state_dimension(ocp)
        m = CTModels.control_dimension(ocp)
        q = CTModels.variable_dimension(ocp)
        grid_size, time_grid = grid_options(discretizer)
        state = hcat([x0[(1 + i * (n + m)):(1 + i * (n + m) + n - 1)] for i in 0:grid_size]...) # grid_size + 1 states
        control = hcat(
        [
            x0[(n + 1 + i * (n + m)):(n + 1 + i * (n + m) + m - 1)] for
            i in 0:(grid_size - 1)
        ]...,
        ) # grid_size controls...
        # +++? todo: pass indeed to grid_size only for euler(_b), trapeze and midpoint
        control = [control control[:, end]] 
        variable = x0[(end - q + 1):end]

        # build Exa model and getters
        # +++ share
        build_exa = CTModels.get_build_examodel(ocp)
        nlp, exa_getter = build_exa(;
            grid_size=grid_size,
            backend=exa_backend,
            scheme=disc_method,
            init=(variable, state, control),
        )
        # remark: nlp.meta.x0[1:docp.dim_NLP_variables] = -vcat(state..., control..., variable) 
        # also work, and supersedes previous init via ExaModels start (itself overridden by init in solve)
    
        return nlp
    end

    # Solution builder for ExaModels
    function build_exa_solution(nlp_solution::SolverCore.AbstractExecutionStats)

        # +++ share
        disc_method = scheme_symbol(discretizer)
        grid_size, time_grid = grid_options(discretizer)
        build_exa = CTModels.get_build_examodel(ocp)
        nlp, exa_getter = build_exa(;
            grid_size=grid_size,
            scheme=disc_method,
        )

        #retrieve data from NLP solver
        #objective, iterations, constraints_violation, message, status, successful = CTModels.extract_solver_infos(nlp_solution)
        objective, iterations, constraints_violation, message, status, successful = CTDirect.SolverInfos(nlp_solution)

        # retrieve time grid
        docp = discretizer.docp
        T = get_time_grid_exa(nlp_solution, docp, exa_getter)

        sol = CTDirect.build_OCP_solution(docp, nlp_solution, objective, iterations, constraints_violation, message, status, successful, T; exa_getter)
        
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
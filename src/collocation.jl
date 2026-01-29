# ---------------------------------------------------------------------------
# Implementation of Collocation discretizer
# ---------------------------------------------------------------------------
__adnlp_backend() = :optimized
__exa_backend() = nothing

function (discretizer::Collocation)(ocp::AbstractOptimalControlProblem)

    # +++ some of these functions could be outside the discretizer ?

    # ==========================================================================================
    # Scheme symbol mapping
    # ==========================================================================================
    #=SchemeSymbol = Dict(MidpointScheme => :midpoint, TrapezoidalScheme => :trapeze, TrapezeScheme => :trapeze)
    function scheme_symbol(discretizer::Collocation)
        scheme = CTModels.get_option_value(discretizer, :scheme)
        return SchemeSymbol[typeof(scheme)]
    end=#
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
    function get_docp()
        
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
    function get_functional_init(initial_guess::Union{CTModels.AbstractOptimalControlInitialGuess,Nothing})
        
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

        # return functional initial guess
        return CTModels.build_initial_guess(ocp, init)
    end
    
    function get_x0(
        initial_guess::Union{CTModels.AbstractOptimalControlInitialGuess,Nothing},
        docp
        )

        # build functional initial guess
        functional_init = get_functional_init(initial_guess)

        # build discretized initial guess
        x0 = DOCP_initial_guess(docp, functional_init)
   
        return x0
    end

    #+++get_x0_exa() see exa model builder below, 
    #return init triplet  

    # construct common data for builders
    discretizer.docp = get_docp()

    # ==========================================================================================
    # The needed builders for the construction of the final DiscretizedOptimalControlProblem
    # ==========================================================================================
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
        x0 = get_x0(initial_guess, docp)

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
            x0,
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

        #retrieve data from NLP solver
        minimize = !docp.flags.max
        objective, iterations, constraints_violation, message, status, successful = CTModels.extract_solver_infos(nlp_solution, minimize)

        # build OCP solution from NLP solution
        sol = CTDirect.build_OCP_solution(docp, nlp_solution, objective, iterations, constraints_violation, message, status, successful)
        
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
        scheme = get_scheme(discretizer)
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
        # +++ later try to call Exa constructor here if possible, reusing existing functions...
        # +++ share
        build_exa = CTModels.get_build_examodel(ocp)
        nlp, discretizer.exa_getter = build_exa(;
            grid_size=grid_size,
            backend=exa_backend,
            scheme=scheme,
            init=(variable, state, control),
        )
        # remark: nlp.meta.x0[1:docp.dim_NLP_variables] = -vcat(state..., control..., variable) 
        # also work, and supersedes previous init via ExaModels start (itself overridden by init in solve)
    
        return nlp
    end

    # Solution builder for ExaModels
    function build_exa_solution(nlp_solution::SolverCore.AbstractExecutionStats)

        docp = discretizer.docp

        #retrieve data from NLP solver
        minimize = !docp.flags.max
        objective, iterations, constraints_violation, message, status, successful = CTModels.extract_solver_infos(nlp_solution, minimize)
  
        # build OCP solution from NLP solution
        sol = CTDirect.build_OCP_solution(docp, nlp_solution, objective, iterations, constraints_violation, message, status, successful; 
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
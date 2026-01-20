# ---------------------------------------------------------------------------
# Implementation of Collocation discretizer
# ---------------------------------------------------------------------------
function (discretizer::Collocation)(ocp::AbstractOptimalControlProblem)

    # ==========================================================================================
    # Scheme symbol mapping: just used for direct_transcription
    # ==========================================================================================
    SchemeSymbol = Dict(MidpointScheme => :midpoint, TrapezoidalScheme => :trapeze, TrapezeScheme => :trapeze)
    function scheme_symbol(discretizer::Collocation)
        scheme = CTModels.get_option_value(discretizer, :scheme)
        return SchemeSymbol[typeof(scheme)]
    end

    # ==========================================================================================
    # Build a DOCP from existing CTDirect.direct_transcription call
    # This will be replace by new implementations of build models / solutions
    # ==========================================================================================
    function get_docp(
        modeler::Symbol;
        kwargs...,
    )
        # recover discretization scheme
        disc_method = scheme_symbol(discretizer)

        # unified grid option: Int => grid_size, Vector => explicit time_grid
        grid = CTModels.get_option_value(discretizer, :grid)
        if grid isa Int
            grid_size = grid
            time_grid = nothing
        else
            grid_size = length(grid)
            time_grid = grid
        end

        #=if modeler == :exa && haskey(kwargs, :backend)
            # Route ExaModeler backend to CTDirect via exa_backend and drop backend from forwarded kwargs.
            exa_backend = kwargs[:backend]
            filtered_kwargs = (; (k => v for (k, v) in pairs(kwargs) if k != :backend)...)
            docp = CTDirect.direct_transcription(
                ocp,
                modeler;
                grid_size=grid_size,
                disc_method=scheme_ctdirect,
                init=init_ctdirect,
                time_grid=time_grid,
                exa_backend=exa_backend,
                filtered_kwargs...,
            )
        else
            docp = CTDirect.direct_transcription(
                ocp,
                modeler;
                grid_size=grid_size,
                disc_method=scheme_ctdirect,
                init=init_ctdirect,
                time_grid=time_grid,
                kwargs...,
            )
        end=#

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
        
    end


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
        docp = get_docp(:adnlp; kwargs...)
        
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

    function build_adnlp_solution(nlp_solution::SolverCore.AbstractExecutionStats)
        
        # build docp (to be renamed later as disc_core ?)
        docp = get_docp(:adnlp)

        #+++ retrieve data from NLP solver
        objective, iterations, constraints_violation, message, status, successful = CTDirect.SolverInfos(nlp_solution)

        # retrieve time grid
        T = get_time_grid(nlp_solution.solution, docp)

        # build OCP solution from NLP solution
        sol = CTDirect.build_OCP_solution(docp, nlp_solution, objective, iterations, constraints_violation, message, status, successful, T)
        
        return sol
    end

    function build_exa_model(
        ::Type{BaseType}, initial_guess::CTModels.AbstractOptimalControlInitialGuess; kwargs...
    )::ExaModels.ExaModel where {BaseType<:AbstractFloat}
        docp = get_docp(initial_guess, :exa; kwargs...)
        return CTDirect.nlp_model(docp)
    end

    function build_exa_solution(nlp_solution::SolverCore.AbstractExecutionStats)
        docp = get_docp(nothing, :exa)
        solu = CTDirect.build_OCP_solution(docp, nlp_solution)
        return solu
    end

    return CTModels.DiscretizedOptimalControlProblem(
        ocp,
        CTModels.ADNLPModelBuilder(build_adnlp_model),
        CTModels.ExaModelBuilder(build_exa_model),
        CTModels.ADNLPSolutionBuilder(build_adnlp_solution),
        CTModels.ExaSolutionBuilder(build_exa_solution),
    )
end
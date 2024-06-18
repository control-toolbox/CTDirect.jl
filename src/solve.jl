# CTDirect interface

# available methods by order of preference: from top to bottom
algorithmes = ()
algorithmes = add(algorithmes, (:adnlp, :ipopt))

"""
$(TYPEDSIGNATURES)

Return the list of available methods to solve the optimal control problem.
"""
function available_methods()::Tuple{Tuple{Vararg{Symbol}}}
    return algorithmes
end


"""
$(TYPEDSIGNATURES)

Discretize an optimal control problem into a nonlinear optimization problem (ie direct transcription)
"""
function directTranscription(ocp::OptimalControlModel,
    description...;
    init=OCPInit(),
    grid_size::Integer=__grid_size_direct(),
    time_grid=nothing)

    # build DOCP
    docp = DOCP(ocp, grid_size, time_grid)

    # set initial guess and bounds
    x0 = DOCP_initial_guess(docp, implicitInit(init))
    docp.var_l, docp.var_u = variables_bounds(docp)
    docp.con_l, docp.con_u = constraints_bounds(docp)

    # call NLP problem constructor
    docp.nlp = ADNLPModel!(x -> DOCP_objective(x, docp), 
                    x0,
                    docp.var_l, docp.var_u, 
                    (c, x) -> DOCP_constraints!(c, x, docp), 
                    docp.con_l, docp.con_u, 
                    backend = :optimized)

return docp

end


"""
$(TYPEDSIGNATURES)

Extract the NLP problem from the DOCP
"""
function getNLP(docp::DOCP)
    return docp.nlp
end


"""
$(TYPEDSIGNATURES)

Extract the NLP problem from the DOCP
"""
function setDOCPInit(docp::DOCP, init)

    nlp = getNLP(docp)
    nlp.meta.x0 .= DOCP_initial_guess(docp, implicitInit(init))

end


"""
$(TYPEDSIGNATURES)

Solve a discretized optimal control problem DOCP
"""
function solve(docp::DOCP;
    init=nothing,
    display::Bool=__display(),
    print_level::Integer=__print_level_ipopt(),
    mu_strategy::String=__mu_strategy_ipopt(),
    linear_solver::String=__linear_solver(),
    kwargs...)

    # Linear solver
    if (linear_solver == "ma27") || (linear_solver == "ma57") || (linear_solver == "ma77") || (linear_solver == "ma86") || (linear_solver == "ma97")
        if !LIBHSL_isfunctional()
            linear_solver = "mumps"
        end
    end
    if linear_solver == "spral"
        if !haskey(ENV, "OMP_CANCELLATION") || !haskey(ENV, "OMP_PROC_BIND")
            linear_solver = "mumps"
        end
    end

    # solve DOCP with NLP solver
    print_level = display ?  print_level : 0
    if init == nothing
        # use initial guess embedded in the DOCP
        docp_solution = ipopt(getNLP(docp), print_level=print_level, mu_strategy=mu_strategy, sb="yes",
                              linear_solver=linear_solver, kwargs...)
    else
        # use given initial guess
        docp_solution = ipopt(getNLP(docp), x0=DOCP_initial_guess(docp, implicitInit(init)), print_level=print_level,
                              mu_strategy=mu_strategy, sb="yes", linear_solver=linear_solver, kwargs...)
    end

    # return DOCP solution
    return docp_solution
end


"""
$(TYPEDSIGNATURES)

Solve an optimal control problem OCP by direct method
"""
function solve(ocp::OptimalControlModel,
    description...;
    init=OCPInit(),
    grid_size::Integer=__grid_size_direct(),
    time_grid=nothing,
    display::Bool=__display(),
    print_level::Integer=__print_level_ipopt(),
    mu_strategy::String=__mu_strategy_ipopt(),
    linear_solver::String=__linear_solver(),
    kwargs...)

    # build discretized OCP
    docp = directTranscription(ocp, description, init=implicitInit(init), grid_size=grid_size, time_grid=time_grid)

    # solve DOCP
    docp_solution = solve(docp, display=display, print_level=print_level, mu_strategy=mu_strategy, linear_solver=linear_solver, kwargs...)

    # build and return OCP solution
    return OCPSolutionFromDOCP(docp, docp_solution)
end


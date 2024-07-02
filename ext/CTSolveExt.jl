module CTSolveExt

using CTDirect
using CTBase
using CommonSolve
using DocStringExtensions

using NLPModelsIpopt
using HSL

"""
$(TYPEDSIGNATURES)

Solve a discretized optimal control problem DOCP
"""
function CommonSolve.solve(docp::DOCP;
    init=nothing,
    display::Bool=CTDirect.__display(),
    print_level::Integer=CTDirect.__print_level_ipopt(),
    mu_strategy::String=CTDirect.__mu_strategy_ipopt(),
    linear_solver::String=CTDirect.__linear_solver(),
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
    nlp = getNLP(docp)
    if isnothing(init)
        # use initial guess embedded in the DOCP
        docp_solution = ipopt(nlp, print_level=print_level, mu_strategy=mu_strategy, sb="yes", linear_solver=linear_solver; kwargs...)
    else
        # use given initial guess
        ocp = docp.ocp
        x0 = CTDirect.DOCP_initial_guess(docp, _OptimalControlInit(init, state_dim=ocp.state_dimension, control_dim=ocp.control_dimension, variable_dim=ocp.variable_dimension))

        docp_solution = ipopt(nlp, x0=x0, print_level=print_level, mu_strategy=mu_strategy, sb="yes", linear_solver=linear_solver; kwargs...)
    end

    # return DOCP solution
    return docp_solution
end


"""
$(TYPEDSIGNATURES)

Solve an optimal control problem OCP by direct method
"""
function CommonSolve.solve(ocp::OptimalControlModel,
    description...;
    init=nothing,
    grid_size::Integer=CTDirect.__grid_size_direct(),
    time_grid=nothing,
    display::Bool=CTDirect.__display(),
    print_level::Integer=CTDirect.__print_level_ipopt(),
    mu_strategy::String=CTDirect.__mu_strategy_ipopt(),
    linear_solver::String=CTDirect.__linear_solver(),
    kwargs...)

    # build discretized OCP
    docp = directTranscription(ocp, description, init=init, grid_size=grid_size, time_grid=time_grid)

    # solve DOCP (NB. init is already embedded in docp)
    docp_solution = solve(docp, display=display, print_level=print_level, mu_strategy=mu_strategy, linear_solver=linear_solver; kwargs...)

    # build and return OCP solution
    return OCPSolutionFromDOCP(docp, docp_solution)
end

end



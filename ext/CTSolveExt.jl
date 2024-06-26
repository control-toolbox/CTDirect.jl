module CTSolveExt

using CTDirect
using CTBase
using CommonSolve
using DocStringExtensions

using NLPModelsIpopt # solve
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
    if init == nothing
        # use initial guess embedded in the DOCP
        docp_solution = ipopt(getNLP(docp), print_level=print_level, mu_strategy=mu_strategy, sb="yes", linear_solver=linear_solver; kwargs...)
    else
        # use given initial guess
        docp_solution = ipopt(getNLP(docp), x0=CTDirect.DOCP_initial_guess(docp, OptimalControlInit(init)), print_level=print_level, mu_strategy=mu_strategy, sb="yes", linear_solver=linear_solver; kwargs...)
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
    init=OptimalControlInit(),
    grid_size::Integer=CTDirect.__grid_size_direct(),
    time_grid=nothing,
    display::Bool=CTDirect.__display(),
    print_level::Integer=CTDirect.__print_level_ipopt(),
    mu_strategy::String=CTDirect.__mu_strategy_ipopt(),
    linear_solver::String=CTDirect.__linear_solver(),
    kwargs...)

    # build discretized OCP
    docp = directTranscription(ocp, description, init=OptimalControlInit(init), grid_size=grid_size, time_grid=time_grid)

    # solve DOCP
    docp_solution = solve(docp, display=display, print_level=print_level, mu_strategy=mu_strategy, linear_solver=linear_solver; kwargs...)

    # build and return OCP solution
    return OCPSolutionFromDOCP(docp, docp_solution)
end

end
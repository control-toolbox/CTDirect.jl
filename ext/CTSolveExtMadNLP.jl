module CTSolveExtMadNLP

using CTDirect
using CTBase
using DocStringExtensions

using MadNLP
using HSL

"""
$(TYPEDSIGNATURES)

Solve a discretized optimal control problem DOCP
"""
function CTDirect.solve_docp(solver::MadNLPSolver, docp::DOCP;
    init=nothing,
    print_level::Integer=CTDirect.__print_level_ipopt(),
    mu_strategy::String=CTDirect.__mu_strategy_ipopt(),
    max_iter::Integer=CTDirect.__max_iter(),
    tol::Real=CTDirect.__tol(),
    linear_solver::String=CTDirect.__linear_solver(),
    kwargs...)

    # check linear solver requirements, default to mumps if needed
    # HSL
    if (linear_solver == "ma27") || (linear_solver == "ma57") || (linear_solver == "ma77") || (linear_solver == "ma86") || (linear_solver == "ma97")
        if !LIBHSL_isfunctional()
            println("WARNING: HSL not available, defaulting linear solver ",linear_solver, " to MUMPS") 
            linear_solver = "mumps"
        end
    end
    # SPRAL
    if linear_solver == "spral"
        if !haskey(ENV, "OMP_CANCELLATION") || !haskey(ENV, "OMP_PROC_BIND")
            println("WARNING: missing required environment variables for SPRAL (OMP_CANCELLATION=TRUE and OMP_PROC_BIND=TRUE), defaulting to MUMPS") 
            linear_solver = "mumps"
        end
    end

    # solve discretized problem with NLP solver
    if isnothing(init)
        # use initial guess embedded in the NLP
        #docp_solution = solve!(solver,
        #print_level=print_level, mu_strategy=mu_strategy, tol=tol, max_iter=max_iter, sb="yes", linear_solver=linear_solver; kwargs...)
        docp_solution = solve!(solver)
    else
        # build initial guess from provided data
        x0 = CTDirect.DOCP_initial_guess(docp, OptimalControlInit(init, state_dim=docp.dim_OCP_x, control_dim=docp.dim_NLP_u, variable_dim=docp.dim_NLP_v))

        # override initial guess embedded in the NLP
        #docp_solution = solve!(solver, x0=x0, print_level=print_level, mu_strategy=mu_strategy, tol=tol, max_iter=max_iter, sb="yes", linear_solver=linear_solver; kwargs...)
        docp_solution = solve!(solver, x0=x0)
    end

    # return DOCP solution
    return docp_solution
end


function CTDirect.build_solution(docp, docp_solution_madnlp::MadNLP.MadNLPExecutionStats)

    # could pass some status info too (get_status ?)
    solution = docp_solution_madnlp.solution

    # time grid
    N = docp.dim_NLP_steps
    T = zeros(N+1)
    for i=1:N+1
        T[i] = CTDirect.get_unnormalized_time(solution, docp, docp.NLP_normalized_time_grid[i])
    end

    # adjust objective sign for maximization problems
    if is_min(docp.ocp)
        objective = docp_solution_madnlp.objective
    else        
        objective = - docp_solution_madnlp.objective
    end

    # recover primal variables
    X, U, v = CTDirect.parse_DOCP_solution_primal(docp, solution)

    # recover costate
    P = CTDirect.parse_DOCP_solution_dual(docp, docp_solution_madnlp.multipliers)

    # build and return OCP solution
    return CTDirect.OCPSolutionFromDOCP_raw(docp, T, X, U, v, P,
    objective=objective, iterations=docp_solution_madnlp.iter,constraints_violation=docp_solution_madnlp.primal_feas, 
    message="MadNLP")
end


end



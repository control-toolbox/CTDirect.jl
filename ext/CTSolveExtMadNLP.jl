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
function CTDirect.solve_docp(tag::CTDirect.MadNLPTag, docp::DOCP, nlp;
    display::Bool=CTDirect.__display(),
    max_iter::Integer=CTDirect.__max_iterations(),
    tol::Real=CTDirect.__tolerance(),
    linear_solver::String=CTDirect.__madnlp_linear_solver(),
    kwargs...)

    # disable output if needed
    print_level = display ?  MadNLP.INFO : MadNLP.ERROR

    # preallocate solver (+++need to pass printlevel here ?)
    solver = MadNLPSolver(nlp, print_level=print_level)

    # solve discretized problem with NLP solver
    docp_solution = solve!(solver, tol=tol, max_iter=max_iter; kwargs...)

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



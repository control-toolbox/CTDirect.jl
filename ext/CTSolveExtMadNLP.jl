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
function CTDirect.solve_docp(
    solver_backend::CTDirect.MadNLPBackend,
    docp::CTDirect.DOCP,
    nlp;
    display::Bool = CTBase.__display(),
    max_iter::Integer = CTDirect.__max_iterations(),
    tol::Real = CTDirect.__tolerance(),
    linear_solver::String = CTDirect.__madnlp_linear_solver(),
    kwargs...,
)

    # disable output if needed
    print_level = display ? MadNLP.INFO : MadNLP.ERROR

    # preallocate solver (NB. need to pass printlevel here)
    solver = MadNLPSolver(nlp, print_level = print_level)

    # solve discretized problem with NLP solver
    docp_solution = solve!(solver, tol = tol, max_iter = max_iter; kwargs...)

    # return DOCP solution
    return docp_solution
end

function CTBase.OptimalControlSolution(docp, docp_solution::MadNLP.MadNLPExecutionStats)

    # adjust objective sign for maximization problems
    if is_min(docp.ocp)
        objective = docp_solution.objective
    else
        objective = -docp_solution.objective
    end

    # call lower level constructor
    return OptimalControlSolution(
        docp,
        primal = docp_solution.solution,
        dual = docp_solution.multipliers,
        objective = objective,
        iterations = docp_solution.iter,
        constraints_violation = docp_solution.primal_feas,
        message = "MadNLP",
        mult_LB = docp_solution.multipliers_L,
        mult_UB = docp_solution.multipliers_U,
    )
end

end

module CTDirectExtMadNLP

using CTDirect

using DocStringExtensions

using MadNLP
using HSL
using MKL

"""
$(TYPEDSIGNATURES)

Solve a discretized optimal control problem DOCP
"""
function CTDirect.solve_docp(
    solver_backend::CTDirect.MadNLPBackend,
    docp::CTDirect.DOCP;
    display::Bool=CTDirect.__display(),
    max_iter::Integer=CTDirect.__max_iterations(),
    tol::Real=CTDirect.__tolerance(),
    kwargs...,
)

    # disable output if needed
    print_level = display ? MadNLP.INFO : MadNLP.ERROR

    # retrieve NLP
    nlp = CTDirect.model(docp)

    # preallocate solver (NB. need to pass printlevel here)
    solver = MadNLPSolver(nlp; print_level=print_level, tol=tol, max_iter=max_iter, kwargs...)

    # solve discretized problem with NLP solver
    docp_solution = solve!(solver)

    # return DOCP solution
    return docp_solution
end

function CTDirect.SolverInfos(docp_solution::MadNLP.MadNLPExecutionStats)

    # info from SolverCore.GenericExecutionStats
    iterations = docp_solution.iter
    constraints_violation = docp_solution.primal_feas
    status = Symbol(docp_solution.status)
    successful = (status == :SOLVE_SUCCEEDED) || (status == :SOLVED_TO_ACCEPTABLE_LEVEL)

    return iterations, constraints_violation, "MadNLP", status, successful
end

end

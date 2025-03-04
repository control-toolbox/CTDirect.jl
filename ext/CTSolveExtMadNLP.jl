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
    linear_solver::String = CTDirect.__madnlp_linear_solver(), # ?
    kwargs...,
)

    # todo: add default print_level, pass kwargs properly

    # disable output if needed
    print_level = display ? MadNLP.INFO : MadNLP.ERROR

    # preallocate solver (NB. need to pass printlevel here)
    solver = MadNLPSolver(nlp, print_level=print_level, tol=tol, max_iter=max_iter)

    # solve discretized problem with NLP solver
    docp_solution = solve!(solver)

    # return DOCP solution
    return docp_solution
end

end

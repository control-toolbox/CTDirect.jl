module CTDirectExtMadNLP

using CTDirect

using DocStringExtensions

using MadNLP
using MadNLPMumps
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
    linear_solver=CTDirect.__madnlp_linear_solver(),
    kwargs...,
)

    # disable output if needed
    print_level = display ? MadNLP.INFO : MadNLP.ERROR

    # retrieve NLP
    nlp = CTDirect.nlp_model(docp)

    # set linear solver
    if linear_solver == "umfpack"
        linear_solver = MadNLP.UmfpackSolver
    elseif linear_solver == "mumps"
        linear_solver = MadNLPMumps.MumpsSolver
    end

    # preallocate solver (NB. need to pass printlevel here)
    solver = MadNLPSolver(nlp; print_level=print_level, tol=tol, max_iter=max_iter, linear_solver=linear_solver, kwargs...
    )

    # solve discretized problem with NLP solver
    docp_solution = solve!(solver)

    # return DOCP solution
    return docp_solution
end

function CTDirect.SolverInfos(nlp_solution::MadNLP.MadNLPExecutionStats)

    objective = nlp_solution.objective # NB sign is incorrct for max problems !
    iterations = nlp_solution.iter
    constraints_violation = nlp_solution.primal_feas
    status = Symbol(nlp_solution.status)
    successful = (status == :SOLVE_SUCCEEDED) || (status == :SOLVED_TO_ACCEPTABLE_LEVEL)

    return objective, iterations, constraints_violation, "MadNLP", status, successful
end

end

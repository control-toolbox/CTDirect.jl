module CTDirectExtMadNLP

using CTDirect

using DocStringExtensions

using MadNLP
using HSL
using MKL
using NLPModels

"""
$(TYPEDSIGNATURES)

Default value for MadNLP print level: `MadNLP.INFO`
Valid values are: MadNLP.{TRACE, DEBUG, INFO, NOTICE, WARN, ERROR}.
"""
__madnlp_print_level() = MadNLP.INFO

"""
$(TYPEDSIGNATURES)

Solve a discretized optimal control problem DOCP
"""
function CTDirect.solve_docp(
    ::CTDirect.MadNLPBackend,
    docp::CTDirect.DOCP;
    display::Bool=CTDirect.__display(),
    max_iter::Integer=CTDirect.__max_iterations(),
    tol::Real=CTDirect.__tolerance(),
    print_level=__madnlp_print_level(),
    kwargs...,
)

    # disable output if needed
    print_level = display ? print_level : MadNLP.ERROR

    # retrieve NLP
    nlp = CTDirect.nlp_model(docp)

    # preallocate solver (NB. need to pass printlevel here)
    solver = MadNLPSolver(
        nlp; print_level=print_level, tol=tol, max_iter=max_iter, kwargs...
    )

    # solve discretized problem with NLP solver
    nlp_solution = solve!(solver)

    # return NLP solution
    return nlp_solution
end

function CTDirect.SolverInfos(
    nlp_solution::MadNLP.MadNLPExecutionStats, nlp::NLPModels.AbstractNLPModel
)
    minimize = NLPModels.get_minimize(nlp)
    objective = minimize ? nlp_solution.objective : -nlp_solution.objective # sign depends on minimization for MadNLP
    iterations = nlp_solution.iter
    constraints_violation = nlp_solution.primal_feas
    status = Symbol(nlp_solution.status)
    successful = (status == :SOLVE_SUCCEEDED) || (status == :SOLVED_TO_ACCEPTABLE_LEVEL)
    return objective, iterations, constraints_violation, "MadNLP", status, successful
end

end

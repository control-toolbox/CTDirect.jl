module CTSolveExtKnitro

using CTDirect
using CTBase
using DocStringExtensions

using NLPModelsKnitro
using HSL

"""
$(TYPEDSIGNATURES)

Solve a discretized optimal control problem with Ipopt
"""
function CTDirect.solve_docp(
    solver_backend::CTDirect.KnitroBackend,
    docp::CTDirect.DOCP,
    nlp;
    display::Bool = CTBase.__display(),
    max_iter::Integer = CTDirect.__max_iterations(),
    tol::Real = CTDirect.__tolerance(),
    print_level::Integer = CTDirect.__ipopt_print_level(),
    kwargs...,
)

    solver = KnitroSolver(nlp)
    stats = solve!(solver, nlp)

    # return DOCP solution
    return stats
end

end

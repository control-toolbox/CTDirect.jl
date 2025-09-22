module CTDirectExtKnitro

using CTDirect

using DocStringExtensions

using NLPModelsKnitro
using HSL
using MKL

"""
$(TYPEDSIGNATURES)

Default value for Knitro print level: `3`
"""
__knitro_print_level() = 3

"""
$(TYPEDSIGNATURES)

Solve a discretized optimal control problem with Ipopt
"""
function CTDirect.solve_docp(
    solver_backend::CTDirect.KnitroBackend,
    docp::CTDirect.DOCP;
    display::Bool=CTDirect.__display(),
    max_iter::Integer=CTDirect.__max_iterations(),
    tol::Real=CTDirect.__tolerance(),
    print_level::Integer=__knitro_print_level(),
    kwargs...,
)

    # todo: pass kwargs properly

    # disable output if needed
    print_level = display ? print_level : 0

    # retrieve NLP
    nlp = CTDirect.nlp_model(docp)

    # preallocate solver (cannot seem to pass kwargs to solve! ...)
    solver = KnitroSolver(
        nlp; outlev=print_level, maxit=max_iter, feastol_abs=tol, opttol_abs=tol
    )

    # solve discretized problem with NLP solver
    docp_solution = solve!(solver, nlp)

    # return DOCP solution
    return docp_solution
end

end

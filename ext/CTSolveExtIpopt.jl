module CTSolveExtIpopt

using CTDirect
using CTBase
using DocStringExtensions

using NLPModelsIpopt
using HSL

"""
$(TYPEDSIGNATURES)

Solve a discretized optimal control problem DOCP
"""
#function CTDirect.solve_docp(solver::IpoptSolver, docp::DOCP, nlp;
function CTDirect.solve_docp(tag::CTDirect.IpoptTag, docp::DOCP, nlp;
    display::Bool=CTDirect.__display(),
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

    # override print_level if needed
    print_level = display ?  print_level : 0

    # preallocate solver
    solver = IpoptSolver(nlp)

    # solve discretized problem with NLP solver
    docp_solution = solve!(solver, nlp, print_level=print_level, mu_strategy=mu_strategy, tol=tol, max_iter=max_iter, sb="yes", linear_solver=linear_solver; kwargs...)

    # return DOCP solution
    return docp_solution
end

end



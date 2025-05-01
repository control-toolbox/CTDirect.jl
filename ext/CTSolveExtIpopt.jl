module CTSolveExtIpopt

using CTDirect

using DocStringExtensions

using NLPModelsIpopt
using HSL
using MKL

"""
$(TYPEDSIGNATURES)

Solve a discretized optimal control problem with Ipopt
"""
function CTDirect.solve_docp(
    solver_backend::CTDirect.IpoptBackend,
    docp::CTDirect.DOCP,
    nlp;
    display::Bool=CTDirect.__display(),
    max_iter::Integer=CTDirect.__max_iterations(),
    tol::Real=CTDirect.__tolerance(),
    print_level::Integer=CTDirect.__ipopt_print_level(),
    mu_strategy::String=CTDirect.__ipopt_mu_strategy(),
    linear_solver::String=CTDirect.__ipopt_linear_solver(),
    kwargs...,
)

    # check HSL requirements
    if linear_solver in ["ma27", "ma57", "ma77", "ma86", "ma97"] && !LIBHSL_isfunctional()
        println(
            "WARNING: HSL not available, defaulting given linear solver ",
            linear_solver,
            " to MUMPS",
        )
        linear_solver = "mumps"
    end

    # check SPRAL requirements
    if linear_solver == "spral" &&
       (!haskey(ENV, "OMP_CANCELLATION") || !haskey(ENV, "OMP_PROC_BIND"))
        println(
            "WARNING: missing required environment variables for SPRAL (OMP_CANCELLATION=TRUE and OMP_PROC_BIND=TRUE), defaulting to MUMPS",
        )
        linear_solver = "mumps"
    end

    # disable output if needed
    print_level = display ? print_level : 0

    # preallocate solver
    solver = IpoptSolver(nlp)

    # solve discretized problem with NLP solver
    docp_solution = solve!(
        solver,
        nlp,
        print_level=print_level,
        mu_strategy=mu_strategy,
        tol=tol,
        max_iter=max_iter,
        sb="yes",
        #check_derivatives_for_naninf = "yes",
        linear_solver=linear_solver;
        kwargs...,
    )

    # return DOCP solution
    return docp_solution
end

end

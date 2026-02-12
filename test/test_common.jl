# CT packages
using CTBase
using CTParser: CTParser, @def, @init
using CTModels: CTModels, objective, time_grid, iterations, state, control, variable, costate
using CTDirect
using CTSolvers

# other
using CommonSolve

# NLP modelers
using ADNLPModels
using ExaModels

# NLP solvers
using NLPModels
using NLPModelsIpopt
using MadNLPMumps

# misc
using SplitApplyCombine # for flatten in some tests

# solve a given OCP problem (given as struct)
# +++ make more similar to high level solve from Optimal Control (options ...)
function solve_problem(prob;
    max_iter=1000,
    tol=1e-6,
    modeler=:adnlp, #+++exa
    solver=:ipopt, #+++madnlp
    display=false,
    graph=false,
    init=nothing,
    adnlp_backend=:optimized,
    exa_backend=nothing,
    kwargs...)

    # discretized problem (model and solution builders)
    discretizer = CTDirect.Collocation(; kwargs...) # kwargs here
    docp = CTDirect.discretize(prob.ocp, discretizer)
    
    # initial guess for model builders
    if isnothing(init)
        my_init = CTModels.build_initial_guess(prob.ocp, prob.init)
    else
        my_init = CTModels.build_initial_guess(prob.ocp, init)
    end

    # NLP modeler
    if modeler == :adnlp
        my_modeler = CTSolvers.ADNLPModeler(; backend=adnlp_backend) # kwargs
    elseif modeler == :exa
        my_modeler = CTSolvers.ExaModeler(; backend=exa_backend) # kwargs
    else
        error("Unknown modeler: ", modeler)
    end

    # NLP solver
    if solver == :ipopt
        ipopt_options = Dict(
            :max_iter => max_iter,
            :tol => tol,
            :print_level => 5,
            :mu_strategy => "adaptive",
            :linear_solver => "mumps",
            :sb => "yes",
        )
        my_solver = CTSolvers.IpoptSolver(; ipopt_options...)
    elseif solver == :madnlp
        madnlp_options = Dict(
            :max_iter => max_iter,
            :tol => tol,
        )
        my_solver = CTSolvers.MadNLPSolver(; madnlp_options...)
    else
        error("Unknown solver: ", solver)
    end

    # solve DOCP and return OCP solution
    sol = CommonSolve.solve(docp, my_init, my_modeler, my_solver; display=display)

    return sol
end

# CT packages
using CTBase
using CTParser: CTParser, @def
using CTModels: CTModels, objective, time_grid, iterations
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

# check a specific example OCP.
# +++ use a solve moved from OptimalControl to CTSolvers ?
max_iter = 1000
tol = 1e-6
ipopt_options = Dict(
    :max_iter => max_iter,
    :tol => tol,
    :print_level => 3,
    :mu_strategy => "adaptive",
    :linear_solver => "Mumps",
    :sb => "yes",
)

function solve_problem(prob;
    modeler=:adnlp,
    solver=:ipopt,
    display=false,
    graph=false,
    kwargs...)

    discretizer = CTDirect.Collocation(; kwargs...) # kwargs here
    docp = CTDirect.discretize(prob.ocp, discretizer)
    init = CTModels.initial_guess(prob.ocp; prob.init...) # check if still needed

    if modeler == :adnlp
        my_modeler = CTModels.ADNLPModeler() # kwargs
    elseif modeler == :exa
        my_modeler = CTModels.ExaModeler() # kwargs
    else
        error("Unknown modeler: ", modeler)
    end

    if solver == :ipopt
        my_solver = CTSolvers.IpoptSolver(; ipopt_options...) # kwargs
    elseif solver == :madnlp
        my_solver = CTSolvers.MadNLPSolver() # kwargs
    else
        error("Unknown solver: ", solver)
    end

    sol = CommonSolve.solve(docp, init, my_modeler, my_solver; display=display)

    return sol
end
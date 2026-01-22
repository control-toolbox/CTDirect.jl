# CT packages
using CTBase
using CTParser: CTParser, @def
using CTModels
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
    modeler = CTModels.ADNLPModeler(),
    solver = CTSolvers.IpoptSolver(; ipopt_options...),
    display=false,
    graph=false)

    docp = CTDirect.discretize(prob.ocp, CTDirect.Collocation())
    init = CTModels.initial_guess(prob.ocp; prob.init...) # check if still needed
    sol = CommonSolve.solve(docp, init, modeler, solver; display=display)

    return sol
end
# runtests.jl
using Test

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
function check_problem(prob; display=false)
    modeler = CTModels.ADNLPModeler()
    solver = CTSolvers.IpoptSolver(; ipopt_options...)
    docp = CTDirect.discretize(prob.ocp, CTDirect.Collocation())
    init = CTModels.initial_guess(prob.ocp; prob.init...)
    sol = CommonSolve.solve(docp, init, modeler, solver; display=display)
    @test CTModels.successful(sol)
    @test CTModels.iterations(sol) <= max_iter
    @test CTModels.constraints_violation(sol) <= tol
    @test sol.objective ≈ prob.obj rtol = 1e-2
end

# check local test suite
macro ignore(e)
    :()
end

# # run either usual test suite on CPU, or GPU tests only 
# @testset verbose = true showtiming = true "Test CTDirect" begin
#     if "GPU" in ARGS
#         # ExaModels tests only (GPU on moonshot workflow)
#         include("./suite/test_exa.jl")
#     else
#         # CPU: run all scripts in subfolder suite/
#         include.(filter(contains(r".jl$"), readdir("./suite"; join=true)))
#     end
# end

# new tests
const VERBOSE = true
const SHOWTIMING = true
@testset verbose = VERBOSE showtiming = SHOWTIMING "New tests for CTDirect" begin
    #include("./ci/test_core_types.jl")
    #test_ctdirect_core_types()
    #include("./ci/test_discretization_api.jl")
    #test_ctdirect_discretization_api()
    #include("./ci/test_collocation.jl")
    #test_ctdirect_collocation()

    include("./ci/test_solve.jl")
    test_solve()
    
    #include("./ci/test_all_ocp.jl")
end

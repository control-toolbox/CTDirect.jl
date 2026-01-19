# runtests.jl
using Test

# OptimalControl
using CTBase
using CTParser: CTParser, @def
using CTModels:
    CTModels, objective, state, control, variable, costate, time_grid, iterations, criterion
using CTDirect:
    CTDirect,
    direct_transcription,
    build_OCP_solution,
    nlp_model,
    ocp_model

# activate NLP modelers
using ADNLPModels
using ExaModels

# activate NLP solvers
using NLPModels
using NLPModelsIpopt
using MadNLPMumps

# misc
using SplitApplyCombine # for flatten in some tests

# check a specific example
function check_problem(prob; kwargs...)
    sol = solve(prob.ocp; init=prob.init, kwargs...)
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
    include("./new/test_ctdirect_core_types.jl")
    test_ctdirect_core_types()
    include("./new/test_ctdirect_discretization_api.jl")
    test_ctdirect_discretization_api()
    include("./new/test_ctdirect_collocation_impl.jl")
    test_ctdirect_collocation_impl()
end
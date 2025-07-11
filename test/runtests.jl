# runtests.jl
using Test

# OptimalControl
using CTBase
using CTParser: CTParser, @def
using CTModels: CTModels, objective, state, control, variable, costate, time_grid, iterations
using CTDirect: CTDirect, solve, direct_transcription, set_initial_guess, build_OCP_solution

# activate NLP modelers
using ADNLPModels
# + using ExaModels (in test_exa for now)

# activate NLP solvers
using NLPModelsIpopt
using MadNLP

# misc
using SplitApplyCombine # for flatten in some tests

# check a specific example
function check_problem(prob; kwargs...)
    sol = solve(prob.ocp; init=prob.init, kwargs...)
    @test sol.objective ≈ prob.obj rtol = 1e-2
end

# check local test suite
macro ignore(e) :() end

# run either usual test suite on CPU, or GPU tests only 
@testset verbose = true showtiming = true "Test CTDirect" begin
    if "GPU" in ARGS
        # GPU tests only (moonshot workflow)
        include("test_gpu.jl")
    else
        # CPU: run all scripts in subfolder suite/
        include.(filter(contains(r".jl$"), readdir("./suite"; join=true)))
    end
end
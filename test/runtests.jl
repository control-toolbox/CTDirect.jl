using Test

using CTDirect: CTDirect, solve, direct_transcription, set_initial_guess, build_OCP_solution
using CTModels: CTModels, objective, state, control, variable, costate, time_grid, iterations
using CTParser: CTParser, @def, set_prefix
set_prefix(:CTModels) # tell CTParser def macro to use CTModels instead of OptimalControl

# NLP solvers
using NLPModelsIpopt
using MadNLP

# misc
using SplitApplyCombine # for flatten in some tests

# check local test suite
@testset verbose = true showtiming = true "Test suite" begin
    # run all scripts in subfolder suite/
    include.(filter(contains(r".jl$"), readdir("./suite"; join=true)))
end

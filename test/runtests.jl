# runtests.jl
using Test

# OptimalControl
using CTBase
using CTParser: CTParser, @def, prefix!, e_prefix!
using CTModels: CTModels, objective, state, control, variable, costate, time_grid, iterations
using CTDirect: CTDirect, solve, direct_transcription, set_initial_guess, build_OCP_solution, nlp
prefix!(:CTModels) # set CTParser def macro to use CTModels instead of OptimalControl
e_prefix!(:CTBase) # set CTParser def macro to use CTBase instead of OptimalControl (errors)

# activate NLP modelers
using ADNLPModels
# + using ExaModels (in test_exa for now)

# activate NLP solvers
using NLPModelsIpopt
using MadNLP

# misc
using SplitApplyCombine # for flatten in some tests

# check local test suite
macro ignore(e) :() end

@testset verbose = true showtiming = true "Test suite" begin
    # run all scripts in subfolder suite/
    include.(filter(contains(r".jl$"), readdir("./suite"; join=true)))
    #include("suite/test_exa.jl")
end

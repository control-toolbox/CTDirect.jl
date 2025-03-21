using Test

using CTDirect
using CTModels
import CTModels: objective, state, control, variable, costate, time_grid, iterations
import CTParser: @def, set_prefix # for abstract formulation 

using NLPModelsIpopt
using MadNLP
using SplitApplyCombine # for flatten in some tests

# tell CTParser def macro to use CTModels instead of OptimalControl
set_prefix(:CTModels)

# check local test suite
@testset verbose = true showtiming = true "Test suite" begin
    # run all scripts in subfolder suite/
    include.(filter(contains(r".jl$"), readdir("./suite"; join = true)))
end

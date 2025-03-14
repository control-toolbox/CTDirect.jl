using Test

using CTDirect
import CTBase
import CTModels

using NLPModelsIpopt
using MadNLP
using SplitApplyCombine # for flatten in some tests

# check local test suite
@testset verbose = true showtiming = true "Test suite" begin
    # run all scripts in subfolder suite/
    include.(filter(contains(r".jl$"), readdir("./suite"; join = true)))
end

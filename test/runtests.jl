using CTDirect
using Test
#using LinearAlgebra
#include("test_utils.jl")

# check local test suite
@testset verbose = true showtiming = true "All problems" begin
    # run all scripts in subfolder suite/
    include.(filter(contains(r".jl$"), readdir("./suite"; join=true)))  
end

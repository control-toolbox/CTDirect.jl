using Test
include("deps.jl")
# +++ later load all problems ?
include("problems/goddard.jl")

# check local test suite
@testset verbose = true showtiming = true "Test suite" begin
    # run all scripts in subfolder suite/
    include.(filter(contains(r".jl$"), readdir("./suite"; join=true)))  
end

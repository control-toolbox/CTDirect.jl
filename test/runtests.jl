using Test
include("common_deps.jl")
# +++ later load all examples ?
include("../examples/goddard.jl")

# check local test suite
@testset verbose = true showtiming = true "Test suite" begin
    # run all scripts in subfolder suite/
    include.(filter(contains(r".jl$"), readdir("./suite"; join=true)))  
end

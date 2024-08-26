using Test
include("deps.jl")

# check local test suite
# +++ add explicit list (cf OC) (keep suite folder)
@testset verbose = true showtiming = true "Test suite" begin
    # run all scripts in subfolder suite/
    include.(filter(contains(r".jl$"), readdir("./suite"; join = true)))
end

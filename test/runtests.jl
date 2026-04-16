# runtests.jl
using Test
include("test_common.jl")

function test_problem(prob; kwargs...)

    sol = solve_problem(prob; kwargs...)
    @test CTModels.successful(sol)
    @test sol.objective ≈ prob.obj rtol = 1e-2

end

# check local test suite
macro ignore(e)
    :()
end

const VERBOSE = true
const SHOWTIMING = true

# run either usual test suite on CPU, or GPU tests only 
@testset verbose = VERBOSE showtiming = SHOWTIMING "Test CTDirect" begin
    if "GPU" in ARGS
        # GPU tests only (specific workflow)
        include("./test_gpu.jl")
    else
        # CPU: run all scripts in subfolder ci/
        include.(filter(contains(r".jl$"), readdir("./ci"; join=true)))
    end
end

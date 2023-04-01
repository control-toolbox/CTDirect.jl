using CTDirect
using Test
using CTBase
# CTProblems
#using LinearAlgebra # +++?

include("test_utils.jl")

# test of optimal control problem from CTProblems
@testset verbose = true showtiming = true "Direct" begin
    for name in (
        "local_integrator",
        "local_goddard",
        )
        @testset "$name" begin
            include("test_$name.jl")
        end
    end
end

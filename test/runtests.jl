using CTDirect
using Test
using CTBase
using CTProblems
using LinearAlgebra # +++?

include("test_utils.jl")

prob = Problem(:integrator, :dim2, :energy); 
ocp = prob.model

# test of optimal control problem from CTProblems
@testset verbose = true showtiming = true "Direct" begin
    for name in (
        "integrator",
        "goddard",
        )
        @testset "$name" begin
            include("test_$name.jl")
        end
    end
end

#generate_coverage(; run_test = true)

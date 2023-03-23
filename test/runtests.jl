using CTDirect
using Test
using CTBase
using CTProblems

# direct_infos funtion tests
include("test_direct_infos.jl")
include("test_constraints.jl")

prob = Problem(:integrator, :dim2, :energy); 
ocp = prob.model

# test of optimal control problem from CTProblems
@testset verbose = true showtiming = true "Direct" begin
    for name in (
        #"integrator",
        #"goddard",
        )
        @testset "$name" begin
            include("test_$name.jl")
        end
    end
end

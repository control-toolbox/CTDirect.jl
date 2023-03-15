using CTDirect
using Test
using CTBase
using CTProblems

#
@testset verbose = true showtiming = true "Direct" begin
    for name in (
        #"integrator",
        "goddard",
        )
        @testset "$name" begin
            include("test_$name.jl")
        end
    end
end

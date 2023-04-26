using CTDirect
using Test
using CTBase
using CTProblems
using LinearAlgebra

include("test_utils.jl")

# test all problems in CTProblems (except consumption ones)
@testset verbose = true showtiming = true "All problems" begin

    problems_list = Problems(:(!:consumption))
    for prob in problems_list
        println("Test: ",prob.description)
        @testset "$(prob.description)" begin
            sol = solve(prob.model, grid_size=20, print_level=0)
            @test sol.objective ≈ prob.solution.objective rtol=1e-2
        end
    end        
end

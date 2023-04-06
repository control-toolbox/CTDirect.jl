using CTDirect
using Test
using CTBase
using CTProblems
using LinearAlgebra

include("test_utils.jl")

# test all problems in CTProblems
@testset verbose = true showtiming = true "All problems" begin
    problems_list = Problems()
    for problem_description in problems_list
        if :consumption in problem_description
            println("Skip: ", problem_description)
        else
            prob = Problem(problem_description)
            println("Test: ",problem_description)
            @testset "$(prob.description)" begin
                sol = solve(prob.model, grid_size=20, print_level=0)
                @test sol.objective ≈ prob.solution.objective rtol=1e-2
            end
        end
    end
end

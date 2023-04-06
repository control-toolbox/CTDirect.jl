using CTDirect
using Test
using CTBase
using CTProblems
using LinearAlgebra

include("test_utils.jl")

#=
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
=#

#generate_coverage(; run_test = true) ??

# test all problems in CTProblems
@testset verbose = true showtiming = true "All problems" begin
    problems_list = Problems()
    for problem_description in problems_list
        if :consumption ∈ problem_description
            nothing # smooth it?
        else
            prob = Problem(problem_description)
            println(prob.description)
            @testset "$(prob.description)" begin
                sol = solve(prob.model, grid_size=20, print_level=0)
                @test sol.objective ≈ prob.solution.objective rtol=1e-2
            end
        end
    end
end

println("Test: constraint types")

@testset verbose = true showtiming = true ":goddard :box_constraints" begin
    sol1 = solve(goddard, print_level=0, tol=1e-8)
    @test sol1.objective ≈ 1.0125 rtol=1e-2
end

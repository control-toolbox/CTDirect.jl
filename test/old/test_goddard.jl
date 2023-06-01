# goddard with state constraint - maximize altitude
println("Goddard test")
prob = Problem(:goddard, :all_constraint)
init = [1.01, 0.05, 0.8, 0.1]
sol = solve(prob.model, grid_size=10, print_level=0, init=init)
@testset verbose = true showtiming = true ":goddard :all_constraints" begin
    @test sol.objective ≈ prob.solution.objective rtol=1e-2
end

# goddard with state constraint - maximize altitude
println("Goddard test")
prob = Problem(:goddard, :state_constraint)
ocp = prob.model

init = [1.01, 0.25, 0.5, 0.4]
sol = solve(ocp, grid_size=10, print_level=0, init=init)

@test sol.objective ≈ -1.0 atol=1e-1
# @test constraints_violation(sol) < 1e-6 # n'existe pas pour une OptimalControlSolution

# goddard with state constraint - maximize altitude
println("Goddard test")
prob = Problem(:goddard, :state_constraint)
ocp = prob.model

# solve
println("Is solvable ? ", CTDirect.is_solvable(ocp))
init = [1.01, 0.25, 0.5, 0.4]
sol = solve(ocp, grid_size=10, print_level=0, init=init)

# check solution
@test sol.objective ≈ -1.0 atol=1e-1
# @test constraints_violation(sol) < 1e-6 # n'existe pas pour une OptimalControlSolution

remove_constraint!(ocp, :state_con2)
vmax = 0.1
constraint!(ocp, :state, 0, vmax, :control_con4)
sol = solve(ocp, grid_size=10, print_level=0, init=init)

# check solution
@test sol.objective ≈ -1.0 atol=1e-1


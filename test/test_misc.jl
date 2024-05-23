using CTDirect

println("Test: misc")

# simple integrator (control split as positive/negative parts)
ocp = Model()
state!(ocp, 1)
control!(ocp, 2)
time!(ocp, [0, 1])
constraint!(ocp, :initial, -1, :initial_constraint)
constraint!(ocp, :final, 0, :final_constraint)
constraint!(ocp, :control, 1:2, [0,0], [Inf, Inf], :positive_controls)
dynamics!(ocp, (x, u) -> -x - u[1] + u[2])
objective!(ocp, :lagrange, (x, u) -> (u[1]+u[2])^2)

println(available_methods())
println(is_solvable(ocp))

sol = solve(ocp, print_level=5)
println("Expected Objective 0.313, found ", sol.objective)
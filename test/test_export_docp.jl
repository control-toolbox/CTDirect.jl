using CTDirect
using CTBase

println("Test simple integrator: old formulation")
# min energy
ocp1 = Model()
state!(ocp1, 1)
control!(ocp1, 1)
time!(ocp1, [0, 1])
constraint!(ocp1, :initial, -1, :initial_constraint)
constraint!(ocp1, :final, 0, :final_constraint)
dynamics!(ocp1, (x, u) -> -x + u)
objective!(ocp1, :lagrange, (x, u) -> u^2)
sol1 = solve(ocp1, grid_size=100, print_level=0, tol=1e-12)
println("Expected Objective 0.313, found ", sol1.objective)

println("Test simple integrator: new formulation with export")
docp = DirectTranscription(ocp1, grid_size=100)
sol2 = solveDOCP(docp, print_level=0, tol=1e-12)
println("Expected Objective 0.313, found ", sol2.objective)


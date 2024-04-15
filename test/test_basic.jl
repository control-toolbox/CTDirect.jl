using CTDirect

# simple integrator min energy
ocp = Model()
state!(ocp, 1)
control!(ocp, 1)
time!(ocp, [0, 1])
constraint!(ocp, :initial, -1, :initial_constraint)
constraint!(ocp, :final, 0, :final_constraint)
dynamics!(ocp, (x, u) -> -x + u)
objective!(ocp, :lagrange, (x, u) -> u^2)

# all-in-one solve call
println("Test simple integrator: all in one solve call")
sol = solveDirect(ocp, grid_size=100, print_level=5, tol=1e-12)
println("Expected Objective 0.313, found ", sol.objective)

# split calls
println("Test simple integrator: split calls")
println("Direct transcription")
docp = directTranscription(ocp, grid_size=100)
nlp = getNLP(docp)
println("Solve discretized problem and retrieve solution")
sol = solveDOCP(docp, print_level=5, tol=1e-12)
println("Expected Objective 0.313, found ", sol.objective)

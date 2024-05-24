using CTDirect

println("Test: misc")

# simple integrator min energy
# control split as positive/negative parts for m=2 tets case
ocp = Model()
state!(ocp, 1)
control!(ocp, 2)
time!(ocp, [0, 1])
constraint!(ocp, :initial, -1, :initial_constraint)
constraint!(ocp, :final, 0, :final_constraint)
constraint!(ocp, :control, 1:2, [0,0], [Inf, Inf], :positive_controls)
dynamics!(ocp, (x, u) -> -x - u[1] + u[2])
objective!(ocp, :lagrange, (x, u) -> (u[1]+u[2])^2)

# all-in-one solve call
println(available_methods())
println(is_solvable(ocp))
println("Test simple integrator: all in one solve call")
sol = solve(ocp, grid_size=100, print_level=0, tol=1e-12)
println("Target 0.313, found ", sol.objective, " at ", sol.iterations, " iterations")

# split calls
println("Test simple integrator: split calls")
println("Direct transcription with default init")
docp = directTranscription(ocp, grid_size=100)
sol, dsol = solve(docp, print_level=0, tol=1e-12)
println("Target 0.313, found ", sol.objective, " at ", sol.iterations, " iterations")

# test NLP getter
nlp = getNLP(docp)

# warm start in directTranscription
println("Direct transcription with warm start (compact syntax)")
docp2 = directTranscription(ocp, grid_size=100, init=sol)
sol2, dsol2 = solve(docp2, print_level=5, tol=1e-12)

# test OCPSolutionFromDOCP_raw
println("\nRebuild OCP solution from raw vector")
sol3 = OCPSolutionFromDOCP_raw(docp2, dsol2.solution)

println("")
# 
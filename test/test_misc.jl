using CTDirect
using CTBase
using Plots

println("Test: misc")

# simple integrator min energy
# control split as positive/negative parts for m=2 tets case
ocp = Model()
state!(ocp, 1)
control!(ocp, 2)
time!(ocp, t0=0, tf=1)
constraint!(ocp, :initial, lb=-1, ub=-1)
constraint!(ocp, :final, lb=0, ub=0)
constraint!(ocp, :control, lb=[0,0], ub=[Inf, Inf])
dynamics!(ocp, (x, u) -> -x - u[1] + u[2])
objective!(ocp, :lagrange, (x, u) -> (u[1]+u[2])^2)
#=
@def ocp begin
    t ∈ [ 0, 1], time
    x ∈ R, state
    u ∈ R^2, control
    u(t) ≥ [0,0]
    x(0) == -1
    x(1) == 0 
    ẋ(t) == -x(t) - u₁(t) + u₂(t)
    int(u₁(t)+u₂(t))^2 → min
end
=#

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
dsol = solve(docp, print_level=0, tol=1e-12)
sol1 = OCPSolutionFromDOCP(docp, dsol)
println("Target 0.313, found ", sol1.objective, " at ", sol1.iterations, " iterations")

# test NLP getter
nlp = getNLP(docp)

# warm start in directTranscription
println("\nDirect transcription with warm start (compact syntax)")
docp2 = directTranscription(ocp, grid_size=100, init=sol)
dsol2 = solve(docp2, print_level=0, tol=1e-12)

# test OCPSolutionFromDOCP_raw
#println("\nRebuild OCP solution from raw vector")
#sol3 = OCPSolutionFromDOCP_raw(docp2, dsol2.solution)

# save / load solution in JLD2 format
save_OCP_solution(sol, filename_prefix="solution_test")
sol4 = load_OCP_solution("solution_test")
plot(sol4, show=true)
println("\nCheck JLD2 solution ", sol.objective == sol4.objective)

# export / read discrete solution in JSON format
# NB. we recover here a JSON Object...
export_OCP_solution(sol, filename_prefix="solution_test")
sol_disc_reloaded = read_OCP_solution("solution_test")
println("\nCheck JSON solution ", sol.objective == sol_disc_reloaded.objective)

# solve with explicit and non uniform time grid
println("\nTest explicit time grid") 
sol5 = solve(ocp, time_grid=LinRange(0,1,101), print_level=0, tol=1e-12)
println("Check default-like grid ", (sol5.objective == sol.objective) && (sol5.iterations == sol.iterations))
 
sol6 = solve(ocp, time_grid=[0,0.1,0.6,0.98,0.99,1], print_level=0, tol=1e-12)
println("Objective with small unbalanced grid ", sol6.objective)
plot(sol6, show=true)

println("")

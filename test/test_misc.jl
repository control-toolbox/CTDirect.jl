include("common_deps.jl")
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
sol = solve(ocp, print_level=0, tol=1e-12)
println("Target 0.313, found ", sol.objective, " at ", sol.iterations, " iterations")

# split calls
println("Test simple integrator: split calls")
println("Direct transcription with default init")
docp = direct_transcription(ocp)
dsol = solve(docp, print_level=0, tol=1e-12)
sol1 = ocp_solution_from_docp(docp, dsol)
println("Target 0.313, found ", sol1.objective, " at ", sol1.iterations, " iterations")

# test NLP getter
nlp = get_nlp(docp)

# warm start in directTranscription
println("\nDirect transcription with warm start (compact syntax)")
docp2 = direct_transcription(ocp, init=sol)
dsol2 = solve(docp2, print_level=0, tol=1e-12)

# test OCPSolutionFromNLP (no costate +++)
println("\nRebuild OCP solution from raw NLP solution")
sol3 = ocp_solution_from_nlp(docp2, dsol2.solution)
#plot(sol3, show=true)

# save / load solution in JLD2 format
save(sol, filename_prefix="solution_test")
sol4 = load("solution_test")
println("\nCheck JLD2 solution ", sol.objective == sol4.objective)

# export / read discrete solution in JSON format
# NB. we recover here a JSON Object...
export_ocp_solution(sol, filename_prefix="solution_test")
sol_disc_reloaded = import_ocp_solution("solution_test")
println("\nCheck JSON solution ", sol.objective == sol_disc_reloaded.objective)

println("")

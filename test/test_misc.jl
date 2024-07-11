include("deps.jl")
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
sol1 = build_solution(docp, dsol)
println("Target 0.313, found ", sol1.objective, " at ", sol1.iterations, " iterations")

# test NLP getter
nlp = get_nlp(docp)

# warm start in directTranscription
println("\nDirect transcription with warm start")
docp2 = direct_transcription(ocp, init=sol)
dsol2 = solve(docp2, print_level=0, tol=1e-12)

# build_solution from NLP solution (primal only)
println("Rebuild OCP solution from NLP solution (primal only)")
sol3 = build_solution(docp2, primal=dsol2.solution)
plot(sol3, show=true)

# build_solution from NLP solution (primal only)
println("Rebuild OCP solution from NLP solution (primal,dual)")
sol4 = build_solution(docp2, primal=dsol2.solution, dual=dsol2.multipliers)
plot(sol4, show=true)

# save / load solution in JLD2 format
println("\nSave / load solution in JLD2 format")
save(sol, filename_prefix="./test/solution_test")
sol_reloaded = load("./test/solution_test")
println("Check JLD2 solution ", sol.objective == sol_reloaded.objective)

# export / read discrete solution in JSON format
# NB. we recover here a JSON Object...
println("Export / import solution in JSON format")
export_ocp_solution(sol, filename_prefix="./test/solution_test")
sol_disc_reloaded = import_ocp_solution("./test/solution_test")
println("Check JSON solution ", sol.objective == sol_disc_reloaded.objective)

println("")

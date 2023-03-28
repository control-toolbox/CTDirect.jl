# double integrator - energy min
println("Double integrator test")
prob = Problem(:integrator, :dim1, :energy) 
ocp = prob.model
u_sol(t) = prob.solution.control(t)[1]

# solve
println("Is solvable ? ", CTDirect.is_solvable(ocp))
sol = solve(ocp, grid_size=100, print_level=0)

# check solution
u = t -> sol.control(t)[1]
T = sol.times
dT = T[2:end]-T[1:end-1]
N = length(T)
@test sum(dT .* [ abs(u(T[i])-u_sol(T[i])) for i ∈ 1:N-1] ) ≈ 0 atol=1e-1
@test sol.objective ≈ prob.solution.objective atol=1e-2
# @test constraints_violation(sol) < 1e-6 # ceci n'existe pas dans la OptimalControlSolution pour le moment

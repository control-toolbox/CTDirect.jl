# double integrator - energy min
prob = Problem(:integrator, :dim2, :energy); ocp = prob.model
u_sol(t) = prob.solution.control(t)[1]

# solve
sol = solve(ocp, print_level=5)

# solution
u = t -> control(sol)(t)[1]
T = time_steps(sol)
dT = T[2:end]-T[1:end-1]
N = length(T)

@test sum(dT .* [ abs(u(T[i])-u_sol(T[i])) for i ∈ 1:N-1] ) ≈ 0 atol=1e-1
@test objective(sol) ≈ 6.0 atol=1e-1 
# @test constraints_violation(sol) < 1e-6 # ceci n'existe pas dans la OptimalControlSolution pour le moment

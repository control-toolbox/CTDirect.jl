# double integrator - energy min
prob = Problem(:integrator, :dim2, :energy)
ocp = prob.model

# solve
sol = solve(ocp, print_level=5)

# solution
u_sol(t) = 6.0-12.0*t
U = sol.U
T = sol.T
dT = T[2:end]-T[1:end-1]

@test sum(dT .* abs.(U[1:end-1] - u_sol.(T[1:end-1]))) ≈ 0 atol=1e-1
@test objective(sol) ≈ 6.0 atol=1e-1
@test constraints_violation(sol) < 1e-6

using CTBase
using CTProblemLibrary
using CTDirect

# double integrator - energy min
prob = Problem(:integrator, :dim2, :energy)
ocp = prob.model

# initial guess (constant state and control functions)
init = [1., 0.5, 0.3]

# solve
#sol = solve(ocp, grid_size=10, print_level=5)
sol = solve(ocp, (:dummy,), grid_size=200, print_level=5, tol=1e-12, mu_strategy="adaptive", init=init)

# plot
plot(sol)

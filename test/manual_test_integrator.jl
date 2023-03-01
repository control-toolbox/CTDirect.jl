using CTDirect
using CTProblemLibrary
using CTBase # for plot

# double integrator - energy min
prob = Problem(:integrator, :dim2, :energy)
ocp = prob.model

# initial guess (constant state and control functions)
init = [1., 0.5, 0.3]

# solve
#sol = solve(ocp, grid_size=10, print_level=5)
sol = solve(ocp, grid_size=10, print_level=5, init=init)

# plot
plot(sol)

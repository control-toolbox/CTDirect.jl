using CTDirect
using CTBase
using CTProblems

# double integrator - energy min
prob = Problem(:integrator, :dim2, :energy)
ocp = prob.model

# initial guess (constant state and control functions)
init = [1.0, 0.5, 0.3]

# solve
#sol = solve(ocp, grid_size=10, print_level=5)
sol = solve(ocp, grid_size = 100, print_level = 5, tol = 1e-12, init = init)

p1 = plot(sol)

# control box
constraint!(ocp, :control, -4.01, 4.01, :control_con1)
sol = solve(ocp, grid_size = 100, print_level = 5, tol = 1e-12, init = init)

# plot
p2 = plot(sol)

plot(p1, p2)

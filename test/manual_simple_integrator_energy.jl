using CTDirect
using CTBase
using CTProblems


# double integrator - energy min

println(":exponential, :dim1, :energy problem")
prob = Problem(:exponential, :dim1, :energy)

ocp = prob.model

# initial guess (constant state and control functions)
#init = [1., 0.5, 0.3]

# solve
#sol = solve(ocp, grid_size=10, print_level=5)
sol = solve(ocp, grid_size=100, print_level=5, tol=1e-12, mu_strategy="adaptive", init=nothing)

p1 = plot(sol)


#= # control box
constraint!(ocp, :control, -4.01, 4.01, :control_con1)
sol = solve(ocp, grid_size=100, print_level=5, tol=1e-12, mu_strategy="adaptive", init=init)

# plot
p2 =  plot(sol) =#

plot(p1)
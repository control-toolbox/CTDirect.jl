using CTDirect
using CTProblems
using CTBase # for plot

# goddard with state constraint - maximize altitude
prob = Problem(:goddard, :state_constraint)
ocp = prob.model

# initial guess (constant state and control functions)
init = [1.01, 0.05, 0.8, 0.1]

#sol = solve(ocp, grid_size=20, print_level=5)
sol = solve(ocp, grid_size=100, print_level=5, tol=1e-12, mu_strategy="adaptive", init=init)

# plot
p1 = plot(sol)

#= remove_constraint!(ocp, :state_con2)
vmax = 0.1
constraint!(ocp, :state, Index(2), 0, vmax, :state_con4)
sol = solve(ocp, grid_size=100, print_level=5, init=init)

p2 = plot(sol) =#

plot(p1, size=(700,900))
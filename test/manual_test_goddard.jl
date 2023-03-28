using CTDirect
using CTProblems
using CTBase # for plot
using ADNLPModels

# goddard with state constraint - maximize altitude
prob = Problem(:goddard, :state_constraint)
ocp = prob.model

# initial guess (constant state and control functions)
init = [1.01, 0.05, 0.8, 0.1]

#sol = solve(ocp, grid_size=20, print_level=5)
sol = solve(ocp, grid_size=100, print_level=5, tol=1e-12, mu_strategy="adaptive", init=init)

# plot
plot(sol)

# goddard with state constraint - maximize altitude
prob = Problem(:goddard, :state_constraint)
ocp = prob.model

# initial guess (constant state and control functions)
init = [1.01, 0.25, 0.5, 0.4]

#sol = solve(ocp, grid_size=20, print_level=5)
sol = direct_solve(ocp, (:dummy,), grid_size=20, print_level=5, init=init)

plot(sol)

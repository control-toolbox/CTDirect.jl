using CTDirect
using CTBase

#=
ocp = Model()
state!(ocp, 1)
control!(ocp, 1)
time!(ocp, [0, 1])
constraint!(ocp, :initial, [-1], :initial_constraint)
constraint!(ocp, :control, [-1], [1], :control_constraint)
dynamics!(ocp, (x, u) ->  -x[1] + u[1])

println("dummy mayer problem to keep scalar state")
objective!(ocp, :mayer, (x0,xf) -> -xf)
sol = solve(ocp, grid_size=20, print_level=3)

println("min energy: vector syntax")
constraint!(ocp, :final, [0], :final_constraint)
objective!(ocp, :lagrange, (x, u) -> u[1]*u[1])
sol = solve(ocp, grid_size=20, print_level=3)

println("min energy: scalar syntax ")
remove_constraint!(ocp,:initial_constraint)
constraint!(ocp, :initial, -1, :initial_constraint)
remove_constraint!(ocp,:final_constraint)
constraint!(ocp, :final, 0, :final_constraint)
remove_constraint!(ocp,:control_constraint)
constraint!(ocp, :control, -1, 1, :control_constraint)
dynamics!(ocp, (x, u) ->  -x + u)
objective!(ocp, :lagrange, (x, u) -> u*u)
#display(ocp)
sol = solve(ocp, grid_size=20, print_level=3)

println("min energy: scalar x, vector u ")
remove_constraint!(ocp,:control_constraint)
constraint!(ocp, :control, [-1], [1], :control_constraint)
objective!(ocp, :lagrange, (x, u) -> u[1]*u[1])
dynamics!(ocp, (x, u) ->  -x + u[1])
sol = solve(ocp, grid_size=20, print_level=3)

println("min energy: scalar u, vector x ")
remove_constraint!(ocp,:initial_constraint)
constraint!(ocp, :initial, [-1], :initial_constraint)
remove_constraint!(ocp,:final_constraint)
constraint!(ocp, :final, [0], :final_constraint)
remove_constraint!(ocp,:control_constraint)
constraint!(ocp, :control, -1, 1, :control_constraint)
objective!(ocp, :lagrange, (x, u) -> u*u)
dynamics!(ocp, (x, u) ->  -x[1] + u)
sol = solve(ocp, grid_size=20, print_level=3)

println("min energy: weird mix")
remove_constraint!(ocp,:initial_constraint)
constraint!(ocp, :initial, [-1], :initial_constraint)
remove_constraint!(ocp,:final_constraint)
constraint!(ocp, :final, 0, :final_constraint)
remove_constraint!(ocp,:control_constraint)
constraint!(ocp, :control, [-1], 1, :control_constraint)
objective!(ocp, :lagrange, (x, u) -> u[1]*u[1])
dynamics!(ocp, (x, u) ->  -x[1] + u)
sol = solve(ocp, grid_size=20, print_level=3)

=#

# try energy with free tf in [0.5, 2]
ocp2 = Model(variable=true)
state!(ocp2, 1)
control!(ocp2, 1)
variable!(ocp2, 1)
time!(ocp2, 0, Index(1))
constraint!(ocp2, :initial, [-1], :initial_constraint)
constraint!(ocp2, :final, [0], :final_constraint)
constraint!(ocp2, :control, [-1], [1], :control_constraint)
constraint!(ocp2, :variable, 0.5, 2, :variable_constraint)
dynamics!(ocp2, (x, u, v) ->  -x[1] + u[1])
objective!(ocp2, :lagrange, (x, u, v) -> u[1]*u[1])
sol = solve(ocp2, grid_size=100, print_level=5, tol=1e-12)

# try min tf lagrange and mayer


plot(sol)


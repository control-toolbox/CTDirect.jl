using CTDirect
using CTBase

n=1
m=1
t0=0
tf=1
ocp = Model()
state!(ocp, n)
control!(ocp, m)
time!(ocp, [t0, tf])
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

plot(sol)


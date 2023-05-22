using CTDirect
using CTBase

n=1
m=1
t0=0
tf=1
x0=-1
xf=0
ocp = Model()
state!(ocp, n)
control!(ocp, m)
time!(ocp, [t0, tf])
constraint!(ocp, :initial, x0, :initial_constraint)
constraint!(ocp, :final, xf, :final_constraint)
constraint!(ocp, :control, -1, 1, :control_constraint)
dynamics!(ocp, (x, u) -> -x + u)
objective!(ocp, :lagrange, (x, u) -> u*u)

sol = solve(ocp, grid_size=20, print_level=5)
p1 = plot(sol)

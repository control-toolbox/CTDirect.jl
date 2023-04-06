using CTDirect
using CTBase
using CTProblems

#prob = Problem(:exponential, :dim1, :consumption)
#ocp = prob.model
#remove_constraint!(ocp,:control_constraint)
#constraint!(ocp, :control, 0, 1, :control_constraint2)

n=1
m=1
t0=0
tf=1
x0=-1
xf=0
ocp = Model()
state!(ocp, n)   # dimension of the state
control!(ocp, m) # dimension of the control
time!(ocp, [t0, tf])
constraint!(ocp, :initial, x0, :initial_constraint)
constraint!(ocp, :final, xf, :final_constraint)
constraint!(ocp, :control, -1, 1, :control_constraint)
#constraint!(ocp, :control, 0, 1, :control_constraint)
constraint!(ocp, :dynamics, (x, u) -> -x + u)
objective!(ocp, :lagrange, (x, u) -> abs(u))
#objective!(ocp, :lagrange, (x, u) -> sqrt(u*u))

sol = solve(ocp, grid_size=20, print_level=5)
p1 = plot(sol)

using CTDirect
using CTBase
using Plots

println("Test: double integrator with several objectives")

# min tf
ocp = Model(variable=true)
state!(ocp, 2)
control!(ocp, 1)
variable!(ocp, 1)
time!(ocp, t0=0, indf=1)
constraint!(ocp, :initial, lb=[0,0], ub=[0,0])
constraint!(ocp, :final, lb=[1,0], ub=[1,0])
constraint!(ocp, :control, lb=-1, ub=1)
constraint!(ocp, :variable, lb=0.1, ub=10)
dynamics!(ocp, (x, u, v) ->  [x[2], u])
objective!(ocp, :mayer, (x0, xf, v) -> v)
sol = solve(ocp, grid_size=100, print_level=0, tol=1e-12)
println("Target 2.0, found ", sol.objective)

# min tf (lagrange)
ocp = Model(variable=true)
state!(ocp, 2)
control!(ocp, 1)
variable!(ocp, 1)
time!(ocp, t0=0, indf=1)
constraint!(ocp, :initial, lb=[0,0], ub=[0,0])
constraint!(ocp, :final, lb=[1,0], ub=[1,0])
constraint!(ocp, :control, lb=-1, ub=1)
constraint!(ocp, :variable, lb=0.1, ub=10)
dynamics!(ocp, (x, u, v) ->  [x[2], u])
objective!(ocp, :lagrange, (x, u, v) -> 1)
sol = solve(ocp, grid_size=100, print_level=0, tol=1e-12)
println("Target 2.0, found ", sol.objective)

# max t0 (free t0 and tf)
ocp = Model(variable=true)
state!(ocp, 2)
control!(ocp, 1)
variable!(ocp, 2)
time!(ocp, ind0=1, indf=2)
constraint!(ocp, :initial, lb=[0,0], ub=[0,0])
constraint!(ocp, :final, lb=[1,0], ub=[1,0])
constraint!(ocp, :control, lb=-1, ub=1)
constraint!(ocp, :variable, lb=[0.1,0.1], ub=[10,10])
constraint!(ocp, :variable, f=v->v[2]-v[1], lb=0.1, ub=Inf)
dynamics!(ocp, (x, u, v) ->  [x[2], u])
objective!(ocp, :mayer, (x0, xf, v) -> v[1], :max)
sol = solve(ocp, grid_size=100, print_level=0, tol=1e-12)
println("Target 8.0, found ", sol.objective)

# with explicit grid
sol = solve(ocp, time_grid=LinRange(0,1,101), print_level=0, tol=1e-12)
println("Target 8.0, found ", sol.objective)

# with non uniform grid
sol = solve(ocp, time_grid=[0,0.1,0.3,0.5,0.6,0.8,0.95,1], print_level=0)
plot(sol, show=true)
println("Target 8.0, coarse grid ", sol.objective)

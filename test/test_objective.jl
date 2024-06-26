include("common_deps.jl")

println("Test: objective")

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
sol = solve(ocp, print_level=0, tol=1e-12)
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
sol = solve(ocp, print_level=0, tol=1e-12)
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
sol = solve(ocp, print_level=0, tol=1e-12)
println("Target 8.0, found ", sol.objective)


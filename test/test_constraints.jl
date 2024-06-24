using CTDirect
using CTBase
using Plots

# goddard max final altitude (all constraint types formulation)
println("Test: all constraint types")

# general definitions
Cd = 310
Tmax = 3.5
β = 500
b = 2
r0 = 1
v0 = 0
vmax = 0.1
m0 = 1
mf = 0.6
x0 = [ r0, v0, m0 ]
function F0(x)
    r, v, m = x
    D = Cd * v^2 * exp(-β*(r - 1))
    return [ v, -D/m - 1/r^2, 0 ]
end
function F1(x)
    r, v, m = x
    return [ 0, Tmax/m, -b*Tmax ]
end

# box constraints
ocp = Model(variable=true)
state!(ocp, 3)
control!(ocp, 1)
variable!(ocp, 1)
time!(ocp, t0=0, indf=1)
constraint!(ocp, :initial, lb=x0, ub=x0)
constraint!(ocp, :final, rg=3, lb=mf, ub=Inf)
constraint!(ocp, :state, lb=[r0,v0,mf], ub=[r0+0.2,vmax,m0])
constraint!(ocp, :control, lb=0, ub=1)
constraint!(ocp, :variable, lb=0.01, ub=Inf)
objective!(ocp, :mayer, (x0, xf, v) -> xf[1], :max)
dynamics!(ocp, (x, u, v) -> F0(x) + u*F1(x) )

sol = solve(ocp, grid_size=100, print_level=0, tol=1e-8)
println("Target 1.0125, found ", sol.objective, " at ", sol.iterations, " iterations")

# functional constraints
ocp1 = Model(variable=true)
state!(ocp1, 3)
control!(ocp1, 1)
variable!(ocp1, 1)
time!(ocp1, t0=0, indf=1)
constraint!(ocp1, :initial, lb=x0, ub=x0)
constraint!(ocp1, :final, rg=3, lb=mf, ub=Inf)
constraint!(ocp1, :state, f=(x,v)->x, lb=[r0,v0,mf], ub=[r0+0.2,vmax,m0]) # fail !
constraint!(ocp1, :control, f=(u,v)->u, lb=0, ub=1)
constraint!(ocp1, :variable, f=v->v, lb=0.01, ub=Inf)
objective!(ocp1, :mayer, (x0, xf, v) -> xf[1], :max)
dynamics!(ocp1, (x, u, v) -> F0(x) + u*F1(x) )

sol = solve(ocp1, grid_size=100, print_level=5, tol=1e-8)
println("Target 1.0125, found ", sol.objective, " at ", sol.iterations, " iterations")
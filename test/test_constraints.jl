using CTDirect
using CTBase
using Plots

println("Test: all constraint types")
# +++ todo: complete with different constraint formulations
# goddard max final altitude (all constraint types formulation)
ocp = Model(variable=true)
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
function FF0(x)
    r, v, m = x
    D = Cd * v^2 * exp(-β*(r - 1))
    return [ v, -D/m - 1/r^2, 0 ]
end
function FF1(x)
    r, v, m = x
    return [ 0, Tmax/m, -b*Tmax ]
end
dynamics!(ocp, (x, u, v) -> FF0(x) + u*FF1(x) )

sol = solve(ocp, grid_size=100, print_level=0, tol=1e-8)
println("Target 1.0125, found ", sol.objective, " at ", sol.iterations, " iterations")

# explicit grid
sol1 = solve(ocp, time_grid=LinRange(0,1,101), print_level=0, tol=1e-8)
println((sol1.objective==sol.objective) && (sol1.iterations==sol.iterations))

# non uniform grid
sol2 = solve(ocp, time_grid=[0,0.1,0.6,0.98,0.99,1], print_level=0, tol=1e-8)
println("Objective with small unbalanced grid ", sol2.objective)
plot(sol2, show=true)
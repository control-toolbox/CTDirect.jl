using CTDirect
using CTBase
using BenchmarkTools

# Goddard with same formulation as bocop3
# max altitude, speed limit, only box constraints
n = 3
m = 1
Cd = 310
Tmax = 3.5
β = 500
b = 2
t0 = 0
r0 = 1
v0 = 0
vmax = 0.1
m0 = 1
mf = 0.6
x0 = [ r0, v0, m0 ]
ocp = Model()
time!(ocp, :initial, t0)
state!(ocp, 3, ["r", "v", "m"])
control!(ocp, 1)
constraint!(ocp, :initial, x0, :initial_constraint)
constraint!(ocp, :state, 1:3, [1,0,0.6], [Inf,0.1,1], :box_state)
constraint!(ocp, :control, Index(1), 0, 1, :box_control)
objective!(ocp, :mayer,  (t0, x0, tf, xf) -> xf[1], :max)
function F0(x)
    r, v, m = x
    D = Cd * v^2 * exp(-β*(r - 1))
    return [ v, -D/m - 1/r^2, 0 ]
end
function F1(x)
    r, v, m = x
    return [ 0, Tmax/m, -b*Tmax ]
end
f(x, u) = F0(x) + u*F1(x)
constraint!(ocp, :dynamics, f)

# settings
init = [1.01, 0.05, 0.8, 0.1]
N = 10

# dummy run then 2 runs
@time sol = solve(ocp, grid_size=N, print_level=0, tol=1e-12, mu_strategy="adaptive", init=init)
println("N=",N, " Obj ",sol.objective," Iter ", sol.iterations)
@time sol = solve(ocp, grid_size=N, print_level=0, tol=1e-12, mu_strategy="adaptive", init=init)
@time sol = solve(ocp, grid_size=N, print_level=0, tol=1e-12, mu_strategy="adaptive", init=init)

# benchmark. Not very practical -_-
#BenchmarkTools.DEFAULT_PARAMETERS.samples = 2
#@benchmark solve(ocp, grid_size=N, print_level=0, tol=1e-12, mu_strategy="adaptive", init=init)

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

# dummy run
init = [1.01, 0.05, 0.8, 0.1]
solve(ocp, grid_size=10, print_level=0, tol=1e-12, mu_strategy="adaptive", init=init)

# benchmark for different grid sizes
BenchmarkTools.DEFAULT_PARAMETERS.samples = 2
println("Benchmark N=10")
@btime sol = solve(ocp, grid_size=10, print_level=0, tol=1e-12, mu_strategy="adaptive", init=init)
println("Obj ",sol.objective," Iter ", sol.iterations)

println("Benchmark N=50")
@btime sol = solve(ocp, grid_size=50, print_level=0, tol=1e-12, mu_strategy="adaptive", init=init)
println("Obj ",sol.objective," Iter ", sol.iterations)

println("Benchmark N=100")
@btime sol = solve(ocp, grid_size=100, print_level=0, tol=1e-12, mu_strategy="adaptive", init=init)
println("Obj ",sol.objective," Iter ", sol.iterations)

println("Benchmark N=150")
@btime sol = solve(ocp, grid_size=150, print_level=0, tol=1e-12, mu_strategy="adaptive", init=init)
println("Obj ",sol.objective," Iter ", sol.iterations)

println("Benchmark N=200")
@btime sol = solve(ocp, grid_size=200, print_level=0, tol=1e-12, mu_strategy="adaptive", init=init)
println("Obj ",sol.objective," Iter ", sol.iterations)
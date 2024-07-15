include("deps.jl")
include("../problems/goddard.jl")

println("Test: constraint types")


# goddard, box constraints
ocp = goddard
sol = solve(ocp, print_level=0, tol=1e-8)
println("Target 1.0125, found ", sol.objective, " at ", sol.iterations, " iterations")

# functional constraints
ocp1 = Model(variable=true)
state!(ocp1, 3)
control!(ocp1, 1)
variable!(ocp1, 1)
time!(ocp1, t0=0, indf=1)
constraint!(ocp1, :initial, lb=x0, ub=x0)
constraint!(ocp1, :final, rg=3, lb=mf, ub=Inf)
constraint!(ocp1, :state, f=(x,v)->x, lb=[r0,v0,mf], ub=[r0+0.2,vmax,m0])
constraint!(ocp1, :control, f=(u,v)->u, lb=0, ub=1)
constraint!(ocp1, :variable, f=v->v, lb=0.01, ub=Inf)
objective!(ocp1, :mayer, (x0, xf, v) -> xf[1], :max)
dynamics!(ocp1, (x, u, v) -> F0(x) + u*F1(x) )

# note: the equations do not handle r<1 well
# without the box constraint on x, the default init (0.1) is not suitable
sol = solve(ocp1, print_level=0, tol=1e-8, init=(state=[1.01,0.05,0.75],))
println("Target 1.0125, found ", sol.objective, " at ", sol.iterations, " iterations")

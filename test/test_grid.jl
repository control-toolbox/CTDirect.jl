include("common_deps.jl")
using Plots

println("Test: grid options")

# simple integrator min energy
# control split as positive/negative parts for m=2 tets case
ocp = Model()
state!(ocp, 1)
control!(ocp, 2)
time!(ocp, t0=0, tf=1)
constraint!(ocp, :initial, lb=-1, ub=-1)
constraint!(ocp, :final, lb=0, ub=0)
constraint!(ocp, :control, lb=[0,0], ub=[Inf, Inf])
dynamics!(ocp, (x, u) -> -x - u[1] + u[2])
objective!(ocp, :lagrange, (x, u) -> (u[1]+u[2])^2)
sol0 = solve(ocp, print_level=0)

# solve with explicit and non uniform time grid
time_grid = LinRange(0,1,CTDirect.__grid_size_direct()+1)
sol5 = solve(ocp, time_grid=time_grid, print_level=0)
println((sol5.objective == sol0.objective) && (sol5.iterations == sol0.iterations))

time_grid = [0,0.1,0.3,0.6,0.98,0.99,1]
sol6 = solve(ocp, time_grid=time_grid, print_level=0)
println(sol6.objective)

# 2. integrator free times
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

sol = solve(ocp, time_grid=LinRange(0,1,CTDirect.__grid_size_direct()+1), print_level=0, tol=1e-12)
println("Target 8.0, found ", sol.objective)

sol = solve(ocp, time_grid=[0,0.1,0.3,0.5,0.6,0.8,0.95,1], print_level=0)
plot(sol, show=true)
println("Target 8.0, coarse grid ", sol.objective)


# 3. parametric ocp
function ocp_T(T)
    @def ocp begin
        t ∈ [ 0, T ], time
        x ∈ R², state
        u ∈ R, control
        q = x₁
        v = x₂
        q(0) == 0
        v(0) == 0
        q(T) == 1
        v(T) == 0
        ẋ(t) == [ v(t), u(t) ]
        ∫(u(t)^2) → min
    end
    return ocp
end

ocpT2 = ocp_T(2)
solT2 = solve(ocpT2, print_level=0)

solT2_exp = solve(ocpT2, time_grid=LinRange(0,1,CTDirect.__grid_size_direct()+1),print_level=0)
println("T=2 Check explicit grid ", (solT2.objective==solT2_exp.objective) && (solT2.iterations==solT2_exp.iterations))

solT2_nonunif = solve(ocpT2, time_grid=[0,0.3,1,1.9,2],print_level=0)
println("T=2 with non-uniform grid ", solT2_nonunif.objective)
plot(solT2_nonunif, show=true)

# +++4 pseudo grid refinement with manual grid input 
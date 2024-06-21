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
constraint!(ocp, :final, rg=3, lb=mf, ub=mf)
constraint!(ocp, :state, lb=[r0,v0,mf], ub=[r0+0.2,vmax,m0])
constraint!(ocp, :control, lb=0, ub=1)
constraint!(ocp, :variable, lb=0.01, ub=Inf)
objective!(ocp, :mayer,  (x0, xf, v) -> xf[1], :max)
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

@testset verbose = true showtiming = true ":goddard :max_rf :all_constraints" begin
    sol1 = solve(ocp, grid_size=100, print_level=0, tol=1e-8)
    @test sol1.objective ≈ 1.0125 rtol=1e-2
end

#=
# explicit grid
@testset verbose = true showtiming = true ":explicit_grid" begin
    sol = solve(ocp, grid_size=100, print_level=0, tol=1e-8)
    sol1 = solve(ocp, time_grid=LinRange(0,1,101), print_level=0, tol=1e-8)
    @test (sol1.objective==sol.objective) && (sol1.iterations==sol.iterations)
end

# non uniform grid
@testset verbose = true showtiming = true ":non_uniform_grid :free_tf" begin
    sol2 = solve(ocp, time_grid=[0,0.1,0.6,0.98,0.99,1], print_level=0, tol=1e-8)
    @test sol2.objective ≈ 1.0094 rtol=1e-2
end
=#
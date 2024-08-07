println("Test: grid options")

# 1. simple integrator min energy (dual control for test)
ocp = Model()
state!(ocp, 1)
control!(ocp, 2)
time!(ocp, t0=0, tf=1)
constraint!(ocp, :initial, lb=-1, ub=-1)
constraint!(ocp, :final, lb=0, ub=0)
constraint!(ocp, :control, lb=[0,0], ub=[Inf, Inf])
dynamics!(ocp, (x, u) -> -x - u[1] + u[2])
objective!(ocp, :lagrange, (x, u) -> (u[1]+u[2])^2)
sol0 = solve(ocp, display=false)

# solve with explicit and non uniform time grid
@testset verbose = true showtiming = true ":explicit_grid" begin
    time_grid = LinRange(0,1,CTDirect.__grid_size()+1)
    sol = solve(ocp, time_grid=time_grid, display=false)
    @test (sol.objective == sol0.objective) && (sol.iterations == sol0.iterations)
end

@testset verbose = true showtiming = true ":non_uniform_grid" begin
    time_grid = [0,0.1,0.3,0.6,0.98,0.99,1]
    sol = solve(ocp, time_grid=time_grid, display=false)
    @test sol.objective ≈ 0.309 rtol=1e-2
end


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
sol0 = solve(ocp, display=false)

@testset verbose = true showtiming = true ":explicit_grid" begin
    sol = solve(ocp, time_grid=LinRange(0,1,CTDirect.__grid_size()+1), display=false)
    @test (sol.objective == sol0.objective) && (sol.iterations == sol0.iterations)
end

@testset verbose = true showtiming = true ":non_uniform_grid" begin
    sol = solve(ocp, time_grid=[0,0.1,0.3,0.5,0.6,0.8,0.95,1],  display=false)
    #plot(sol, show=true)
    @test sol.objective ≈ 7.96 rtol=1e-2
end

# 3. parametric ocp (T=2) with explicit / non-uniform grid
@def ocpT2 begin
    t ∈ [ 0, 2 ], time
    x ∈ R², state
    u ∈ R, control
    q = x₁
    v = x₂
    q(0) == 0
    v(0) == 0
    q(2) == 1
    v(2) == 0
    ẋ(t) == [ v(t), u(t) ]
    ∫(u(t)^2) → min
end
sol0 = solve(ocpT2, display=false)

@testset verbose = true showtiming = true ":explicit_grid" begin
    sol = solve(ocpT2, time_grid=LinRange(0,1,CTDirect.__grid_size()+1), display=false)
    @test (sol.objective == sol0.objective) && (sol.iterations == sol0.iterations)
end

@testset verbose = true showtiming = true ":non_uniform_grid" begin
    sol = solve(ocpT2, time_grid=[0,0.3,1,1.9,2], display=false)
    @test sol.objective ≈ 2.43 rtol=1e-2
end

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

@testset verbose = true showtiming = true ":max_t0 :explicit_grid" begin
    sol = solve(ocp, time_grid=LinRange(0,1,CTDirect.__grid_size()+1), display=false)
    @test sol.objective ≈ 8.0 rtol=1e-2
end

@testset verbose = true showtiming = true ":max_t0 :non_uniform_grid" begin
    sol = solve(ocp, time_grid=[0,0.1,0.6,0.95,1],  display=false)
    @test sol.objective ≈ 7.48 rtol=1e-2
end
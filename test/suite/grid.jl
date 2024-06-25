
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
@testset verbose = true showtiming = true ":explicit_grid" begin
    time_grid = LinRange(0,1,CTDirect.__grid_size_direct()+1)
    sol5 = solve(ocp, time_grid=time_grid, print_level=0)
    @test (sol5.objective == sol0.objective) && (sol5.iterations == sol0.iterations)
end

@testset verbose = true showtiming = true ":non_uniform_grid" begin
    time_grid = [0,0.1,0.3,0.6,0.98,0.99,1]
    sol6 = solve(ocp, time_grid=time_grid, print_level=0)
    @test sol6.objective ≈ 0.309 rtol=1e-2
end



# recheck solution (T=2) with explicit / non-uniform grid
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

@testset verbose = true showtiming = true ":explicit_grid" begin
    solT2_exp = solve(ocpT2, time_grid=LinRange(0,1,CTDirect.__grid_size_direct()+1),print_level=0)
    @test (solT2_exp.objective == solT2.objective) && (solT2_exp.iterations == solT2.iterations)
end

@testset verbose = true showtiming = true ":non_uniform_grid" begin
    solT2_nonunif = solve(ocpT2, time_grid=[0,0.3,1,1.9,2],print_level=0)
    @test solT2_nonunif.objective ≈ 2.43 rtol=1e-2
end
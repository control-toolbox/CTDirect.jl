using CTDirect

println("Test: simple integrator")

# min energy
ocp1 = Model()
state!(ocp1, 1)
control!(ocp1, 1)
time!(ocp1, [0, 1])
constraint!(ocp1, :initial, -1, :initial_constraint)
constraint!(ocp1, :final, 0, :final_constraint)
dynamics!(ocp1, (x, u) -> -x + u)
objective!(ocp1, :lagrange, (x, u) -> u^2)

# with default initial guess (ie 0.1)
@testset verbose = true showtiming = true ":simple_integrator :min_tf" begin
    sol1 = solve(ocp1, grid_size=100, print_level=0, tol=1e-12)
    @test sol1.objective ≈ 0.313 rtol=1e-2
end

# with initial guess (using both vector and scalar syntax, no optimization variables)
@testset verbose = true showtiming = true ":simple_integrator :min_tf :init_constant" begin
    init_constant = OptimalControlInit(x_init=[-0.5], u_init=0)
    sol2 = solve(ocp1, grid_size=100, print_level=0, tol=1e-12, init=init_constant)
    @test sol2.objective ≈ 0.313 rtol=1e-2
end

# with initial guess from solution
sol = solve(ocp1, grid_size=100, print_level=0, tol=1e-12)
@testset verbose = true showtiming = true ":simple_integrator :min_tf :init_sol" begin
    init_sol = OptimalControlInit(sol)
    sol3 = solve(ocp1, grid_size=100, print_level=0, tol=1e-12, init=init_sol)
    @test sol3.objective ≈ 0.313 rtol=1e-2
end
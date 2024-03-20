using CTDirect
using CTBase

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
#sol1 = solve(ocp1, grid_size=100, print_level=0, tol=1e-12)
docp1 = DirectTranscription(ocp1, grid_size=100)
sol1 = solveDOCP(docp1, print_level=0, tol=1e-12)
@testset verbose = true showtiming = true ":double_integrator :min_tf" begin
    @test sol1.objective ≈ 0.313 rtol=1e-2
end

# with initial guess (using both vector and scalar syntax, no optimization variables)
init_constant = OptimalControlInit(x_init=[-0.5], u_init=0)
#sol2 = solve(ocp1, grid_size=100, print_level=0, tol=1e-12, init=init_constant)
docp2 = DirectTranscription(ocp1, grid_size=100, init=init_constant)
sol2 = solveDOCP(docp2, print_level=0, tol=1e-12)
@testset verbose = true showtiming = true ":double_integrator :min_tf :init_constant" begin
    @test sol2.objective ≈ 0.313 rtol=1e-2
end

# with initial guess from solution
init_sol = OptimalControlInit(sol2)
println(init_sol)
#sol3 = solve(ocp1, grid_size=100, print_level=0, tol=1e-12, init=init_sol)
# tis one fails due to state dimension mismatch somehow ?
docp3 = DirectTranscription(ocp1, grid_size=100, init=init_sol)
sol3 = solveDOCP(docp3, print_level=0, tol=1e-12)
@testset verbose = true showtiming = true ":double_integrator :min_tf :init_sol" begin
    @test sol3.objective ≈ 0.313 rtol=1e-2
end
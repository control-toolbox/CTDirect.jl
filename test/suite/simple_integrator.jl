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
init_constant = [-0.5, 0]
sol1 = solve(ocp1, grid_size=100, print_level=0, tol=1e-12, init=init_constant)
@testset verbose = true showtiming = true ":double_integrator :min_tf" begin
    @test sol1.objective ≈ 0.313 rtol=1e-2
end

println("Test: constraint types")

# box constraints
@testset verbose = true showtiming = true ":goddard :box_constraints" begin
    ocp = goddard()
    sol = solve(ocp.ocp, print_level=0, tol=1e-8)
    @test sol.objective ≈ 1.0125 rtol=1e-2
end

# functional constraints
@testset verbose = true showtiming = true ":goddard :functional_constraints" begin
    ocp1 = goddard(functional_constraints=true)
    sol1 = solve(ocp1.ocp, print_level=0, tol=1e-8, init=ocp1.init)
    @test sol1.objective ≈ 1.0125 rtol=1e-2
end
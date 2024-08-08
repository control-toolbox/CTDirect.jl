println("Test: constraint types")

# box constraints
@testset verbose = true showtiming = true ":goddard :box_constraints" begin
    ocp = goddard()
    sol = direct_solve(ocp.ocp, display=false)
    @test sol.objective ≈ 1.0125 rtol=1e-2
end

# functional constraints
@testset verbose = true showtiming = true ":goddard :functional_constraints" begin
    ocp1 = goddard(functional_constraints=true)
    sol1 = direct_solve(ocp1.ocp, display=false, init=ocp1.init)
    @test sol1.objective ≈ 1.0125 rtol=1e-2
end
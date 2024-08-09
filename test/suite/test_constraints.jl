println("Test: constraint types")

if !isdefined(Main, :goddard)
    include("../problems/goddard.jl")
end

# box constraints
@testset verbose = true showtiming = true ":goddard :box_constraints" begin
    ocp = goddard()
    sol = direct_solve(ocp.ocp, display=false)
    @test sol.objective ≈ ocp.obj rtol=1e-2
end

# functional constraints
@testset verbose = true showtiming = true ":goddard :functional_constraints" begin
    ocp = goddard(functional_constraints=true)
    sol = direct_solve(ocp.ocp, display=false, init=ocp.init)
    @test sol.objective ≈ ocp.obj rtol=1e-2
end
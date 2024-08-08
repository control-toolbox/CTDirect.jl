println("Test: abstract OCP definition")

# double integrator min tf, abstract definition
@testset verbose = true showtiming = true ":double_integrator :min_tf :abstract" begin
    prob = double_integrator_a()
    sol = direct_solve(prob.ocp, display=false)    
    @test sol.objective ≈ prob.obj rtol=1e-2
end

@testset verbose = true showtiming = true ":goddard :max_rf :abstract :constr" begin
    prob = goddard_a()
    sol = direct_solve(prob.ocp, display=false)    
    @test sol.objective ≈ prob.obj rtol=1e-2
end

println("Test: OCP definition")

# + beam

# double integrator min tf
if !isdefined(Main, :double_integrator_mintf)
    include("../problems/double_integrator.jl")
end

@testset verbose = true showtiming = true ":double_integrator :min_tf" begin
    prob = double_integrator_mintf()
    sol = direct_solve(prob.ocp, display = false)
    @test sol.objective ≈ prob.obj rtol = 1e-2
end

# + fuller

if !isdefined(Main, :goddard)
    include("../problems/goddard.jl")
end
@testset verbose = true showtiming = true ":goddard :max_rf" begin
    prob = goddard()
    sol = direct_solve(prob.ocp, display = false)
    @test sol.objective ≈ prob.obj rtol = 1e-2
end

# + jackson

# + robbins

# + simple integrator

# + vanderpol
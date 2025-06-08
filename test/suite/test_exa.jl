println("testing: OCP definitions (:exa)")


if CUDA.functional()
    exa_backend = CUDA.functional() ? CUDABackend() : nothing
else
    println("********** CUDA not available")
    exa_backend = nothing 
end

# beam
if !isdefined(Main, :beam2)
    include("../problems/beam2.jl")
end
@testset verbose = true showtiming = true ":beam2" begin
    prob = beam2()
    sol = solve(prob.ocp, :madnlp; nlp_model = :exa, disc_method = :euler, exa_backend = exa_backend) #, display=false)
    @test sol.objective ≈ prob.obj rtol = 1e-2
end

@ignore begin # debug
    
# double integrator min tf
if !isdefined(Main, :double_integrator_mintf)
    include("../problems/double_integrator.jl")
end
@testset verbose = true showtiming = true ":double_integrator :min_tf" begin
    prob = double_integrator_mintf()
    sol = solve(prob.ocp, display=false)
    @test sol.objective ≈ prob.obj rtol = 1e-2
end

# fuller
if !isdefined(Main, :fuller)
    include("../problems/fuller.jl")
end
@testset verbose = true showtiming = true ":fuller" begin
    prob = fuller()
    sol = solve(prob.ocp, display=false)
    @test sol.objective ≈ prob.obj rtol = 1e-2
end

# goddard max rf
if !isdefined(Main, :goddard)
    include("../problems/goddard.jl")
end
@testset verbose = true showtiming = true ":goddard :max_rf" begin
    prob = goddard()
    sol = solve(prob.ocp, display=false)
    @test sol.objective ≈ prob.obj rtol = 1e-2
end

# jackson
if !isdefined(Main, :jackson)
    include("../problems/jackson.jl")
end
@testset verbose = true showtiming = true ":jackson" begin
    prob = jackson()
    sol = solve(prob.ocp, display=false)
    @test sol.objective ≈ prob.obj rtol = 1e-2
end

# robbins
if !isdefined(Main, :robbins)
    include("../problems/robbins.jl")
end
@testset verbose = true showtiming = true ":robbins" begin
    prob = robbins()
    sol = solve(prob.ocp, display=false)
    @test sol.objective ≈ prob.obj rtol = 1e-2
end

# simple integrator
if !isdefined(Main, :simple_integrator)
    include("../problems/simple_integrator.jl")
end
@testset verbose = true showtiming = true ":simple_integrator" begin
    prob = simple_integrator()
    sol = solve(prob.ocp, display=false)
    @test sol.objective ≈ prob.obj rtol = 1e-2
end

# vanderpol
if !isdefined(Main, :vanderpol)
    include("../problems/vanderpol.jl")
end
@testset verbose = true showtiming = true ":vanderpol" begin
    prob = vanderpol()
    sol = solve(prob.ocp, display=false)
    @test sol.objective ≈ prob.obj rtol = 1e-2
end

end # debug
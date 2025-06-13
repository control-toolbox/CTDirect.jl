import ExaModels
using MadNLPGPU
using CUDA

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
# cpu case
@testset verbose = true showtiming = true ":examodel :cpu :beam2" begin
    prob = beam2()
    sol = solve(prob.ocp, :madnlp, :exa; disc_method = :euler, display=false)
    @test sol.objective ≈ prob.obj rtol = 1e-2
end
@testset verbose = true showtiming = true ":examodel :cpu :beam2 :init" begin
    prob = beam2()
    sol = solve(prob.ocp, :madnlp, :exa; disc_method = :euler, display=false, init=(control=0.0,), max_iter = 0)
    # NB. default init would be 0.1
    @test control(sol)(0) == 0e0
end

# gpu case
if !isnothing(exa_backend) 
@testset verbose = true showtiming = true ":examodel :GPU :beam2" begin
    prob = beam2()
    sol = solve(prob.ocp, :madnlp, :exa; disc_method = :euler, exa_backend = exa_backend, display=false)
    @test sol.objective ≈ prob.obj rtol = 1e-2
end
else
    println("skip GPU test")
end
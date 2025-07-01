# Test on GPU +++ split to test_gpu separate script

import ExaModels
using MadNLPGPU
using CUDA

# beam2
if !isdefined(Main, :beam2)
    include("../problems/beam2.jl")
end

println("testing: GPU with ExaModels / MadNLPGPU / CUDA")

CUDA.functional()

if CUDA.functional()
    println("********** CUDA is available")
    exa_backend = CUDA.functional() ? CUDABackend() : nothing
else
    println("********** CUDA is not available")
    exa_backend = nothing 
end
 
if !isnothing(exa_backend)
println("GPU tests with ", exa_backend)
@testset verbose = true showtiming = true ":examodel :GPU :euler" begin
    prob = beam2()
    sol = solve(prob.ocp, :madnlp, :exa; disc_method = :euler, exa_backend = exa_backend, display = display)
    @test sol.objective ≈ prob.obj rtol = 1e-2
end
@testset verbose = true showtiming = true ":examodel :GPU :trapeze :init" begin
    prob = beam2()
    sol = solve(prob.ocp, :madnlp, :exa; disc_method = :trapeze, exa_backend = exa_backend, display = display, init=(control=6.66,), max_iter = 0)
    @test control(sol)(0.5) == 6.66
end
else
    println("skipping GPU tests")
end

@ignore begin # debug
# swimmer2
if !isdefined(Main, :swimmer2)
    include("../problems/swimmer2.jl")
end

@testset verbose = true showtiming = true ":examodel :cpu :swimmer2" begin
    prob = swimmer2()
    sol = solve(prob.ocp, :madnlp, :exa; disc_method = :trapeze, display = display)
    @test sol.objective ≈ prob.obj rtol = 1e-2
end
end # debug
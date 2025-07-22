# ExaModels tests, CPU + GPU (when available)
import ExaModels
using MadNLP
using MadNLPGPU
using CUDA

# beam2

if !isdefined(Main, :beam2)
    include("./problems/beam2.jl")
end

# test_exa for all backends (CPU + GPU)

function test_exa(exa_backend)

    @testset verbose = true showtiming = true ":examodel :euler" begin
        prob = beam2()
        sol = solve(prob.ocp, :madnlp, :exa; disc_method = :euler, exa_backend = exa_backend, display = display)
        @test sol.objective ≈ prob.obj rtol = 1e-2
    end

    @testset verbose = true showtiming = true ":examodel :trapeze" begin
        prob = beam2()
        sol = solve(prob.ocp, :madnlp, :exa; disc_method = :trapeze, exa_backend = exa_backend, display = display)
        @test sol.objective ≈ prob.obj rtol = 1e-2
    end

    @testset verbose = true showtiming = true ":examodel :trapeze :grid_size" begin
        prob = beam2()
        sol = solve(prob.ocp, :madnlp, :exa; disc_method = :trapeze, exa_backend = exa_backend, display = display, grid_size = 1000)
        @test sol.objective ≈ prob.obj rtol = 1e-2
    end

    @testset verbose = true showtiming = true ":examodel :trapeze :init" begin
        prob = beam2()
        sol = solve(prob.ocp, :madnlp, :exa; disc_method = :trapeze, exa_backend = exa_backend, display = display, init=(control=6.66,), max_iter = 0)
        @test control(sol)(0.5) == 6.66
    end
    
end

# CPU tests

display = true
println("testing: CPU with ExaModels + MadNLP")
test_exa(nothing) # CPU tests

# GPU tests

display = true
if CUDA.functional()
    println("testing: GPU with ExaModels + MadNLPGPU + CUDA")
    test_exa(CUDABackend()) # GPU tests
else
    println("********** CUDA is not available")
end
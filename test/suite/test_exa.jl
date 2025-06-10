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
@testset verbose = true showtiming = true ":examodel :beam2" begin
    prob = beam2()
    sol = solve(prob.ocp, :madnlp; nlp_model = :exa, disc_method = :euler, exa_backend = exa_backend)
    #sol = solve(prob.ocp, :madnlp; nlp_model = :exa, disc_method = :euler) # debug: no exa_backend (copy from GPUarray issue...) #, display=false)
    @test sol.objective ≈ prob.obj rtol = 1e-2
end

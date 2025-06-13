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

# beam2
if !isdefined(Main, :beam2)
    include("../problems/beam2.jl")
end

@testset verbose = true showtiming = true ":examodel :cpu :beam2" begin
    prob = beam2()
    sol = solve(prob.ocp, :madnlp, :exa; disc_method = :trapeze, display = true)
    @test sol.objective ≈ prob.obj rtol = 1e-2
end

## @testset verbose = true showtiming = true ":examodel :cpu :beam2 :init" begin
##     prob = beam2()
##     sol = solve(prob.ocp, :madnlp, :exa; disc_method = :euler, display=false, init=(control=0.0,), max_iter = 0)
##     # NB. default init would be 0.1
##     @test control(sol)(0) == 0e0
## end
 
if !isnothing(exa_backend) 
@testset verbose = true showtiming = true ":examodel :GPU :beam2" begin
    prob = beam2()
    sol = solve(prob.ocp, :madnlp, :exa; disc_method = :trapeze, exa_backend = exa_backend, display = true)
    @test sol.objective ≈ prob.obj rtol = 1e-2
end
else
    println("skip GPU test")
end

@ignore begin # debug
# swimmer2
if !isdefined(Main, :swimmer2)
    include("../problems/swimmer2.jl")
end

@testset verbose = true showtiming = true ":examodel :cpu :swimmer2" begin
    prob = swimmer2()
    sol = solve(prob.ocp, :madnlp, :exa; disc_method = :trapeze, display = true)
    @test sol.objective ≈ prob.obj rtol = 1e-2
end
end # debug
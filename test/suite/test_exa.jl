# ExaModels tests, CPU + GPU (when available)
# +++ later move these to runtest with the others
using ExaModels: ExaModels
using MadNLPMumps
using MadNLPGPU
using CUDA
using AMDGPU

# load problems
if !isdefined(Main, :beam2)
    include("../problems/beam.jl")
end
if !isdefined(Main, :goddard2)
    include("../problems/goddard.jl")
end
if !isdefined(Main, :double_integrator_nobounds)
    include("../problems/double_integrator.jl")
end

# test_exa for all backends (CPU + GPU)
function test_exa(exa_backend, display)

    # beam2
    @testset verbose = true showtiming = true "beam2 :examodel :euler" begin
        prob = beam2()
        sol = solve(
            prob.ocp,
            :madnlp,
            :exa;
            disc_method=:euler,
            exa_backend=exa_backend,
            display=display,
        )
        @test objective(sol) ≈ prob.obj rtol = 1e-2
    end

    @testset verbose = true showtiming = true "beam2 :examodel :trapeze" begin
        prob = beam2()
        sol = solve(
            prob.ocp,
            :madnlp,
            :exa;
            disc_method=:trapeze,
            exa_backend=exa_backend,
            display=display,
        )
        @test objective(sol) ≈ prob.obj rtol = 1e-2
    end

    @testset verbose = true showtiming = true "beam2 :examodel :trapeze :grid_size" begin
        prob = beam2()
        sol = solve(
            prob.ocp,
            :madnlp,
            :exa;
            disc_method=:trapeze,
            exa_backend=exa_backend,
            display=display,
            grid_size=1000,
        )
        @test objective(sol) ≈ prob.obj rtol = 1e-2
    end

    @testset verbose = true showtiming = true "beam2 :examodel :trapeze :init" begin
        xi1 = 0.05
        xi2 = 2
        ui = 5
        prob = beam2()
        sol = solve(
            prob.ocp,
            :madnlp,
            :exa;
            disc_method=:trapeze,
            exa_backend=exa_backend,
            grid_size=4,
            display=display,
            init=(state=[xi1, xi2], control=ui),
            max_iter=0,
        )
        @test control(sol)(0.5) == ui
    end

    # no bounds
    @testset verbose = true showtiming = true "nobounds :examodel :madnlp :ipopt" begin
        prob = double_integ_nobounds()
        sol = solve(
            prob.ocp,
            :madnlp,
            :exa;
            disc_method=:trapeze,
            exa_backend=exa_backend,
            display=display,
        )
        @test objective(sol) ≈ prob.obj rtol = 1e-2
            sol = solve(
            prob.ocp,
            :ipopt,
            :exa;
            disc_method=:trapeze,
            exa_backend=exa_backend,
            display=display,
        )
        @test objective(sol) ≈ prob.obj rtol = 1e-2
    end

    # goddard2
    @testset verbose = true showtiming = true "goddard2 :examodel :trapeze :grid_size" begin
        prob = goddard2()
        sol = solve(
            prob.ocp,
            :madnlp,
            :exa;
            disc_method=:trapeze,
            exa_backend=exa_backend,
            display=display,
            grid_size=1000,
        )
        @test time_grid(sol)[end] ≈ 0.201965 rtol = 1e-2  # check time grid
        @test objective(sol) ≈ prob.obj rtol = 1e-2
    end

    @testset verbose = true showtiming = true ":examodel :cpu :transcription :grid_size" begin
        prob = beam2()
        docp = direct_transcription(
            prob.ocp, 
            :madnlp, 
            :exa; 
            display=display,
            disc_method=:trapeze, 
            grid_size=100)
        @test docp.dim_NLP_variables == 303
    end
end

# CPU tests
display = false
println("testing: ExaModels on CPU (MadNLP)")
test_exa(nothing, display) # CPU tests

# GPU tests (moonshot workflow)
display = false
if CUDA.functional()
    println("testing: ExaModels on GPU (MadNLPGPU / CUDA)")
    test_exa(CUDABackend(), display) # GPU tests
elseif AMDGPU.functional()
    println("testing: ExaModels on GPU (MadNLPGPU / AMDGPU)")
    test_exa(ROCBackend(), display) # GPU tests
else
    println("********** No GPU available, skipping GPU tests")
end

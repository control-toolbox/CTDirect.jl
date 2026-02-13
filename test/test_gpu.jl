# GPU only tests (do not run on CPU)
using MadNLPGPU
using CUDA
using AMDGPU

# load problems
if !isdefined(Main, :beam2)
    include("./problems/beam.jl")
end
if !isdefined(Main, :goddard2)
    include("./problems/goddard.jl")
end
if !isdefined(Main, :double_integrator_nobounds)
    include("./problems/double_integrator.jl")
end


function test_exa(exa_backend, display; linear_solver=CUDSSSolver)

    # beam2
    @testset verbose = true showtiming = true "beam2 :examodel :euler" begin
        test_problem(beam2(); solver=:madnlp, modeler=:exa, scheme=:euler,
                exa_backend=exa_backend, display=display, linear_solver=linear_solver)
    end

    @testset verbose = true showtiming = true "beam2 :examodel :trapeze" begin
        test_problem(beam2(); solver=:madnlp, modeler=:exa, scheme=:trapeze,
                exa_backend=exa_backend, display=display, linear_solver=linear_solver)
    end

    @testset verbose = true showtiming = true "beam2 :examodel :trapeze :grid_size" begin
        test_problem(beam2(); solver=:madnlp, modeler=:exa, scheme=:trapeze,
                exa_backend=exa_backend, display=display, grid_size=1000, linear_solver=linear_solver)
    end

    @testset verbose = true showtiming = true "beam2 :examodel :trapeze :init" begin
        init = (state=[0.05, 2], control=5)
        sol = solve_problem(beam2(); solver=:madnlp, modeler=:exa, scheme=:trapeze,
                exa_backend=exa_backend, display=display, init=init, max_iter=0, linear_solver=linear_solver)
        @test control(sol)(0.5) == 5
    end

    # no bounds
    @testset verbose = true showtiming = true "nobounds :examodel :madnlp" begin
        test_problem(double_integrator_nobounds(); solver=:madnlp, modeler=:exa, scheme=:trapeze,
                exa_backend=exa_backend, display=display, linear_solver=linear_solver)
    end

    # goddard2
    @testset verbose = true showtiming = true "goddard2 :examodel :trapeze :grid_size :objective" begin
        prob = goddard2()
        sol = solve_problem(prob; solver=:madnlp, modeler=:exa, scheme=:trapeze,
                exa_backend=exa_backend, display=display, linear_solver=linear_solver, tol=1e-8)
        @test time_grid(sol)[end] ≈ 0.2014 rtol = 1e-2  # check time grid
        @test objective(sol) ≈ prob.obj rtol = 1e-2
    end
end

# GPU tests (kkt workflow)
display = true
if CUDA.functional()
    println("testing: ExaModels on GPU (MadNLPGPU / CUDA)")
    test_exa(CUDABackend(), display) # GPU tests
elseif AMDGPU.functional()
    println("testing: ExaModels on GPU (MadNLPGPU / AMDGPU)")
    test_exa(ROCBackend(), display) # GPU tests
else
    println("********** No GPU available, skipping GPU tests")
end

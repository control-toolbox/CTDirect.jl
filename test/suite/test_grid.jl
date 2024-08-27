println("Test: grid options")

# 1. simple integrator min energy (dual control for test)
if !isdefined(Main, :simple_integrator)
    include("../problems/simple_integrator.jl")
end
ocp = simple_integrator().ocp
sol0 = direct_solve(ocp, display = false)

# solve with explicit and non uniform time grid
@testset verbose = true showtiming = true ":explicit_grid" begin
    time_grid = LinRange(0, 1, CTDirect.__grid_size() + 1)
    sol = direct_solve(ocp, time_grid = time_grid, display = false)
    @test (sol.objective == sol0.objective) && (sol.iterations == sol0.iterations)
end

@testset verbose = true showtiming = true ":non_uniform_grid" begin
    time_grid = [0, 0.1, 0.3, 0.6, 0.98, 0.99, 1]
    sol = direct_solve(ocp, time_grid = time_grid, display = false)
    @test sol.objective ≈ 0.309 rtol = 1e-2
end

# 2. integrator free times
if !isdefined(Main, :double_integrator_freet0tf)
    include("../problems/double_integrator.jl")
end
ocp = double_integrator_freet0tf().ocp
sol0 = direct_solve(ocp, display = false)

@testset verbose = true showtiming = true ":explicit_grid" begin
    sol = direct_solve(ocp, time_grid = LinRange(0, 1, CTDirect.__grid_size() + 1), display = false)
    @test (sol.objective == sol0.objective) && (sol.iterations == sol0.iterations)
end

@testset verbose = true showtiming = true ":max_t0 :non_uniform_grid" begin
    sol = direct_solve(ocp, time_grid = [0, 0.1, 0.6, 0.95, 1], display = false)
    @test sol.objective ≈ 7.48 rtol = 1e-2
end

# 3. parametric ocp (T=2) with explicit / non-uniform grid
if !isdefined(Main, :double_integrator_T)
    include("../problems/double_integrator.jl")
end
ocp = double_integrator_T(2).ocp
sol0 = direct_solve(ocp, display = false)

@testset verbose = true showtiming = true ":explicit_grid" begin
    sol = direct_solve(ocp, time_grid = LinRange(0, 1, CTDirect.__grid_size() + 1), display = false)
    @test (sol.objective == sol0.objective) && (sol.iterations == sol0.iterations)
end

@testset verbose = true showtiming = true ":non_uniform_grid" begin
    sol = direct_solve(ocp, time_grid = [0, 0.3, 1, 1.9, 2], display = false)
    @test sol.objective ≈ 2.43 rtol = 1e-2
end

println("Test: discretization options")

normalize_grid(t) = return (t .- t[1]) ./ (t[end] - t[1])

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
    @test sol.time_grid ≈ time_grid
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
    time_grid = [0, 0.1, 0.6, 0.95, 1]
    sol = direct_solve(ocp, time_grid = time_grid, display = false)
    @test normalize_grid(sol.time_grid) ≈ time_grid
end

# 3. parametric ocp (T=2) with explicit / non-uniform grid
if !isdefined(Main, :double_integrator_minenergy)
    include("../problems/double_integrator.jl")
end
ocp = double_integrator_minenergy(2).ocp
sol0 = direct_solve(ocp, display = false)

@testset verbose = true showtiming = true ":explicit_grid" begin
    sol = direct_solve(ocp, time_grid = LinRange(0, 1, CTDirect.__grid_size() + 1), display = false)
    @test (sol.objective == sol0.objective) && (sol.iterations == sol0.iterations)
end

@testset verbose = true showtiming = true ":non_uniform_grid" begin
    time_grid = [0, 0.3, 1, 1.9, 2]
    sol = direct_solve(ocp, time_grid = time_grid, display = false)
    @test sol.time_grid ≈ time_grid
end

# implicit midpoint scheme
@testset verbose = true showtiming = true ":implicit_midpoint" begin
    ocp = simple_integrator().ocp
    sol_t = direct_solve(ocp, display = false)
    sol_m = direct_solve(ocp, display = false, disc_method = "midpoint")
    @test sol_m.objective ≈ sol_t.objective rtol = 1e-2 
end

@testset verbose = true showtiming = true ":implicit_midpoint" begin
    ocp = double_integrator_freet0tf().ocp
    sol_t = direct_solve(ocp, display = false)
    sol_m = direct_solve(ocp, display = false, disc_method = :midpoint)
    @test sol_m.objective ≈ sol_t.objective rtol = 1e-2 
end

println("testing: discretization options")

normalize_grid(t) = return (t .- t[1]) ./ (t[end] - t[1])
N = 250

# 1. simple integrator min energy (dual control for test)
if !isdefined(Main, :simple_integrator)
    include("../problems/simple_integrator.jl")
end
ocp = simple_integrator().ocp
sol0 = solve(ocp; display=false)

# solve with explicit and non uniform time grid
@testset verbose = true showtiming = true ":explicit_grid" begin
    sol = solve(ocp, time_grid=LinRange(0, 1, N + 1), display=false)
    @test (objective(sol) == objective(sol0)) && (iterations(sol) == iterations(sol0))
end

@testset verbose = true showtiming = true ":non_uniform_grid" begin
    grid = [0, 0.1, 0.3, 0.6, 0.98, 0.99, 1]
    sol = solve(ocp, time_grid=grid, display=false)
    @test time_grid(sol) ≈ grid
end

# 2. integrator free times
if !isdefined(Main, :double_integrator_freet0tf)
    include("../problems/double_integrator.jl")
end
ocp = double_integrator_freet0tf().ocp
sol0 = solve(ocp; display=false)

@testset verbose = true showtiming = true ":explicit_grid" begin
    sol = solve(ocp, time_grid=LinRange(0, 1, N + 1), display=false)
    @test (objective(sol) == objective(sol0)) && (iterations(sol) == iterations(sol0))
end

@testset verbose = true showtiming = true ":max_t0 :non_uniform_grid" begin
    grid = [0, 0.1, 0.6, 0.95, 1]
    sol = solve(ocp, time_grid=grid, display=false)
    @test normalize_grid(time_grid(sol)) ≈ grid
end

# 3. double integrator min energy ocp (T=2) with explicit / non-uniform grid
if !isdefined(Main, :double_integrator_minenergy)
    include("../problems/double_integrator.jl")
end
ocp = double_integrator_minenergy(2).ocp
sol0 = solve(ocp; display=false)

@testset verbose = true showtiming = true ":explicit_grid" begin
    sol = solve(ocp, time_grid=LinRange(0, 1, N + 1), display=false)
    @test (objective(sol) == objective(sol0)) && (iterations(sol) == iterations(sol0))
end

@testset verbose = true showtiming = true ":non_uniform_grid" begin
    grid = [0, 0.3, 1, 1.9, 2]
    sol = solve(ocp, time_grid=grid, display=false)
    @test time_grid(sol) ≈ grid
end

# discretization methods: lagrange with constraint
@testset verbose = true showtiming = true ":beam :disc_method" begin
    for disc_method in
        [:trapeze, :midpoint, :euler, :euler_implicit, :gauss_legendre_2, :gauss_legendre_3]
        check_problem(beam(), display=false, disc_method=disc_method)
    end
end

# discretization methods: mayer with free t0 and tf
@testset verbose = true showtiming = true ":double_integrator :disc_method" begin
    for disc_method in
        [:trapeze, :midpoint, :euler, :euler_implicit, :gauss_legendre_2, :gauss_legendre_3]
        check_problem(
            double_integrator_freet0tf(),
            display=false,
            disc_method=disc_method,
            grid_size=50,
        )
    end
end

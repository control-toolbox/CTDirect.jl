println("testing: discretization options")

normalize_grid(t) = return (t .- t[1]) ./ (t[end] - t[1])
N = 250

# 1. simple integrator min energy (dual control for test) with explicit / non-uniform grid
if !isdefined(Main, :simple_integrator)
    include("../problems/simple_integrator.jl")
end
prob = simple_integrator()
sol0 = solve_problem(prob; display=false)

@testset verbose = true showtiming = true ":explicit_grid" begin
    sol = solve_problem(prob; time_grid=LinRange(0, 1, N + 1))
    @test (objective(sol) == objective(sol0)) && (iterations(sol) == iterations(sol0))
end

@testset verbose = true showtiming = true ":non_uniform_grid" begin
    grid = [0, 0.1, 0.3, 0.6, 0.98, 0.99, 1]
    sol = solve_problem(prob; time_grid=grid)
    @test time_grid(sol) ≈ grid
end

# 2. integrator free times with explicit / non-uniform grid
if !isdefined(Main, :double_integrator_freet0tf)
    include("../problems/double_integrator.jl")
end
prob = double_integrator_freet0tf()
sol0 = solve_problem(prob; display=false)

@testset verbose = true showtiming = true ":explicit_grid" begin
    sol = solve_problem(prob; time_grid=LinRange(0, 1, N + 1))
    @test (objective(sol) == objective(sol0)) && (iterations(sol) == iterations(sol0))
end

@testset verbose = true showtiming = true ":max_t0 :non_uniform_grid" begin
    grid = [0, 0.1, 0.6, 0.95, 1]
    sol = solve_problem(prob; time_grid=grid)
    @test normalize_grid(time_grid(sol)) ≈ grid
end

# 3. double integrator min energy (T=2) with explicit / non-uniform grid
if !isdefined(Main, :double_integrator_minenergy)
    include("../problems/double_integrator.jl")
end
prob = double_integrator_minenergy(2)
sol0 = solve_problem(prob; display=false)

@testset verbose = true showtiming = true ":explicit_grid" begin
    sol = solve_problem(prob; time_grid=LinRange(0, 1, N + 1))
    @test (objective(sol) == objective(sol0)) && (iterations(sol) == iterations(sol0))
end

@testset verbose = true showtiming = true ":non_uniform_grid" begin
    grid = [0, 0.3, 1, 1.9, 2]
    sol = solve_problem(prob; time_grid=grid)
    @test time_grid(sol) ≈ grid
end


# discretization methods: 
scheme_list = [:trapeze, :midpoint, :euler, :euler_implicit, 
            :gauss_legendre_2, :gauss_legendre_3]

# lagrange with constraint
if !isdefined(Main, :beam)
    include("../problems/beam.jl")
end
@testset verbose = true showtiming = true ":beam :scheme" begin
    for scheme in scheme_list
        test_problem(beam(); scheme=scheme)
    end
end

# mayer with free t0 and tf
@testset verbose = true showtiming = true ":double_integrator :scheme" begin
    for scheme in scheme_list
        test_problem(double_integrator_freet0tf(); scheme=scheme, grid_size=50)
    end
end

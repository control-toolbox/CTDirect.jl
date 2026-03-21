println("testing: zero control dimension - parameter estimation")

# Load test problems
if !isdefined(Main, :estimate_initial_condition)
    include("../problems/autonomous_system.jl")
end

# Test 1: Estimate initial condition with all schemes
@testset verbose = true showtiming = true ":param_estimation :all_schemes" begin
    prob = estimate_initial_condition()
    
    @testset ":euler" begin
        sol = solve_problem(prob; scheme=:euler, grid_size=50, display=false)
        @test CTModels.successful(sol)
        @test sol.objective >= 0
    end
    
    @testset ":midpoint" begin
        sol = solve_problem(prob; scheme=:midpoint, grid_size=50, display=false)
        @test CTModels.successful(sol)
        @test sol.objective >= 0
    end
    
    @testset ":trapeze" begin
        sol = solve_problem(prob; scheme=:trapeze, grid_size=50, display=false)
        @test CTModels.successful(sol)
        @test sol.objective >= 0
    end
end

# Test 2: Estimate parameter in dynamics
@testset verbose = true showtiming = true ":param_estimation :rotation_rate" begin
    prob = estimate_rotation_rate()
    sol = solve_problem(prob; scheme=:midpoint, grid_size=50, display=false)
    @test CTModels.successful(sol)
    @test sol.objective >= 0
    # Verify variable is properly retrieved
    @test length(variable(sol)) == 1
end

# Test 3: Least squares with path constraint
@testset verbose = true showtiming = true ":param_estimation :with_constraint" begin
    prob = least_squares_with_constraint()
    sol = solve_problem(prob; scheme=:midpoint, grid_size=50, display=false)
    @test CTModels.successful(sol)
    @test sol.objective >= 0
end

# Test 4: Verify solution dimensions (zero control)
@testset verbose = true showtiming = true ":param_estimation :solution_dimensions" begin
    prob = estimate_initial_condition()
    sol = solve_problem(prob; scheme=:midpoint, grid_size=50, display=false)
    
    T = time_grid(sol, :state)
    
    # State must have 2 dimensions
    @test length(state(sol)(0.5)) == 2
    
    # Control must be empty
    @test length(control(sol)(0.5)) == 0
    
    # Verify control() doesn't crash
    @test control(sol) isa Function
    for t in T
        u = control(sol)(t)
        @test u isa AbstractVector
        @test length(u) == 0
    end
end

# Test 5: Initial guess with empty control
@testset verbose = true showtiming = true ":param_estimation :initial_guess" begin
    prob = estimate_initial_condition()
    
    # Functional initial guess
    x_init = t -> [cos(π*t/4), sin(π*t/4)]
    u_init = t -> Float64[]  # empty vector
    
    sol = solve_problem(prob; 
        scheme=:midpoint, 
        grid_size=50, 
        init=(state=x_init, control=u_init),
        display=false
    )
    @test CTModels.successful(sol)
end

# Test 6: ADNLP manual backend (sparsity patterns)
@testset verbose = true showtiming = true ":param_estimation :adnlp_manual" begin
    prob = estimate_initial_condition()
    sol = solve_problem(prob; 
        scheme=:midpoint, 
        grid_size=20, 
        adnlp_backend=:manual,
        display=false
    )
    @test CTModels.successful(sol)
end

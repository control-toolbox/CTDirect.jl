println("testing: zero control allocations")

using CTDirect
using CTModels

# Load test problems
if !isdefined(Main, :estimate_initial_condition)
    include("../problems/autonomous_system.jl")
end

# Create DOCP with zero control dimension
@testset verbose = true showtiming = true ":zero_control_allocations" begin
    
    @testset "DOCP creation" begin
        prob = estimate_initial_condition()
        ocp = prob.ocp
        
        # Build DOCP directly like in get_docp
        grid_size = 10
        control_steps = 1
        scheme = :midpoint
        time_grid = nothing
        docp = CTDirect.DOCP(ocp, grid_size, control_steps, scheme, time_grid)
        
        # Verify NLP_u = 0
        @test docp.dims.NLP_u == 0
        
        # Verify dim_NLP_variables is consistent
        # For midpoint: steps * (NLP_x + NLP_u) + NLP_x + NLP_v
        # = 10 * (2 + 0) + 2 + 2 = 24 (2 variables for initial condition)
        @test docp.dim_NLP_variables == 10 * 2 + 2 + 2
        
        # Verify dim_NLP_constraints is consistent
        @test docp.dim_NLP_constraints > 0
    end
    
    @testset "Control getters return empty views" begin
        prob = estimate_initial_condition()
        ocp = prob.ocp
        
        # Build DOCP directly
        grid_size = 10
        control_steps = 1
        scheme = :midpoint
        time_grid = nothing
        docp = CTDirect.DOCP(ocp, grid_size, control_steps, scheme, time_grid)
        
        # Create dummy xu vector
        xu = zeros(docp.dim_NLP_variables)
        
        # get_OCP_control_at_time_step should return empty view
        u = CTDirect.get_OCP_control_at_time_step(xu, docp, 1)
        @test u isa AbstractVector
        @test length(u) == 0
        @test eltype(u) == Float64
    end
    
    @testset "Control setters are no-ops" begin
        prob = estimate_initial_condition()
        ocp = prob.ocp
        
        # Build DOCP directly
        grid_size = 10
        control_steps = 1
        scheme = :midpoint
        time_grid = nothing
        docp = CTDirect.DOCP(ocp, grid_size, control_steps, scheme, time_grid)
        
        xu = zeros(docp.dim_NLP_variables)
        xu_copy = copy(xu)
        
        # set_control_at_time_step! should not modify anything
        CTDirect.set_control_at_time_step!(xu, Float64[], docp, 1)
        @test xu == xu_copy
    end
    
    @testset "All schemes support zero control" begin
        prob = estimate_initial_condition()
        ocp = prob.ocp
        
        for scheme in [:euler, :midpoint, :trapeze]
            docp = CTDirect.DOCP(ocp, 10, 1, scheme, nothing)
            @test docp.dims.NLP_u == 0
            @test docp.dim_NLP_variables > 0
        end
    end
    
    @testset "Variables bounds with zero control" begin
        prob = estimate_initial_condition()
        ocp = prob.ocp
        
        docp = CTDirect.DOCP(ocp, 10, 1, :midpoint, nothing)
        CTDirect.__variables_bounds!(docp)
        
        # Verify bounds exist and have correct dimensions
        @test length(docp.bounds.var_l) == docp.dim_NLP_variables
        @test length(docp.bounds.var_u) == docp.dim_NLP_variables
        # Bounds can contain -Inf/Inf for unbounded variables
        @test all(docp.bounds.var_l .<= docp.bounds.var_u)
    end
    
    @testset "Initial guess with zero control" begin
        prob = estimate_initial_condition()
        ocp = prob.ocp
        
        docp = CTDirect.DOCP(ocp, 10, 1, :midpoint, nothing)
        init = CTModels.build_initial_guess(ocp, ())
        
        x0 = CTDirect.__initial_guess(docp, init)
        
        @test length(x0) == docp.dim_NLP_variables
        @test all(isfinite, x0)
    end
    
    @testset "Control getters at all time steps" begin
        prob = estimate_initial_condition()
        ocp = prob.ocp
        
        docp = CTDirect.DOCP(ocp, 10, 1, :midpoint, nothing)
        xu = zeros(docp.dim_NLP_variables)
        
        # Test at all time steps
        for i in 1:(docp.time.steps + 1)
            u = CTDirect.get_OCP_control_at_time_step(xu, docp, i)
            @test length(u) == 0
            @test eltype(u) == Float64
        end
    end
    
    @testset "Zero control with optimization variable" begin
        prob = estimate_rotation_rate()
        ocp = prob.ocp
        
        docp = CTDirect.DOCP(ocp, 10, 1, :midpoint, nothing)
        
        @test docp.dims.NLP_u == 0
        @test docp.dims.NLP_v == 1  # one optimization variable
        @test docp.dim_NLP_variables == 10 * 2 + 2 + 1  # +1 for variable
    end
    
    @testset "Sparsity patterns with zero control" begin
        prob = estimate_initial_condition()
        ocp = prob.ocp
        
        docp = CTDirect.DOCP(ocp, 10, 1, :midpoint, nothing)
        
        # Jacobian pattern should be constructible
        Is, Js = CTDirect.DOCP_Jacobian_pattern(docp)
        @test length(Is) == length(Js)
        @test all(1 .<= Is .<= docp.dim_NLP_constraints)
        @test all(1 .<= Js .<= docp.dim_NLP_variables)
        
        # Hessian pattern should be constructible
        Is_h, Js_h = CTDirect.DOCP_Hessian_pattern(docp)
        @test length(Is_h) == length(Js_h)
        @test all(1 .<= Is_h .<= docp.dim_NLP_variables)
        @test all(1 .<= Js_h .<= docp.dim_NLP_variables)
    end
end

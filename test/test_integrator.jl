# double integrator - energy min
println("Double integrator test")
prob = Problem(:integrator, :dim2, :energy) 
ocp = prob.model


# solve
println("Is solvable ? ", CTDirect.is_solvable(ocp))




@testset verbose = true showtiming = true ":integrator :dim2 :energy" begin
    N = 1
    ctd = CTDirect.CTDirect_data(ocp, N, nothing)
    @testset verbose = true showtiming = true "CTDirect_data" begin
        @test ctd.dim_control_constraints == 0
        @test ctd.dim_state_constraints == 0
        @test ctd.dim_mixed_constraints == 0
        @test ctd.dim_boundary_conditions == 4
        @test ctd.dim_control_box == 0
        @test ctd.dim_state_box == 0
        @test ctd.has_control_constraints == false
        @test ctd.has_state_constraints   == false
        @test ctd.has_mixed_constraints   == false
        @test ctd.has_boundary_conditions == true
        @test ctd.has_control_box         == false
        @test ctd.has_state_box           == false
        @test ctd.has_lagrange_cost         == true
        @test ctd.has_mayer_cost            == false
        @test ctd.dim_NLP_state           == ocp.state_dimension + 1
        @test ctd.mayer                   == nothing
        f_Mayer_test(t,x,u)=[ocp.dynamics(t,x,u);ocp.lagrange(t,x,u)]
        @test ctd.dynamics_lagrange_to_mayer(0,[0;2],1) == f_Mayer_test(0,[0;2],1)
        @test ctd.has_free_final_time     == false
        @test ctd.dim_NLP_variables == (N+1)*4
        println("ctd.dim_NLP_constraints = ", ctd.dim_NLP_constraints)
        @test ctd.dim_NLP_constraints == 3*N+5
    end
    @testset verbose = true showtiming = true "constraints" begin
        lb, ub = CTDirect.constraints_bounds(ctd)
        l_var, u_var = CTDirect.variables_bounds(ctd)
        true_lb = zeros(3*N) # test without boundary conditions because of the dictionary
        # true_lb[end-4:end] = [-1,0,0,0,0]
        @test lb[1:3*N] == true_lb
        @test ub[1:3*N] == true_lb
        @test l_var == -Inf*ones((N+1)*4)
        @test u_var == -l_var
    end

    # Tests on the numerical solution
    @testset verbose = true showtiming = true "numerical solution" begin
        #using LinearAlgebra
        sol = solve(ocp, grid_size=100, print_level=0) 
        # check solution
        u_sol(t) = prob.solution.control(t)[1]
        u = t -> sol.control(t)[1]
        x_sol(t) = prob.solution.state(t)
        x = t -> sol.state(t)[1:2]
        T = sol.times
        dT = T[2:end]-T[1:end-1]
        N = length(T)
        # test on the infty norm of the state
        @test maximum([ norm(x(T[i])-x_sol(T[i])) for i in 1:N] ) ≈ 0 atol = 0.01
        #@test maximum([ sqrt(sum((x(T[i])-x_sol(T[i])).^2)) for i in 1:N]) ≈ 0 atol = 0.01
        # test on the L_2 norm of the control
        @test sum(dT .* [ abs(u(T[i])-u_sol(T[i])) for i ∈ 1:N-1] ) ≈ 0 atol=1e-1
        # test on the objectif
        @test sol.objective ≈ prob.solution.objective atol=1e-2
        # @test constraints_violation(sol) < 1e-6 # ceci n'existe pas dans la OptimalControlSolution pour le moment
    end
    @testset verbose = true showtiming = true "box_constraints" begin
      umax = 5.
      constraint!(ocp, :control, -umax, umax, :control_con1)
      ctd = CTDirect.CTDirect_data(ocp, 3, nothing)
      @test ctd.dim_control_constraints == 0
      @test ctd.has_control_box         == true
      @test ctd.has_control_constraints == false
      @test ctd.dim_control_box == 1
      @test ctd.dim_NLP_state == 3
   end

   @testset verbose = true showtiming = true "contol_constraints" begin
      umax = 5
      remove_constraint!(ocp, :control_con1)
      constraint!(ocp, :control, u -> u, -umax, umax, :control_con2)
      N  = 3
      ctd = CTDirect.CTDirect_data(ocp, N, nothing)
      @test ctd.dim_control_constraints == 1
      @test ctd.has_control_box         == false
      @test ctd.has_control_constraints == true
      @test ctd.dim_control_box == 0
      @test ctd.dim_control_constraints == 1
   end

   @testset verbose = true showtiming = true "mixed_constraints" begin
      umax = 5
      remove_constraint!(ocp, :control_con2)
      constraint!(ocp, :mixed, (x,u) -> u, -umax, umax, :control_con3)
      N  = 3
      ctd = CTDirect.CTDirect_data(ocp, N, nothing)
      @test ctd.dim_control_constraints == 0
      @test ctd.has_control_box         == false
      @test ctd.has_control_constraints == false
      @test ctd.has_mixed_constraints   == true
      @test ctd.dim_control_constraints == 0
      @test ctd.dim_mixed_constraints   == 1
   end
end



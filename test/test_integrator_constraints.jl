# double integrator - energy min
println("Double integrator test")
prob = Problem(:integrator, :dim2, :energy) 
ocp = prob.model

# solve
println("Is solvable ? ", CTDirect.is_solvable(ocp))

@testset verbose = true showtiming = true ":integrator :dim2 :energy with constraints" begin
    
    @testset verbose = true showtiming = true "control_box_constraints" begin
      umax = 5.
      N = 2
      constraint!(ocp, :control, -umax, umax, :control_con1)
      ctd = CTDirect.CTDirect_data(ocp, N, nothing)
      @test ctd.dim_control_constraints == 0
      @test ctd.has_control_box         == true
      @test ctd.has_control_constraints == false
      @test ctd.dim_control_box == 1
      @test ctd.dim_NLP_state == 3

      lb, ub = CTDirect.constraints_bounds(ctd)
      l_var, u_var = CTDirect.variables_bounds(ctd)
      true_lb = zeros(3*N) # test without boundary conditions because of the dictionary
      # true_lb[end-4:end] = [-1,0,0,0,0]
      true_l_var =  -Inf*ones((N+1)*4)
      true_l_var[3*(N+1)+1:4*(N+1)] .= -umax
      @test lb[1:end-5] == true_lb
      @test ub[1:end-5] == -true_lb
      @test l_var == true_l_var
      @test u_var == -l_var
   end

   @testset verbose = true showtiming = true "control_constraints" begin
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

      lb, ub = CTDirect.constraints_bounds(ctd)
      l_var, u_var = CTDirect.variables_bounds(ctd)
      true_lb = zeros(4*N+1) # test without boundary conditions because of the dictionary
      true_lb[4:4:end] .= -umax
      true_lb[end] = -umax
      # true_lb[end-4:end] = [-1,0,0,0,0]
      true_l_var =  -Inf*ones((N+1)*4)
      @test lb[1:end-5] == true_lb
      @test ub[1:end-5] == -true_lb
      @test l_var == true_l_var
      @test u_var == -l_var
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

      lb, ub = CTDirect.constraints_bounds(ctd)
      l_var, u_var = CTDirect.variables_bounds(ctd)
      true_lb = zeros(4*N+1) # test without boundary conditions because of the dictionary
      true_lb[4:4:end] .= -umax
      true_lb[end] = -umax
      # true_lb[end-4:end] = [-1,0,0,0,0]
      true_l_var =  -Inf*ones((N+1)*4)
      @test lb[1:end-5] == true_lb
      @test ub[1:end-5] == -true_lb
      @test l_var == true_l_var
      @test u_var == -l_var
   end

   @testset verbose = true showtiming = true "state_box_constraints" begin
    x_min1 = 1; x_min2 = -1.; x_max1 = 2; x_max2 = 5.
    N = 2
    remove_constraint!(ocp, :control_con3)
    constraint!(ocp, :state, [x_min1, x_min2], [x_max1, x_max2], :state_con1)
    ctd = CTDirect.CTDirect_data(ocp, N, nothing)
    @test ctd.dim_control_constraints == 0
    @test ctd.has_control_box         == false
    @test ctd.has_control_constraints == false
    @test ctd.dim_control_box == 0
    @test ctd.has_state_box == true


    lb, ub = CTDirect.constraints_bounds(ctd)
    l_var, u_var = CTDirect.variables_bounds(ctd)
    true_lb = zeros(3*N) # test without boundary conditions because of the dictionary
    # true_lb[end-4:end] = [-1,0,0,0,0]
    true_l_var =  -Inf*ones((N+1)*4)
    true_l_var[1:3:end-(N+1)] .= x_min1
    true_l_var[2:3:end-(N+1)] .= x_min2
    true_u_var =  Inf*ones((N+1)*4)
    true_u_var[1:3:end-(N+1)] .= x_max1
    true_u_var[2:3:end-(N+1)] .= x_max2
    @test lb[1:end-5] == true_lb
    @test ub[1:end-5] == -true_lb
    @test l_var == true_l_var
    @test u_var == true_u_var
 end
end



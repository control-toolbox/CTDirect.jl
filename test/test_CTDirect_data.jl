using CTDirect
using CTProblems
using CTBase # for the functions
using Test

println("CTDirect_data function tests")
# 
prob = Problem(:integrator, :dim2, :energy); 
ocp = prob.model

ctd = CTDirect.CTDirect_data(ocp, 3, nothing)

@testset verbose = true showtiming = true "CTDirect_data" begin

    @testset verbose = true showtiming = true "data" begin
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





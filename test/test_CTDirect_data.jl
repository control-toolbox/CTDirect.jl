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

   
#   @test control_constraints[2] == ocp.
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
#   @test criterion               == :min
end

constraint!(ocp, :control, -4.01, 4.01, :control_con1)

#constraint!(ocp, :control, -5, 5, :control_con1)

ctd = CTDirect.CTDirect_data(ocp, 3, nothing)

@testset verbose = true showtiming = true "box_constraints" begin
   @test ctd.dim_control_constraints == 0
   @test ctd.has_control_box         == true
   @test ctd.dim_control_box == 1
   @test ctd.dim_NLP_state == 3
end


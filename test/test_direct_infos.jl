using CTDirect
using CTProblems
using CTBase # for the functions
using Test

println("direct_infos function tests")
# 
prob = Problem(:integrator, :dim2, :energy); 
ocp = prob.model

N = 10
t0, tf, n_x, m, f, control_constraints, state_constraints, mixed_constraints, boundary_conditions, control_box, state_box, dim_control_constraints, dim_state_constraints, dim_mixed_constraints, dim_boundary_conditions, dim_control_box, dim_state_box,has_control_constraints, has_state_constraints, has_mixed_constraints, has_boundary_conditions, has_control_box, has_state_box, hasLagrangeCost, hasMayerCost, dim_x, nc, dim_xu, g, f_Mayer, has_free_final_time, criterion = CTDirect.direct_infos(ocp, N)

println("control_constraints = ", control_constraints)
println("state_constraints = ", state_constraints)
println("mixed_constraints = ", mixed_constraints)
println("boundary_conditions = ", boundary_conditions)
println("control_box = ", control_box)
println("state_box[1] = ", state_box[1])
println("state_box[2] = ", state_box[2])
println("state_box[3] = ", state_box[3])
println("control_constraints[1] = ", control_constraints[1])
println("control_constraints[2] = ", control_constraints[2])
println("control_constraints[3] = ", control_constraints[3])

println("state_box = ", state_box)

@testset verbose = true showtiming = true "direct_infos" begin
   @test t0 == ocp.initial_time
   @test tf == ocp.final_time
   @test n_x == ocp.state_dimension
   @test m == ocp.control_dimension
   @test f == ocp.dynamics

   @test dim_control_constraints == 0
   @test dim_state_constraints == 0
   @test dim_mixed_constraints == 0
   @test dim_boundary_conditions == 4
   @test dim_control_box == 0
   @test dim_state_box == 0
   @test has_control_constraints == false
   @test has_state_constraints   == false
   @test has_mixed_constraints   == false
   @test has_boundary_conditions == true
   @test has_control_box         == false
   @test has_state_box           == false
   @test hasLagrangeCost         == true
   @test hasMayerCost            == false
   @test dim_x                   == ocp.state_dimension + 1
   @test g                       == nothing
   f_Mayer_test(t,x,u)=[f(t,x,u);ocp.lagrange(t,x,u)]
   @test f_Mayer(0,[0;2],1)      == f_Mayer_test(0,[0;2],1)
   @test has_free_final_time     == false
   @test criterion               == :min
end

constraint!(ocp, :control, -4.01, 4.01, :control_con1)

#constraint!(ocp, :control, -5, 5, :control_con1)


t0, tf, n_x, m, f, control_constraints, state_constraints, mixed_constraints, boundary_conditions, control_box, state_box, dim_control_constraints, dim_state_constraints, dim_mixed_constraints, dim_boundary_conditions, dim_control_box, dim_state_box,has_control_constraints, has_state_constraints, has_mixed_constraints, has_boundary_conditions, has_control_box, has_state_box, hasLagrangeCost, hasMayerCost, dim_x, nc, dim_xu, g, f_Mayer, has_free_final_time, criterion = CTDirect.direct_infos(ocp, N)

@testset verbose = true showtiming = true "direct_infos_box_constraints" begin
   @test dim_control_constraints == 0
   @test has_control_box         == true
   @test dim_control_box == 1
   @test dim_x == 3
end


using CTDirect
using CTProblems
using CTBase # for the functions
using Test
#
println("constraints")
# 
prob = Problem(:integrator, :dim2, :energy); 
ocp = prob.model
umax = 5
constraint!(ocp, :control, -umax, umax, :control_con1)

@testset verbose = true showtiming = true "constraints" begin
    @testset verbose = true showtiming = true "box-control" begin
        N  = 3
        ctd = CTDirect.CTDirect_data(ocp, N, nothing)
        l_var, u_var = CTDirect.variables_bounds(ctd)
        true_u_var = Inf*ones((N+1)*(ocp.state_dimension+1+ocp.control_dimension))
        true_u_var[end-N:1:end] .= umax
        @test u_var == true_u_var
        true_l_var = -true_u_var
        @test l_var == true_l_var
    end
    @testset verbose = true showtiming = true "control" begin
        remove_constraint!(ocp, :control_con1)
        constraint!(ocp, :control, u -> u, -umax, umax, :control_con2)
        N  = 3
        ctd = CTDirect.CTDirect_data(ocp, N, nothing)
        lb, ub = CTDirect.constraints_bounds(ctd)
        println("lb = ", lb)
        println("ub = ", ub)
        #= true_u_var = Inf*ones((N+1)*(ocp.state_dimension+1+ocp.control_dimension))
        @test u_var == true_u_var
        true_l_var = -true_u_var
        @test l_var == true_l_var =#

    end
    @testset verbose = true showtiming = true "state" begin
        remove_constraint!(ocp, :control_con2)
        umax = 10*ones(ocp.state_dimension)
        constraint!(ocp, :state, -umax, umax, :state_con3)
        N  = 3
        ctd = CTDirect.CTDirect_data(ocp, N, nothing)

        println("ctd.state_box[2] = ", ctd.state_box[2])
        @test ctd.has_state_box == true
        @test ctd.dim_state_box == 2
        
        l_var, u_var = CTDirect.variables_bounds(ctd)
        println("l_var = ", l_var)
        println("u_var = ", u_var)
        #= true_u_var = Inf*ones((N+1)*(ocp.state_dimension+1+ocp.control_dimension))
        @test u_var == true_u_var
        true_l_var = -true_u_var
        @test l_var == true_l_var =#

    end
end


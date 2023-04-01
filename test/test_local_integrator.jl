# double integrator - energy min
println("Local Double integrator test")

# OCP model
n=2
m=1
t0=0.0
tf=1.0
x0=[-1.0, 0.0]
xf=[0.0, 0.0]
ocp = Model()
state!(ocp, n)   # dimension of the state
control!(ocp, m) # dimension of the control
time!(ocp, [t0, tf])
constraint!(ocp, :initial, x0)
constraint!(ocp, :final,   xf)
A = [ 0.0 1.0
    0.0 0.0 ]
B = [ 0.0
    1.0 ]
constraint!(ocp, :dynamics, (x, u) -> A*x + B*u[1])
objective!(ocp, :lagrange, (x, u) -> 0.5u[1]^2) # default is to minimise

# OCP solution
a = x0[1]
b = x0[2]
C = [-(tf-t0)^3/6.0 (tf-t0)^2/2.0
    -(tf-t0)^2/2.0 (tf-t0)]
D = [-a-b*(tf-t0), -b]+xf
p0 = C\D
α = p0[1]
β = p0[2]
x(t) = [a+b*(t-t0)+β*(t-t0)^2/2.0-α*(t-t0)^3/6.0, b+β*(t-t0)-α*(t-t0)^2/2.0]
p(t) = [α, -α*(t-t0)+β]
u(t) = [p(t)[2]]
objective = 0.5*(α^2*(tf-t0)^3/3+β^2*(tf-t0)-α*β*(tf-t0)^2)
N=201
times = range(t0, tf, N)
sol_ocp = OptimalControlSolution()
sol_ocp.state_dimension = n
sol_ocp.control_dimension = m
sol_ocp.times = times
sol_ocp.state = x
sol_ocp.state_labels = [ "x" * ctindices(i) for i ∈ range(1, n)]
sol_ocp.adjoint = p
sol_ocp.control = u
sol_ocp.control_labels = [ "u" ]
sol_ocp.objective = objective
sol_ocp.iterations = 0
sol_ocp.stopping = :dummy
sol_ocp.message = "analytical solution"
sol_ocp.success = true


# test retrieved OCP parameters
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
        @test ctd.dim_NLP_constraints == 3*N+5
    end

    # test retrieved OCP constraints
    @testset verbose = true showtiming = true "constraints" begin
        lb, ub = CTDirect.constraints_bounds(ctd)
        l_var, u_var = CTDirect.variables_bounds(ctd)
        true_lb = zeros(3*N)
        @test lb[1:3*N] == true_lb
        @test ub[1:3*N] == true_lb
        @test l_var == -Inf*ones((N+1)*4)
        @test u_var == -l_var
    end


    # test solution
    @testset verbose = true showtiming = true "numerical solution" begin

        sol_direct = solve(ocp, grid_size=100, print_level=0) 

        u_sol(t) = sol_ocp.control(t)[1]
        u = t -> sol_direct.control(t)[1]
        x_sol(t) = sol_ocp.state(t)
        x = t -> sol_direct.state(t)[1:2]
        T = sol_direct.times
        dT = T[2:end]-T[1:end-1]
        N = length(T)
        @test distance_infty(x,x_sol,T) ≈ 0 atol = 0.01
        @test distance_L2(u,u_sol,T) ≈ 0 atol=1e-1
        @test sol_direct.objective ≈ sol_ocp.objective atol=1e-2

    end


end



using CTDirect
using CTBase

println("Test: goddard all constraints")

# goddard max final altitue (all constraint types formulation)
ocp = Model(variable=true)
Cd = 310
Tmax = 3.5
β = 500
b = 2
r0 = 1
v0 = 0
vmax = 0.1
m0 = 1
mf = 0.6
x0 = [ r0, v0, m0 ]
state!(ocp, 3)
control!(ocp, 1)
variable!(ocp, 1)
time!(ocp, 0, Index(1))
# use all possible types of constraints
# initial condition
constraint!(ocp, :initial, x0, :initial_constraint)
# final condition
constraint!(ocp, :final, Index(3), mf, :final_constraint)
# state constraint
constraint!(ocp, :state, (x,v)->x[2], -Inf, vmax, :state_con_v_ub)
# control constraint
constraint!(ocp, :control, (u,v)->u, -Inf, 1, :control_con_u_ub)
# mixed constraint
constraint!(ocp, :mixed, (x,u,v)->x[3], mf, Inf, :mixed_con_m_lb)
# variable constraint
constraint!(ocp, :variable, v->v, -Inf, 10, :variable_con_tf_ubx)
# state box
constraint!(ocp, :state, 1:2, [r0,v0], [r0+0.2, Inf], :state_box_rv)
#constraint!(ocp, :state, Index(2), v0, Inf, :state_box_vmin)
# control box
constraint!(ocp, :control, Index(1), 0, Inf, :control_box_lb)
# variable box
constraint!(ocp, :variable, Index(1), 0.01, Inf, :variable_box_tfmin)
objective!(ocp, :mayer,  (x0, xf, v) -> xf[1], :max)
function F0(x)
    r, v, m = x
    D = Cd * v^2 * exp(-β*(r - 1))
    return [ v, -D/m - 1/r^2, 0 ]
end
function F1(x)
    r, v, m = x
    return [ 0, Tmax/m, -b*Tmax ]
end
dynamics!(ocp, (x, u, v) -> F0(x) + u*F1(x) )
sol1 = solve(ocp, grid_size=30, print_level=0, tol=1e-8)
@testset verbose = true showtiming = true ":goddard :max_rf :all_constraints_types" begin
    @test sol1.objective ≈ 1.0125 rtol=1e-2
end

# with initial guess
x_init = [1.05, 0.1, 0.8]
u_init = 0.5
v_init = 0.1
init_constant = OptimalControlInit(x_init, u_init, v_init)
sol2 = solve(ocp, grid_size=30, print_level=0, tol=1e-8, init=init_constant)
@testset verbose = true showtiming = true ":goddard :max_rf :all_constraints_types :init" begin
    @test sol2.objective ≈ 1.0125 rtol=1e-2
end

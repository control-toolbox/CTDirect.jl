using CTDirect

println("Test: all constraint types")
# +++ todo: redo this one with different constraint formulations, use default initial guess only
# goddard max final altitude (all constraint types formulation)
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
function FF0(x)
    r, v, m = x
    D = Cd * v^2 * exp(-β*(r - 1))
    return [ v, -D/m - 1/r^2, 0 ]
end
function FF1(x)
    r, v, m = x
    return [ 0, Tmax/m, -b*Tmax ]
end
dynamics!(ocp, (x, u, v) -> FF0(x) + u*FF1(x) )

@testset verbose = true showtiming = true ":goddard :max_rf :all_constraints" begin
    sol1 = solve(ocp, grid_size=100, print_level=0, tol=1e-8)
    @test sol1.objective ≈ 1.0125 rtol=1e-2
end

# goddard with speed constraint - maximize altitude
println("Local Goddard test")

# OCP model
Cd = 310
Tmax = 3.5
β = 500
b = 2
t0 = 0
r0 = 1
v0 = 0
vmax = 0.1
m0 = 1
mf = 0.6
x0 = [ r0, v0, m0 ]
ocp = Model()
time!(ocp, :initial, t0)
state!(ocp, 3, ["r", "v", "m"])
control!(ocp, 1)
constraint!(ocp, :initial, x0, :initial_constraint)
constraint!(ocp, :final, Index(3), mf, :final_constraint)
constraint!(ocp, :state, x->x[2], -Inf, 0.1, :state_con_vmax)
constraint!(ocp, :control, u->u, -Inf, 1, :control_con_umax)
constraint!(ocp, :mixed, (x,u)->x[3], mf, Inf, :mixed_con_mmin)
constraint!(ocp, :state, Index(1), r0, Inf, :state_box_rmin)
constraint!(ocp, :state, Index(2), v0, Inf, :state_box_vmin)
constraint!(ocp, :control, Index(1), 0, Inf, :control_box_umin)
objective!(ocp, :mayer,  (t0, x0, tf, xf) -> xf[1], :max)
function F0(x)
    r, v, m = x
    D = Cd * v^2 * exp(-β*(r - 1))
    return [ v, -D/m - 1/r^2, 0 ]
end
function F1(x)
    r, v, m = x
    return [ 0, Tmax/m, -b*Tmax ]
end
f(x, u) = F0(x) + u*F1(x)
constraint!(ocp, :dynamics, f)


# test
init = [1.01, 0.05, 0.8, 0.1]
sol = solve(ocp, grid_size=10, print_level=0, init=init)
@testset verbose = true showtiming = true ":goddard :all_constraints" begin
    @test sol.objective ≈ 1.012 atol=5e-3
end
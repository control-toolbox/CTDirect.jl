# goddard max final altitude (all constraint types formulation)
println("Test: initial guess options")
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
# control box
constraint!(ocp, :control, Index(1), 0, Inf, :control_box_lb)
# variable box
constraint!(ocp, :variable, Index(1), 0.01, Inf, :variable_box_tfmin)
objective!(ocp, :mayer,  (x0, xf, v) -> xf[1], :max)
function FFF0(x)
    r, v, m = x
    D = Cd * v^2 * exp(-β*(r - 1))
    return [ v, -D/m - 1/r^2, 0 ]
end
function FFF1(x)
    r, v, m = x
    return [ 0, Tmax/m, -b*Tmax ]
end
dynamics!(ocp, (x, u, v) -> FFF0(x) + u*FFF1(x) )
sol0 = solve(ocp, print_level=0)

# default init
@testset verbose = true showtiming = true ":default_init" begin
    sol = solve(ocp, print_level=0)
    @test sol.objective ≈ 1.0125 rtol=1e-2
end

# constant initial guess
x_const = [1.05, 0.2, 0.8]
u_const = 0.5
v_const = 0.15

@testset verbose = true showtiming = true ":constant_init_x" begin
    sol = solve(ocp, print_level=0, init=OCPInit(state=x_const))
    @test sol.objective ≈ 1.0125 rtol=1e-2
end
@testset verbose = true showtiming = true ":constant_init_u" begin
    sol = solve(ocp, print_level=0, init=OCPInit(control=u_const))
    @test sol.objective ≈ 1.0125 rtol=1e-2
end
@testset verbose = true showtiming = true ":constant_init_v" begin
    sol = solve(ocp, print_level=0, init=OCPInit(variable=v_const))
    @test sol.objective ≈ 1.0125 rtol=1e-2
end
@testset verbose = true showtiming = true ":constant_init_xu" begin
    sol = solve(ocp, print_level=0, init=OCPInit(state=x_const, control=u_const))
    @test sol.objective ≈ 1.0125 rtol=1e-2
end
@testset verbose = true showtiming = true ":constant_init_xv" begin
    sol = solve(ocp, print_level=0, init=OCPInit(state=x_const, variable=v_const))
    @test sol.objective ≈ 1.0125 rtol=1e-2
end
@testset verbose = true showtiming = true ":constant_init_uv" begin
    sol = solve(ocp, print_level=0, init=OCPInit(control=u_const, variable=v_const))
    @test sol.objective ≈ 1.0125 rtol=1e-2
end
@testset verbose = true showtiming = true ":constant_init_xuv" begin
    sol = solve(ocp, print_level=0, init=OCPInit(state=x_const, control=u_const, variable=v_const))
    @test sol.objective ≈ 1.0125 rtol=1e-2
end
@testset verbose = true showtiming = true ":constant_init_xuv :compact" begin
    sol = solve(ocp, print_level=0, init=(state=x_const, control=u_const, variable=v_const))
    @test sol.objective ≈ 1.0125 rtol=1e-2
end

# functional initial guess
x_func = t->[1+t^2, sqrt(t), 1-t]
u_func = t->(cos(t)+1)*0.5

@testset verbose = true showtiming = true ":functional_init_x" begin
    sol = solve(ocp, print_level=0, init=OCPInit(state=x_func))
    @test sol.objective ≈ 1.0125 rtol=1e-2
end
@testset verbose = true showtiming = true ":functional_init_u" begin
    sol = solve(ocp, print_level=0, init=OCPInit(control=u_func))
    @test sol.objective ≈ 1.0125 rtol=1e-2
end
@testset verbose = true showtiming = true ":functional_init_xu" begin
    sol = solve(ocp, print_level=0, init=OCPInit(state=x_func, control=u_func))
    @test sol.objective ≈ 1.0125 rtol=1e-2
end
@testset verbose = true showtiming = true ":functional_init_xu :compact" begin
    sol = solve(ocp, print_level=0, init=(state=x_func, control=u_func))
    @test sol.objective ≈ 1.0125 rtol=1e-2
end
@testset verbose = true showtiming = true ":mixed_init" begin
    sol = solve(ocp, print_level=0, init=OCPInit(state=x_func, control=u_const))
    @test sol.objective ≈ 1.0125 rtol=1e-2
end
@testset verbose = true showtiming = true ":mixed_init :compact" begin
    sol = solve(ocp, print_level=0, init=(state=x_func, control=u_const))
    @test sol.objective ≈ 1.0125 rtol=1e-2
end

# warm start
@testset verbose = true showtiming = true ":warm_start" begin
    sol = solve(ocp, print_level=0, init=OCPInit(sol0))
    @test sol.objective ≈ 1.0125 rtol=1e-2
end
@testset verbose = true showtiming = true ":warm_start :compact" begin
    sol = solve(ocp, print_level=0, init=sol0)
    @test sol.objective ≈ 1.0125 rtol=1e-2
end

# set initial guess in DOCP
docp = directTranscription(ocp)
@testset verbose = true showtiming = true ":DOCPInit_mixed" begin
    setDOCPInit(docp, OCPInit(state=x_func, control=u_const))
    dsol = solve(docp, print_level=0)
    sol = OCPSolutionFromDOCP(docp, dsol)
    @test sol.objective ≈ 1.0125 rtol=1e-2
end
@testset verbose = true showtiming = true ":DOCPInit_mixed :compact" begin
    setDOCPInit(docp, (state=x_func, control=u_const))
    dsol = solve(docp, print_level=0)
    sol = OCPSolutionFromDOCP(docp, dsol)
    @test sol.objective ≈ 1.0125 rtol=1e-2
end
@testset verbose = true showtiming = true ":DOCPInit_warm_start" begin
    setDOCPInit(docp, OCPInit(sol0))
    dsol = solve(docp, print_level=0)
    sol = OCPSolutionFromDOCP(docp, dsol)
    @test sol.objective ≈ 1.0125 rtol=1e-2
end
@testset verbose = true showtiming = true ":DOCPInit_warm_start :compact" begin
    setDOCPInit(docp, sol0)
    dsol = solve(docp, print_level=0)
    sol = OCPSolutionFromDOCP(docp, dsol)
    @test sol.objective ≈ 1.0125 rtol=1e-2
end

# pass initial guess to solve
setDOCPInit(docp, OCPInit()) # reset init in docp
@testset verbose = true showtiming = true ":solve_mixed_init" begin
    dsol = solve(docp, init=OCPInit(state=x_func, control=u_const), print_level=0)
    sol = OCPSolutionFromDOCP(docp, dsol)
    @test sol.objective ≈ 1.0125 rtol=1e-2
end
@testset verbose = true showtiming = true ":solve_mixed_init :compact" begin
    dsol = solve(docp, init=(state=x_func, control=u_const), print_level=0)
    sol = OCPSolutionFromDOCP(docp, dsol)
    @test sol.objective ≈ 1.0125 rtol=1e-2
end
@testset verbose = true showtiming = true ":solve_warm_start" begin
    dsol = solve(docp, init=OCPInit(sol0), print_level=0)
    sol = OCPSolutionFromDOCP(docp, dsol)
    @test sol.objective ≈ 1.0125 rtol=1e-2
end
@testset verbose = true showtiming = true ":solve_warm_start :compact" begin
    dsol = solve(docp, init=sol0, print_level=0)
    sol = OCPSolutionFromDOCP(docp, dsol)
    @test sol.objective ≈ 1.0125 rtol=1e-2
end
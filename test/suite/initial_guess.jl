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
time!(ocp, t0=0, indf=1)
constraint!(ocp, :initial, lb=x0, ub=x0)
constraint!(ocp, :final, rg=3, lb=mf, ub=mf)
constraint!(ocp, :state, lb=[r0,v0,mf], ub=[r0+0.2,vmax,m0])
constraint!(ocp, :control, lb=0, ub=1)
constraint!(ocp, :variable, lb=0.01, ub=Inf)
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

@testset verbose = true showtiming = true ":constant_init_x :compact" begin
    sol = solve(ocp, print_level=0, init=(state=x_const,))
    @test sol.objective ≈ 1.0125 rtol=1e-2
end
@testset verbose = true showtiming = true ":constant_init_u :compact" begin
    sol = solve(ocp, print_level=0, init=(control=u_const,))
    @test sol.objective ≈ 1.0125 rtol=1e-2
end
@testset verbose = true showtiming = true ":constant_init_v :compact" begin
    sol = solve(ocp, print_level=0, init=(variable=v_const,))
    @test sol.objective ≈ 1.0125 rtol=1e-2
end
@testset verbose = true showtiming = true ":constant_init_xu :compact" begin
    sol = solve(ocp, print_level=0, init=(state=x_const, control=u_const))
    @test sol.objective ≈ 1.0125 rtol=1e-2
end
@testset verbose = true showtiming = true ":constant_init_xv :compact" begin
    sol = solve(ocp, print_level=0, init=(state=x_const, variable=v_const))
    @test sol.objective ≈ 1.0125 rtol=1e-2
end
@testset verbose = true showtiming = true ":constant_init_uv :compact" begin
    sol = solve(ocp, print_level=0, init=(control=u_const, variable=v_const))
    @test sol.objective ≈ 1.0125 rtol=1e-2
end

@testset verbose = true showtiming = true ":constant_init_xuv :compact" begin
    sol = solve(ocp, print_level=0, init=(state=x_const, control=u_const, variable=v_const))
    @test sol.objective ≈ 1.0125 rtol=1e-2
end

# functional initial guess
x_func = t->[1+t^2, sqrt(t), 1-t]
u_func = t->(cos(t)+1)*0.5

@testset verbose = true showtiming = true ":functional_init_x :compact" begin
    sol = solve(ocp, print_level=0, init=(state=x_func,))
    @test sol.objective ≈ 1.0125 rtol=1e-2
end
@testset verbose = true showtiming = true ":functional_init_u :compact" begin
    sol = solve(ocp, print_level=0, init=(control=u_func,))
    @test sol.objective ≈ 1.0125 rtol=1e-2
end

@testset verbose = true showtiming = true ":functional_init_xu :compact" begin
    sol = solve(ocp, print_level=0, init=(state=x_func, control=u_func))
    @test sol.objective ≈ 1.0125 rtol=1e-2
end

@testset verbose = true showtiming = true ":mixed_init :compact" begin
    sol = solve(ocp, print_level=0, init=(state=x_func, control=u_const))
    @test sol.objective ≈ 1.0125 rtol=1e-2
end

@testset verbose = true showtiming = true ":warm_start :compact" begin
    sol = solve(ocp, print_level=0, init=sol0)
    @test sol.objective ≈ 1.0125 rtol=1e-2
end

# set initial guess in DOCP
docp = directTranscription(ocp)

@testset verbose = true showtiming = true ":DOCPInit_mixed :compact" begin
    setInitialGuess(docp, (state=x_func, control=u_const))
    dsol = solve(docp, print_level=0)
    sol = OCPSolutionFromDOCP(docp, dsol)
    @test sol.objective ≈ 1.0125 rtol=1e-2
end

@testset verbose = true showtiming = true ":DOCPInit_warm_start :compact" begin
    setInitialGuess(docp, sol0)
    dsol = solve(docp, print_level=0)
    sol = OCPSolutionFromDOCP(docp, dsol)
    @test sol.objective ≈ 1.0125 rtol=1e-2
end

# pass initial guess to solve
setInitialGuess(docp, OptimalControlInit()) # reset init in docp

@testset verbose = true showtiming = true ":solve_mixed_init :compact" begin
    dsol = solve(docp, init=(state=x_func, control=u_const), print_level=0)
    sol = OCPSolutionFromDOCP(docp, dsol)
    @test sol.objective ≈ 1.0125 rtol=1e-2
end

@testset verbose = true showtiming = true ":solve_warm_start :compact" begin
    dsol = solve(docp, init=sol0, print_level=0)
    sol = OCPSolutionFromDOCP(docp, dsol)
    @test sol.objective ≈ 1.0125 rtol=1e-2
end
# tests some objective options, variable tf
println("Test: objective")

# min tf
ocp = Model(variable=true)
state!(ocp, 2)
control!(ocp, 1)
variable!(ocp, 1)
time!(ocp, t0=0, indf=1)
constraint!(ocp, :initial, lb=[0,0], ub=[0,0])
constraint!(ocp, :final, lb=[1,0], ub=[1,0])
constraint!(ocp, :control, lb=-1, ub=1)
constraint!(ocp, :variable, lb=0.1, ub=10)
dynamics!(ocp, (x, u, v) ->  [x[2], u])
objective!(ocp, :mayer, (x0, xf, v) -> v)

@testset verbose = true showtiming = true ":min_tf :mayer" begin
    sol = solve(ocp, display=false, tol=1e-12)
    @test sol.objective ≈ 2.0 rtol=1e-2
end


# min tf (lagrange)
ocp = Model(variable=true)
state!(ocp, 2)
control!(ocp, 1)
variable!(ocp, 1)
time!(ocp, t0=0, indf=1)
constraint!(ocp, :initial, lb=[0,0], ub=[0,0])
constraint!(ocp, :final, lb=[1,0], ub=[1,0])
constraint!(ocp, :control, lb=-1, ub=1)
constraint!(ocp, :variable, lb=0.1, ub=10)
dynamics!(ocp, (x, u, v) ->  [x[2], u])
objective!(ocp, :lagrange, (x, u, v) -> 1)

@testset verbose = true showtiming = true ":min_tf :lagrange" begin
    sol = solve(ocp, display=false, tol=1e-12)
    @test sol.objective ≈ 2.0 rtol=1e-2
end


# max t0 (free t0 and tf)
ocp = Model(variable=true)
state!(ocp, 2)
control!(ocp, 1)
variable!(ocp, 2)
time!(ocp, ind0=1, indf=2)
constraint!(ocp, :initial, lb=[0,0], ub=[0,0])
constraint!(ocp, :final, lb=[1,0], ub=[1,0])
constraint!(ocp, :control, lb=-1, ub=1)
constraint!(ocp, :variable, lb=[0.1,0.1], ub=[10,10])
constraint!(ocp, :variable, f=v->v[2]-v[1], lb=0.1, ub=Inf)
dynamics!(ocp, (x, u, v) ->  [x[2], u])
objective!(ocp, :mayer, (x0, xf, v) -> v[1], :max)

@testset verbose = true showtiming = true ":max_t0" begin
    sol = solve(ocp, display=false, tol=1e-12)
    @test sol.objective ≈ 8.0 rtol=1e-2
end



# bolza, non-autonomous mayer term, tf in dynamics
@def ocp begin
    tf ∈ R, variable
    t ∈ [0, tf], time
    x ∈ R, state
    u ∈ R, control
    ẋ(t) == tf * u(t)
    x(0) == 0
    x(tf) == 1
    tf + 0.5∫(u(t)^2) → min
end
@testset verbose = true showtiming = true ":bolza :tf_in_dyn_and_cost" begin
    sol = solve(ocp, display=false)
    @test sol.objective ≈ 1.476 rtol=1e-2
    @test sol.variable[1] ≈ 1.107 rtol=1e-2
end

#=
@def ocp begin
    v = (t0, tf) ∈ R^2, variable
    t ∈ [t0, tf], time
    x ∈ R, state
    u ∈ R, control
    ẋ(t) == tf * u(t) + t0
    x(t0) == 0
    x(tf) == 1
    0 ≤ t0 ≤ 10
    0.01 ≤ tf - t0 ≤ 10
    (t0^2 + tf) + 0.5∫(u(t)^2) → min
end
#@testset verbose = true showtiming = true ":bolza :t0_tf_in_dyn_and_cost" begin
    sol = solve(ocp, print_level=5)
#    @test sol.variable[1] ≈ 1.107 rtol=1e-2
#end
=#

#=
@def ocp2 begin
    s ∈ [0, 1], time
    y ∈ R^2, state
    u ∈ R, control
    ẏ(s) == [u(s), 0]
    y[1](0) == 0
    y[1](1) == 0
    ∫(u(s)^2) → min
end
sol = solve(ocp2)
=#
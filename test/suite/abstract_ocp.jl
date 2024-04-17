using CTDirect

println("Test: abstract OCP definition")

# double integrator min tf, abstract definition
@def ocp1 begin
    tf ∈ R, variable
    t ∈ [ 0, tf ], time
    x ∈ R², state
    u ∈ R, control
    -1 ≤ u(t) ≤ 1
    x(0) == [ 0, 0 ]
    x(tf) == [ 1, 0 ]
    0.1 ≤ tf ≤ Inf 
    ẋ(t) == [ x₂(t), u(t) ] 
    tf → min
end

@testset verbose = true showtiming = true ":double_integrator :min_tf :abstract" begin
    sol1 = solve(ocp1, grid_size=100, print_level=0, tol=1e-12)
    @test sol1.objective ≈ 2.0 rtol=1e-2
end

# same with some random constraints
@def ocp2 begin
    tf ∈ R, variable
    t ∈ [ 0, tf ], time
    x ∈ R², state
    u ∈ R, control
    tf ≥ 0.1
    -1 ≤ u(t) ≤ 1
    q = x₁
    v = x₂
    q(0) == 1
    v(0) == 2
    q(tf) == 0
    v(tf) == 0
    0 ≤ q(t) ≤ 5
    -2 ≤ v(t) ≤ 3
    (u^2)(t) ≤ 100
    ẋ(t) == [ v(t), u(t) ]
    tf → min
end

@testset verbose = true showtiming = true ":double_integrator :min_tf :abstract :constr" begin
    sol2 = solve(ocp2, grid_size=100, print_level=0, tol=1e-12)
    @test sol2.objective ≈ 5.46 rtol=1e-2
end
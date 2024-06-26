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
    sol1 = solve(ocp1, print_level=0, tol=1e-12)
    @test sol1.objective ≈ 2.0 rtol=1e-2
end

# goddard
# NB. the ≤ is not the same as <= (parse error for <=)
@def ocp3 begin
    tf ∈ R, variable
    t ∈ [ 0, tf ], time
    x ∈ R^3, state
    u ∈ R, control
    0.1 ≤ tf ≤ Inf
    r = x[1]
    v = x[2]
    m = x[3]
    x(0) == [1, 0, 1]
    m(tf) == 0.6
    1 ≤ r(t) ≤ 1.1
    0 ≤ v(t) ≤ vmax
    0 ≤ u(t) ≤ 1
    ẋ(t) == F0(x(t)) + u(t)*F1(x(t))
    r(tf) → max
end

@testset verbose = true showtiming = true ":goddard :max_rf :abstract :constr" begin
    sol3 = solve(ocp3, print_level=0, tol=1e-12)    
    @test sol3.objective ≈ 1.0125 rtol=1e-2
end
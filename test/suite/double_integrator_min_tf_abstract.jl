using CTDirect
using CTBase

# min tf, abstract definition
@def ocp begin
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
sol = solve(ocp, grid_size=100, print_level=0, tol=1e-12)
@testset verbose = true showtiming = true ":double_integrator :min_tf :abstract" begin
    @test sol.objective ≈ 2.0 rtol=1e-2
end

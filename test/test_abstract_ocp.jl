using CTDirect
using CTBase

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

sol1 = solve(ocp1, print_level=0, tol=1e-12)
println("Target 2.0, found ", sol1.objective)


# goddard
include("problems/goddard.jl")
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

sol3 = solve(ocp3, print_level=0, tol=1e-12)    
println("Target 1.0125, found ", sol3.objective)
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

sol2 = solve(ocp2, print_level=0, tol=1e-12, init=sol1)
println("Target 5.46, found ", sol2.objective)

# goddard
Cd = 310
β = 500
Tmax = 3.5
b = 2
vmax = 0.1
function F0(x)
    r, v, m = x
    D = Cd * v^2 * exp(-β*(r - 1))
    return [ v, -D/m - 1/r^2, 0 ]
end
function F1(x)
    r, v, m = x
    return [ 0, Tmax/m, -b*Tmax ]
end
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
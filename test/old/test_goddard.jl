# goddard with state constraint - maximize altitude
println("Goddard test")
prob = Problem(:goddard, :all_constraint)
init = [1.01, 0.05, 0.8, 0.1]
sol = solve(prob.model, grid_size = 10, print_level = 0, init = init)
@testset verbose = true showtiming = true ":goddard :all_constraints" begin
    @test sol.objective ≈ prob.solution.objective rtol = 1e-2
end

# goddard, abstract def
const Cd = 310
const Tmax = 3.5
const β = 500
const b = 2
const t0 = 0
r0 = 1
v0 = 0
vmax = 0.1
m0 = 1
mf = 0.6
x0 = [r0, v0, m0]

@def ocp begin

    tf, variable
    t ∈ [t0, tf], time
    x ∈ R³, state
    u ∈ R, control

    r = x₁
    v = x₂
    m = x₃

    x(t0) == [r0, v0, m0]
    0 ≤ u(t) ≤ 1
    r0 ≤ r(t) ≤ Inf, (1)
    0 ≤ v(t) ≤ vmax, (2)
    mf ≤ m(t) ≤ m0, (3)

    ẋ(t) == F0(x(t)) + u(t) * F1(x(t))

    r(tf) → max

end

F0(x) = begin
    r, v, m = x
    D = Cd * v^2 * exp(-β * (r - 1))
    F = [v, -D / m - 1 / r^2, 0]
    return F
end

F1(x) = begin
    r, v, m = x
    F = [0, Tmax / m, -b * Tmax]
    return F
end

sol = solve(ocp, grid_size = 10, print_level = 5)
@testset verbose = true showtiming = true ":goddard :all_constraints" begin
    @test sol.objective ≈ prob.solution.objective rtol = 1e-2
end

# double integrator

function double_integrator_a()
    @def ocp begin
        tf ∈ R, variable
        t ∈ [0, tf], time
        x ∈ R², state
        u ∈ R, control
        -1 ≤ u(t) ≤ 1
        x(0) == [0, 0]
        x(tf) == [1, 0]
        0.05 ≤ tf ≤ Inf
        ẋ(t) == [x₂(t), u(t)]
        tf → min
    end

    return ((ocp = ocp, obj = 2.0, name = "double_integrator_a", init = nothing))
end

function double_integrator_T(T)
    @def ocp begin
        t ∈ [0, T], time
        x ∈ R², state
        u ∈ R, control
        q = x₁
        v = x₂
        q(0) == 0
        v(0) == 0
        q(T) == 1
        v(T) == 0
        ẋ(t) == [v(t), u(t)]
        ∫(u(t)^2) → min
    end

    return ((ocp = ocp, obj = nothing, name = "double_integrator_T", init = nothing))
end

# min tf
function double_integrator_mintf(; lagrange = false)
    ocp = Model(variable = true)
    state!(ocp, 2)
    control!(ocp, 1)
    variable!(ocp, 1)
    time!(ocp, t0 = 0, indf = 1)
    constraint!(ocp, :initial, val = [0, 0])
    constraint!(ocp, :final, val = [1, 0])
    constraint!(ocp, :control, lb = -1, ub = 1)
    constraint!(ocp, :variable, lb = 0.1, ub = 10)
    dynamics!(ocp, (x, u, v) -> [x[2], u])
    if lagrange
        objective!(ocp, :lagrange, (x, u, v) -> 1)
        name = "double_integrator_lagrange"
    else
        objective!(ocp, :mayer, (x0, xf, v) -> v)
        name = "double_integrator_mayer"
    end

    return ((ocp = ocp, obj = 2.0, name = name, init = nothing))
end

# max t0 with free t0,tf
function double_integrator_freet0tf(lagrange = false)
    ocp = Model(variable = true)
    state!(ocp, 2)
    control!(ocp, 1)
    variable!(ocp, 2)
    time!(ocp, ind0 = 1, indf = 2)
    constraint!(ocp, :initial, val = [0, 0])
    constraint!(ocp, :final, val = [1, 0])
    constraint!(ocp, :control, lb = -1, ub = 1)
    constraint!(ocp, :variable, lb = [0.1, 0.1], ub = [10, 10])
    constraint!(ocp, :variable, f = v -> v[2] - v[1], lb = 0.1, ub = Inf)
    dynamics!(ocp, (x, u, v) -> [x[2], u])
    objective!(ocp, :mayer, (x0, xf, v) -> v[1], :max)

    return ((ocp = ocp, obj = 8.0, name = "double_integrator_freet0tf", init = nothing))
end

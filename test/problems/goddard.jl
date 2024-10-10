# Goddard problem
# free final time, fixed final mass, max altitude
# constraint on max speed

# aux functions
# NB defining these inside the problem function does not seem to change the allocations
function F0(x, Cd, beta)
    r, v, m = x
    D = Cd * v^2 * exp(-beta * (r - 1))
    return [v, -D / m - 1 / r^2, 0]
end
function F1(x, Tmax, b)
    r, v, m = x
    return [0, Tmax / m, -b * Tmax]
end

function goddard(; vmax = 0.1, Tmax = 3.5, functional_constraints = false)
    # constants
    Cd = 310
    beta = 500
    b = 2
    r0 = 1
    v0 = 0
    m0 = 1
    mf = 0.6
    x0 = [r0, v0, m0]

    #ocp
    goddard = Model(variable = true)
    state!(goddard, 3)
    control!(goddard, 1)
    variable!(goddard, 1)
    time!(goddard, t0 = 0, indf = 1)
    constraint!(goddard, :initial, val = x0)
    constraint!(goddard, :final, rg = 3, val = mf)
    if functional_constraints
        # note: the equations do not handle r<1 well
        # without the box constraint on x, the default init (0.1) is not suitable
        constraint!(goddard, :state, f = (x, v) -> x, lb = [r0, v0, mf], ub = [r0 + 0.2, vmax, m0])
        constraint!(goddard, :control, f = (u, v) -> u, lb = 0, ub = 1)
    else
        constraint!(goddard, :state, lb = [r0, v0, mf], ub = [r0 + 0.1, vmax, m0])
        constraint!(goddard, :control, lb = 0, ub = 1)
    end
    constraint!(goddard, :variable, lb = 0.01, ub = Inf)
    objective!(goddard, :mayer, (x0, xf, v) -> xf[1], :max)
    dynamics!(goddard, (x, u, v) -> F0(x, Cd, beta) + u * F1(x, Tmax, b))

    return ((ocp = goddard, obj = 1.01257, name = "goddard", init = (state = [1.01, 0.05, 0.8],)))
end

# abstract definition
function goddard_a(; vmax = 0.1, Tmax = 3.5)
    # constants
    Cd = 310
    beta = 500
    b = 2
    r0 = 1
    v0 = 0
    m0 = 1
    mf = 0.6
    x0 = [r0, v0, m0]

    @def goddard_a begin
        tf ∈ R, variable
        t ∈ [0, tf], time
        x ∈ R^3, state
        u ∈ R, control
        0.01 ≤ tf ≤ Inf
        r = x[1]
        v = x[2]
        m = x[3]
        x(0) == x0
        m(tf) == mf
        r0 ≤ r(t) ≤ r0 + 0.1
        v0 ≤ v(t) ≤ vmax
        mf ≤ m(t) ≤ m0
        0 ≤ u(t) ≤ 1
        ẋ(t) == F0(x(t), Cd, beta) + u(t) * F1(x(t), Tmax, b)
        r(tf) → max
    end

    return ((
        ocp = goddard_a,
        obj = 1.01257,
        name = "goddard_a",
        init = (state = [1.01, 0.05, 0.8],),
    ))
end



# all constraint types formulation
function goddard_all()
    # constants
    Cd = 310
    beta = 500
    b = 2
    r0 = 1
    v0 = 0
    m0 = 1
    mf = 0.6
    x0 = [r0, v0, m0]
    vmax = 0.1
    Tmax = 3.5

    #ocp
    goddard = Model(variable = true)
    state!(goddard, 3)
    control!(goddard, 1)
    variable!(goddard, 1)
    time!(goddard, t0 = 0, indf = 1)
    # initial constraint
    constraint!(goddard, :initial, val = x0)
    # final constraint
    constraint!(goddard, :final, rg = 3, val = mf)
    # state box (active at t0 and tf)
    constraint!(goddard, :state, lb = [r0, v0, 0], ub = [Inf, Inf, m0])
    # control box (active on last bang arc)
    constraint!(goddard, :control, lb = 0, ub = Inf)
    # variable box (inactive)
    constraint!(goddard, :variable, lb = 0.01, ub = Inf)
    # state constraint (active on constrained arc)
    constraint!(goddard, :state, f = (x, v) -> x[2], lb = -Inf, ub = vmax)
    # control constraint (active on first bang arc)
    constraint!(goddard, :control, f = (u, v) -> u, lb = -Inf, ub = 1)
    # 'mixed' constraint (active at tf)
    constraint!(goddard, :mixed, f = (x, u, v) -> x[3], lb = mf, ub = Inf)
    objective!(goddard, :mayer, (x0, xf, v) -> xf[1], :max)
    #dynamics!(goddard, (x, u, v) ->  F0(x, Cd, beta) .+ u .* F1(x, Tmax, b)) slightly better
    dynamics!(goddard, (x, u, v) ->  F0(x, Cd, beta) + u * F1(x, Tmax, b))

    return ((
        ocp = goddard,
        obj = 1.01257,
        name = "goddard_all",
        init = (state = [1.01, 0.05, 0.8],),
    ))
end

# inplace version
function goddard_all_inplace()
    # constants
    Cd = 310
    beta = 500
    b = 2
    r0 = 1
    v0 = 0
    m0 = 1
    mf = 0.6
    x0 = [r0, v0, m0]
    vmax = 0.1
    Tmax = 3.5

    function v_fun!(c, x, v)
        @views c[:] .= x[2]
        return
    end
    function u_fun!(c, u, v)
        @views c[:] .= u
        return
    end
    function m_fun!(c, x, u, v) 
        #@views c[:] .= x[3]
        c[1] = x[3]
        return
    end
    function rf_fun!(c, x0, xf, v)
        @views c[:] .= xf[1]
        return
    end
    function f_fun!(f, x, u, v)
        @views f[:] .= F0(x, Cd, beta) + u * F1(x, Tmax, b)
        return
    end

    #ocp
    goddard = Model(variable = true, in_place=true)
    state!(goddard, 3)
    control!(goddard, 1)
    variable!(goddard, 1)
    time!(goddard, t0 = 0, indf = 1)
    # initial constraint
    constraint!(goddard, :initial, val = x0)
    # final constraint
    constraint!(goddard, :final, rg = 3, val = mf)
    # state box (active at t0 and tf)
    constraint!(goddard, :state, lb = [r0, v0, 0], ub = [Inf, Inf, m0])
    # control box (active on last bang arc)
    constraint!(goddard, :control, lb = 0, ub = Inf)
    # variable box (inactive)
    constraint!(goddard, :variable, lb = 0.01, ub = Inf)
    # state constraint (active on constrained arc)
    constraint!(goddard, :state, f = v_fun!, lb = -Inf, ub = vmax)
    # control constraint (active on first bang arc)
    constraint!(goddard, :control, f = u_fun!, lb = -Inf, ub = 1)
    # 'mixed' constraint (active at tf)
    constraint!(goddard, :mixed, f = m_fun!, lb = mf, ub = Inf)
    objective!(goddard, :mayer, rf_fun!, :max)
    dynamics!(goddard, f_fun!)

    return ((
        ocp = goddard,
        obj = 1.01257,
        name = "goddard_all_inplace",
        init = (state = [1.01, 0.05, 0.8],),
    ))
end

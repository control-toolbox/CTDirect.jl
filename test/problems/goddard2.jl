# Goddard2 problem (version for :exa)
# free final time, fixed final mass, max altitude
# constraint on max speed

function goddard2(; vmax=0.1, Tmax=3.5)
    # constants
    Cd = 310
    β = 500
    b = 2
    r0 = 1
    v0 = 0
    m0 = 1
    mf = 0.6
    x0 = [r0, v0, m0]

    goddard2 = @def begin
        tf ∈ R, variable
        t ∈ [0, tf], time
        x = (r, v, m) ∈ R^3, state
        u ∈ R, control
        tf ≥ 0.01
        x(0) == x0
        m(tf) == mf
        r0 ≤ r(t) ≤ r0 + 0.1
        v0 ≤ v(t) ≤ vmax
        mf ≤ m(t) ≤ m0
        0 ≤ u(t) ≤ 1

        ∂(r)(t) == v(t)
        ∂(v)(t) ==
        -Cd * v(t)^2 * exp(-β * (r(t) - 1)) / m(t) - 1 / r(t)^2 + u(t) * Tmax / m(t)
        ∂(m)(t) == -b * Tmax * u(t)

        r(tf) → max
    end

    return ((ocp=goddard2, obj=1.01257, name="goddard2", init=(state=[1.01, 0.05, 0.8],)))
end

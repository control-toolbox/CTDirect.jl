# Schlogl problem

function schlogl()

    # constants
    t0 = 0
    x0 = [1]
    k0 = 6
    k1 = 11
    k2 = 6
    k3 = 1

    @def schlogl begin
        T ∈ R, variable
        t ∈ [t0, T], time
        x ∈ R, state
        u = (u0, u1, u2, u3) ∈ R⁴, control

        x(t0) == x0
        x(T) == [2]
        u0(t) ≥ 0.1
        u1(t) ≥ 0.1
        u2(t) ≥ 0.1
        u3(t) ≥ 0.1
        1 ≥ T ≥ 0.02
        x(t) ≥ 0.5

        ẋ(t) == u0(t) - u1(t) + u2(t) - u3(t)

        ∫(
            u0(t)*log(abs(u0(t)/k0)) - (u0(t)-k0) +
            u1(t)*log(abs(u1(t)/(k1*x(t)))) - (u1(t)-k1*x(t)) +
            u2(t)*log(abs(u2(t)/(k2*x(t)^2))) - (u0(t)-k2*x(t)^2) +
            u3(t)*log(abs(u3(t)/(k3*x(t)^3))) - (u3(t)-k3*x(t)^3)
        ) → min
    end

    return ((ocp=schlogl, obj=nothing, name="schlogl", init=nothing))
end

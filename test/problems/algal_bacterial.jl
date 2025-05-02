# Optimal control of an algal-bacterial consortium system. Original code from Rand Asswad.

function algal_bacterial()

    # parameters
    s_in = 0.5
    β = 23e-3
    γ = 0.44
    dmax = 1.5
    ϕmax = 6.48;
    ks = 0.09;
    ρmax = 27.3e-3;
    kv = 0.57e-3;
    μmax = 1.0211;
    qmin = 2.7628e-3;
    ϕ(s) = ϕmax * s / (ks + s)
    ρ(v) = ρmax * v / (kv + v)
    μ(q) = μmax * (1 - qmin / q)
    t0 = 0;
    tf = 20
    x0 = [0.1629, 0.0487, 0.0003, 0.0177, 0.035, 0]

    # to be closer to jump formulation; similar performance to inlined expression in def
    function f(x, α, d)
        return [
            d*(s_in - x[1]) - ϕ(x[1])*x[2]/γ,           # s
            ((1-α)*ϕ(x[1]) - d) * x[2],                   # e
            α * β * ϕ(x[1]) * x[2] - ρ(x[3])*x[5] - d*x[3],   # v
            ρ(x[3]) - μ(x[4])*x[4],                     # q
            (μ(x[4]) - d) * x[5],                         # c
            d * x[5]                                    # obj = d*c
        ]
    end

    @def algal_bacterial begin
        t ∈ [t0, tf], time
        x ∈ R⁶, state
        u ∈ R², control

        x(t0) == x0
        x(t) ≥ [0, 0, 0, qmin, 0, 0]
        [0, 0] ≤ u(t) ≤ [1, dmax]

        ẋ(t) == f(x(t), u₁(t), u₂(t))

        x₆(tf) → max
    end

    return ((ocp=algal_bacterial, obj=5.45, name="algal_bacterial", init=nothing))
end

"""
Space Shuttle Reentry Trajectory Problem:
    We want to find the optimal trajectory of a space shuttle reentry.
    The objective is to maximize the latitude (cross range) at the terminal point.
    The original problem formulated as a JuMP model can be found [here](https://jump.dev/JuMP.jl/stable/tutorials/nonlinear/space_shuttle_reentry_trajectory/)
    Note: no heating limit path constraint
"""
function space_shuttle()
    
    ## Global variables
    w = 203000.0  # weight (lb)
    g₀ = 32.174    # acceleration (ft/sec^2)
    m = w / g₀    # mass (slug)

    ## Aerodynamic and atmospheric forces on the vehicle
    ρ₀ = 0.002378
    hᵣ = 23800.0
    Rₑ = 20902900.0
    μ = 0.14076539e17
    S = 2690.0
    a₀ = -0.20704
    a₁ = 0.029244
    b₀ = 0.07854
    b₁ = -0.61592e-2
    b₂ = 0.621408e-3

    ## Initial conditions
    h_s = 2.6          # altitude (ft) / 1e5
    ϕ_s = deg2rad(0)   # longitude (rad)
    θ_s = deg2rad(0)   # latitude (rad)
    v_s = 2.56         # velocity (ft/sec) / 1e4
    γ_s = deg2rad(-1)  # flight path angle (rad)
    ψ_s = deg2rad(90)  # azimuth (rad)
    α_s = deg2rad(0)   # angle of attack (rad)
    β_s = deg2rad(0)   # bank angle (rad)

    ## Final conditions, the so-called Terminal Area Energy Management (TAEM)
    h_t = 0.8          # altitude (ft) / 1e5
    v_t = 0.25         # velocity (ft/sec) / 1e4
    γ_t = deg2rad(-5)  # flight path angle (rad)

    ## dynamics
    function dynamics(x, u)
        scaled_h, ϕ, θ, scaled_v, γ, ψ = x
        α, β = u
        h = scaled_h * 1e5
        v = scaled_v * 1e4
        ## Helper functions
        c_D = b₀ + b₁ * rad2deg(α) + b₂ * (rad2deg(α)^2)
        c_L = a₀ + a₁ * rad2deg(α)
        ρ = ρ₀ * exp(-h / hᵣ)
        D = (1 / 2) * c_D * S * ρ * (v^2)
        L = (1 / 2) * c_L * S * ρ * (v^2)
        r = Rₑ + h
        g = μ / (r^2)

        ## dynamics  
        h_dot = v * sin(γ)
        ϕ_dot = (v / r) * cos(γ) * sin(ψ) / cos(θ)
        θ_dot = (v / r) * cos(γ) * cos(ψ)
        v_dot = -(D / m) - g * sin(γ)
        γ_dot = (L / (m * v)) * cos(β) + cos(γ) * ((v / r) - (g / v))
        ψ_dot =
            (1 / (m * v * cos(γ))) * L * sin(β) +
            (v / (r * cos(θ))) * cos(γ) * sin(ψ) * sin(θ)

        return [h_dot / 1e5, ϕ_dot, θ_dot, v_dot / 1e4, γ_dot, ψ_dot]
    end

    ocp = @def begin
  
        ## define the problem
        tf ∈ R¹, variable 
        t ∈ [0, tf], time
        x ∈ R⁶, state
        u ∈ R², control

        ## state variables
        scaled_h = x₁
        ϕ = x₂
        θ = x₃
        scaled_v = x₄
        γ = x₅
        ψ = x₆

        ## control variables
        α = u₁
        β = u₂

        ## constraints
        1750 ≤ tf ≤ 2250 # NB jump with 503 steps between 3.5 and 4.5
        # state constraints
        0 ≤ scaled_h(t) ≤ Inf, (scaled_h_con)
        deg2rad(-89) ≤ θ(t) ≤ deg2rad(89), (θ_con)
        0 ≤ scaled_v(t) ≤ Inf, (scaled_v_con)
        deg2rad(-89) ≤ γ(t) ≤ deg2rad(89), (γ_con)
        # control constraints
        deg2rad(-90) ≤ α(t) ≤ deg2rad(90), (α_con)
        deg2rad(-89) ≤ β(t) ≤ deg2rad(1), (β_con)

        # initial conditions
        scaled_h(0) == h_s, (scaled_h0_con)
        ϕ(0) == ϕ_s, (ϕ0_con)
        θ(0) == θ_s, (θ0_con)
        scaled_v(0) == v_s, (scaled_v0_con)
        γ(0) == γ_s, (γ0_con)
        ψ(0) == ψ_s, (ψ0_con)
        # final conditions
        scaled_h(tf) == h_t, (scaled_hf_con)
        scaled_v(tf) == v_t, (scaled_vf_con)
        γ(tf) == γ_t, (γf_con)

        ## dynamics  
        ẋ(t) == dynamics(x(t), u(t))

        ## objective
        θ(tf) → max
    end

    # initial guess: linear interpolation for h, v, gamma (NB. t0 = 0), constant for the rest
    # variable time step seems to be initialized at 1 in jump
    # note that ipopt will project the initial guess inside the bounds anyway.
    tf_init = 500
    x_init = t -> [ h_s + t / tf_init * (h_t - h_s) ,
    ϕ_s,
    θ_s,
    v_s + t / tf_init * (v_t - v_s),
    γ_s + t / tf_init * (γ_t - γ_s),
    ψ_s]
    init = (state=x_init, control=[α_s, β_s], variable=[tf_init])

    # objective should be 34.18deg (0.5966rad) with tf = 2009
    return ((ocp=ocp, obj=deg2rad(34.18), name="space_shuttle", init=init))
end

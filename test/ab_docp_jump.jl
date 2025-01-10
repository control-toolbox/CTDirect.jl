# Jump version

import Pkg
Pkg.add.(["JuMP", "Ipopt"])
using JuMP, Ipopt

function algal_bacterial_jump(;N=5000, disc_method=:trapeze)

    #initialize JuMP model with Ipopt solver backend
    sys = JuMP.Model(Ipopt.Optimizer)
    set_optimizer_attribute(sys, "print_level", 4)
    set_optimizer_attribute(sys, "tol", 1e-8)
    set_optimizer_attribute(sys, "max_iter", 1500)
    set_optimizer_attribute(sys, "mu_strategy", "adaptive")

    # Discretization parameters
    t0 = 0; tf = 20
    Δt = (tf - t0) / N

    # System parameters
    s_in = 0.5
    β = 23e-3
    γ = 0.44
    dmax = 1.5

    ϕmax = 6.48; ks = 0.09;
    ρmax = 27.3e-3; kv = 0.57e-3;
    μmax = 1.0211; qmin = 2.7628e-3;

    # Variables
    @variables(sys, begin
        s[1:N+1] ≥ 0                  # s
        e[1:N+1] ≥ 0                  # e
        v[1:N+1] ≥ 0                  # v
        q[1:N+1] ≥ qmin               # q
        c[1:N+1] ≥ 0                  # c
        g[1:N+1] ≥ 0
        0.0 ≤ α[1:N+1] ≤ 1.0          # α
        0.0 ≤ d[1:N+1] ≤ dmax         # d
    end)

    # Objective
    @objective(sys, Max, g[N+1])

    # Boundary constraints: initial condition
    x0 = [0.1629, 0.0487, 0.0003, 0.0177, 0.035]
    @constraints(sys, begin
        con_s0, s[1] == x0[1]
        con_e0, e[1] == x0[2]
        con_v0, v[1] == x0[3]
        con_q0, q[1] == x0[4]
        con_c0, c[1] == x0[5]
        con_obj, g[1] == 0
    end)

    # Dynamics
    @NLexpressions(sys, begin
        # rate functions
        ϕ_se[i = 1:N+1], ϕmax * s[i] * e[i] / (ks + s[i])
        ρ_v[i = 1:N+1], ρmax * v[i] / (kv + v[i])
        μ_q[i = 1:N+1], μmax * (1 - qmin / q[i])
        # dynamics
        ds[i = 1:N+1], d[i]*(s_in - s[i]) - ϕ_se[i]/γ
        de[i = 1:N+1], (1-α[i])*ϕ_se[i] - d[i]*e[i]
        dv[i = 1:N+1], α[i]*β*ϕ_se[i] - ρ_v[i]*c[i] - d[i]*v[i]
        dq[i = 1:N+1], ρ_v[i] - μ_q[i]*q[i]
        dc[i = 1:N+1], (μ_q[i] - d[i])*c[i]
        # objective dynamics
        dg[i = 1:N+1], d[i] * c[i]
    end)

    # Crank-Nicolson scheme (aka implicit trapezoidal rule)
    if disc_method == :trapeze
        @NLconstraints(sys, begin
            con_ds[i = 1:N], s[i+1] == s[i] + Δt * (ds[i] + ds[i+1])/2.0
            con_de[i = 1:N], e[i+1] == e[i] + Δt * (de[i] + de[i+1])/2.0
            con_dv[i = 1:N], v[i+1] == v[i] + Δt * (dv[i] + dv[i+1])/2.0
            con_dq[i = 1:N], q[i+1] == q[i] + Δt * (dq[i] + dq[i+1])/2.0
            con_dc[i = 1:N], c[i+1] == c[i] + Δt * (dc[i] + dc[i+1])/2.0
            con_dg[i = 1:N], g[i+1] == g[i] + Δt * (dg[i] + dg[i+1])/2.0
        end)
    else
        error("unknown disc method: ", disc_method)
    end

    # Optimization 
    println("Solving...")
    optimize!(sys)
    println("Objective value = ", objective_value(sys), "\n")
    
end

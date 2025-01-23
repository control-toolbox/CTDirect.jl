# Jump version for algal bacterial problem
# original code from Rand Asswad 
# https://gist.github.com/rand-asswad/6443f65cd1e8704221925b886eebd12b
# https://gist.github.com/rand-asswad/7120f1ab39bf28fffd414e780b4c0196
using JuMP, Ipopt


struct rk_method
    name::Symbol
    s::Integer
    a::Matrix{<:Real}
    b::Vector{<:Real}
    c::Vector{<:Real}
end

rk = rk_method(:gauss2, 2, 
    [0.25 (0.25 - sqrt(3)/6); (0.25 + sqrt(3)/6) 0.25],
    [0.5, 0.5],
    [(0.5 - sqrt(3)/6), (0.5 + sqrt(3)/6)]
)

function algal_bacterial_jump(;grid_size=1000, disc_method=:trapeze, print_level=5)

    #initialize JuMP model with Ipopt solver backend
    sys = JuMP.Model(Ipopt.Optimizer)
    set_optimizer_attribute(sys, "print_level", print_level)
    set_optimizer_attribute(sys, "tol", 1e-8)
    set_optimizer_attribute(sys, "max_iter", 1500)
    set_optimizer_attribute(sys, "mu_strategy", "adaptive")
    set_optimizer_attribute(sys, "sb", "yes")

    # Discretization parameters
    N = grid_size
    t0 = 0; tf = 20
    Δt = (tf - t0) / N

    # System parameters
    s_in = 0.5
    β = 23e-3
    γ = 0.44
    dmax = 1.5
    ϕmax = 6.48
    ks = 0.09
    ρmax = 27.3e-3
    kv = 0.57e-3
    μmax = 1.0211
    qmin = 2.7628e-3
    ϕ(s) = ϕmax * s / (ks + s)
    ρ(v) = ρmax * v / (kv + v)
    μ(q) = μmax * (1 - qmin / q)

    x_lower = [0, 0, 0, qmin, 0, 0]     # lower bound for x
    x0 = [0.1629, 0.0487, 0.0003, 0.0177, 0.035, 0.]

    if disc_method == :trapeze
        # Crank-Nicolson scheme (aka implicit trapezoidal rule)

        # Variables
        @variables(sys, begin
            s[1:N+1] ≥ x_lower[1]       
            e[1:N+1] ≥ x_lower[2]   
            v[1:N+1] ≥ x_lower[3]      
            q[1:N+1] ≥ x_lower[4]     
            c[1:N+1] ≥ x_lower[5]       
            g[1:N+1] ≥ x_lower[6]
            0.0 ≤ α[1:N+1] ≤ 1.0 
            0.0 ≤ d[1:N+1] ≤ dmax  
        end)

        # Objective
        @objective(sys, Max, g[N+1])

        # Boundary constraints: initial condition
        @constraints(sys, begin
            con_s0, s[1] == x0[1]
            con_e0, e[1] == x0[2]
            con_v0, v[1] == x0[3]
            con_q0, q[1] == x0[4]
            con_c0, c[1] == x0[5]
            con_obj, g[1] == x0[6]
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

        @NLconstraints(sys, begin
            con_ds[i = 1:N], s[i+1] == s[i] + Δt * (ds[i] + ds[i+1])/2.0
            con_de[i = 1:N], e[i+1] == e[i] + Δt * (de[i] + de[i+1])/2.0
            con_dv[i = 1:N], v[i+1] == v[i] + Δt * (dv[i] + dv[i+1])/2.0
            con_dq[i = 1:N], q[i+1] == q[i] + Δt * (dq[i] + dq[i+1])/2.0
            con_dc[i = 1:N], c[i+1] == c[i] + Δt * (dc[i] + dc[i+1])/2.0
            con_dg[i = 1:N], g[i+1] == g[i] + Δt * (dg[i] + dg[i+1])/2.0
        end)

    elseif disc_method == :gauss_legendre_2
        # Gauss Legendre 2

        # Variables
        n = 5
        @variables(sys, begin
            x[1:N, i=1:n+1] ≥ x_lower[i]    # x
            0 ≤ α[1:N] ≤ 1                  # α
            0 ≤ d[1:N] ≤ dmax               # d
            k[1:rk.s, 1:N, 1:n+1]           # k (for Runge-Kutta)
        end)

        # Objective function
        @objective(sys, Max, x[end, end])

        # Initial condition
        @constraint(sys, initial, x[1,:] == x0[:])

        # Autonomous dynamics function
        function f(x, α, d)
            return [
                d*(s_in - x[1]) - ϕ(x[1])*x[2]/γ,           # s
                ((1-α)*ϕ(x[1]) - d)*x[2],                   # e
                α*β*ϕ(x[1])*x[2] - ρ(x[3])*x[5] - d*x[3],   # v
                ρ(x[3]) - μ(x[4])*x[4],                     # q
                (μ(x[4]) - d)*x[5],                         # c
                d * x[5]                                    # obj = d*c
            ]
        end

        # Runge-Kutta methods for autonomous systems (as a nonlinear constraint)
        # x[i+1] = x[i] + Δt Σ_j b[j]k[j,i]
        # k[j,i] = f( x[i] + Δt Σ_s A[j,s]k[s,i] )
        @constraints(sys, begin
            rk_nodes[j=1:rk.s, i=1:N], k[j,i,:] == f(x[i,:] + Δt*sum(rk.a[j,s]*k[s,i,:] for s in 1:rk.s), α[i], d[i])
            rk_scheme[i = 1:N-1], x[i+1,:] == x[i,:] + Δt*sum(rk.b[j] * k[j,i,:] for j in 1:rk.s)
        end)

    else
        error("unknown disc method: ", disc_method)
    end

    # Optimization 
    print_level > 0 && println("Solving...")
    optimize!(sys)
    print_level > 0 && println("Objective value = ", objective_value(sys), "\n")

    # todo: add optional plot if graphs=true is passed

end

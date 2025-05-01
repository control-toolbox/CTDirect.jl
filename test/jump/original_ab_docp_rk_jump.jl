# Standalone script to calculate optimal control of an algal-bacterial consortium system
# using the direct method, using Julia's JuMP package with Ipopt backend
# for nonlinear optimization.
# This script is based on the article: https://doi.org/10.48550/arXiv.2212.03157
# Caillau et al. (2022) - An An algorithmic guide for finite-dimensional optimal control problems

using JuMP, Ipopt, Plots, LaTeXStrings

# Define Runge-Kutta Butcher tabeau data structure
struct rk_method
    name::Symbol
    s::Integer
    a::Matrix{<:Real}
    b::Vector{<:Real}
    c::Vector{<:Real}
end

# Define Gauss-Legendre method (of order 2s) with s=2
# https://en.wikipedia.org/wiki/Gauss-Legendre_method
rk = rk_method(:gauss2, 2,
    [0.25 (0.25 - sqrt(3)/6); (0.25 + sqrt(3)/6) 0.25],
    [0.5, 0.5],
    [(0.5 - sqrt(3)/6), (0.5 + sqrt(3)/6)]
)

# Initialize JuMP model with Ipopt solver backend 
docp = JuMP.Model(Ipopt.Optimizer)
set_optimizer_attribute(docp, "print_level", 5)
set_optimizer_attribute(docp, "tol", 1e-8)
set_optimizer_attribute(docp, "max_iter", 1500)
set_optimizer_attribute(docp, "mu_strategy", "adaptive")

# Discretization parameters
N = 5000    # time steps
t0 = 0;
tf = 20
Δt = (tf - t0) / N

# System definition
n = 5       # dim(x)

# System parameters
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

# Declare constrained variables
x_lower = [0, 0, 0, qmin, 0, 0]     # lower bound for x
@variables(docp, begin
    x[1:(N+1), i=1:(n+1)] ≥ x_lower[i]    # x
    0 ≤ α[1:N] ≤ 1                  # α
    0 ≤ d[1:N] ≤ dmax               # d
    k[1:rk.s, 1:N, 1:(n+1)]           # k (for Runge-Kutta)
end)

# Objective function
@objective(docp, Max, x[end, end])

# Initial condition
x0 = [0.1629, 0.0487, 0.0003, 0.0177, 0.035, 0]
@constraint(docp, initial, x[1, :] == x0[:])

# Autonomous dynamics function
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

# Runge-Kutta methods for autonomous systems (as a nonlinear constraint)
# x[i+1] = x[i] + Δt Σ_j b[j]k[j,i]
# k[j,i] = f( x[i] + Δt Σ_s A[j,s]k[s,i] )
@constraints(docp, begin
    rk_nodes[j=1:rk.s, i=1:N], k[j, i, :] == f(x[i, :] + Δt*sum(rk.a[j, s]*k[s, i, :] for s in 1:rk.s), α[i], d[i])
    rk_scheme[i=1:N], x[i+1, :] == x[i, :] + Δt*sum(rk.b[j] * k[j, i, :] for j in 1:rk.s)
end)

# Optimization --------------------------------------------------------------

# Solution data structure
struct OCSolution
    t::Vector{Real}
    x::Matrix{Real}
    λ::Matrix{Real}
    u::Matrix{Real}
    obj::Real
end

@timev begin
    # Solve DOCP as nonlinear optimization problem
    optimize!(docp)
    x = value.(docp.obj_dict[:x])
    a = value.(docp.obj_dict[:α])
    d = value.(docp.obj_dict[:d])
    λ = dual.(docp.obj_dict[:rk_scheme])

    # Retrieve OCP solution data
    t = (0:(N-1))*Δt
    u = transpose(vcat(transpose.((a, d))...));
    λ = -reduce(vcat, transpose.(λ))
    sol = OCSolution(t, x[1:N, :], λ, u, x[end, end])
end

# Plots ---------------------------------------------------------------------

latexify(tab) = latexstring.(L"$" .* tab .* L"$")
x = ["s", "e", "v", "q", "c"]
u = [raw"\alpha", "d"]
ind = [1, 3, 4, 2, 5]
l = @layout [° ° °; ° °]

function plot_state(sol::OCSolution)
    plot(sol.t, sol.x[:, ind], label=latexify(reshape(x[ind], (1, 5))), layout=l)
    plot!(fontfamily="sans-serif")
    xlabel!("Time")
end

function plot_costate(sol::OCSolution)
    plot(sol.t, sol.λ[:, ind], label=latexify(reshape("\\lambda_" .* x[ind], (1, 5))), layout=l)
    xlabel!("Time")
    plot!(fontfamily="sans-serif")
end

function plot_control(sol::OCSolution)
    p1 = plot(sol.t, sol.u[:, 1], label=latexstring(u[1]))
    p2 = plot(sol.t, sol.u[:, 2], label=latexstring(u[2]))
    plot(p1, p2, layout=(2, 1))
end

function plot_objective(sol::OCSolution)
    plot(sol.t, sol.x[:, end], fillrange=zeros(length(sol.t)), alpha=0.5, legend=false)
    mid = Integer(ceil(length(sol.t) * 0.75))
    xpos = sol.t[mid]
    ypos = sol.x[mid, end] / 2
    plot!(annotations=(xpos, ypos, Plots.text(round(sol.obj; digits=3), :hcenter)))
end

# save plots
fname(s) = "test/docs/jump_gl2_" * string(N) * "_" * s * ".pdf"
function save_plots(sol)
    savefig(plot_state(sol), fname("state"))
    savefig(plot_costate(sol), fname("costate"))
    savefig(plot_control(sol), fname("control"))
    savefig(plot_objective(sol), fname("objective"))
end

save_plots(sol)

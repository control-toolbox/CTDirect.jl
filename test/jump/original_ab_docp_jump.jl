# Standalone script to calculate optimal control of an algal-bacterial consortium system
# using the direct method, using Julia's JuMP package with Ipopt backend
# for nonlinear optimization.
# This script is based on the article: https://doi.org/10.48550/arXiv.2212.03157
# Caillau et al. (2022) - An An algorithmic guide for finite-dimensional optimal control problems

using JuMP, Ipopt, Plots, LaTeXStrings

#initialize JuMP model with Ipopt solver backend
sys = JuMP.Model(Ipopt.Optimizer)
set_optimizer_attribute(sys, "print_level", 5)
set_optimizer_attribute(sys, "tol", 1e-8)
set_optimizer_attribute(sys, "max_iter", 1500)
set_optimizer_attribute(sys, "mu_strategy", "adaptive")

# Discretization parameters
N = 5000
t0 = 0;
tf = 20
Δt = (tf - t0) / N

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

# Variables
@variables(sys, begin
    s[1:(N+1)] ≥ 0                  # s
    e[1:(N+1)] ≥ 0                  # e
    v[1:(N+1)] ≥ 0                  # v
    q[1:(N+1)] ≥ qmin               # q
    c[1:(N+1)] ≥ 0                  # c
    g[1:(N+1)] ≥ 0
    0.0 ≤ α[1:(N+1)] ≤ 1.0          # α
    0.0 ≤ d[1:(N+1)] ≤ dmax         # d
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
    ϕ_se[i=1:(N+1)], ϕmax * s[i] * e[i] / (ks + s[i])
    ρ_v[i=1:(N+1)], ρmax * v[i] / (kv + v[i])
    μ_q[i=1:(N+1)], μmax * (1 - qmin / q[i])
    # dynamics
    ds[i=1:(N+1)], d[i]*(s_in - s[i]) - ϕ_se[i]/γ
    de[i=1:(N+1)], (1-α[i])*ϕ_se[i] - d[i]*e[i]
    dv[i=1:(N+1)], α[i]*β*ϕ_se[i] - ρ_v[i]*c[i] - d[i]*v[i]
    dq[i=1:(N+1)], ρ_v[i] - μ_q[i]*q[i]
    dc[i=1:(N+1)], (μ_q[i] - d[i])*c[i]
    # objective dynamics
    dg[i=1:(N+1)], d[i] * c[i]
end)

# Crank-Nicolson scheme (aka implicit trapezoidal rule)
@NLconstraints(sys, begin
    con_ds[i=1:N], s[i+1] == s[i] + Δt * (ds[i] + ds[i+1])/2.0
    con_de[i=1:N], e[i+1] == e[i] + Δt * (de[i] + de[i+1])/2.0
    con_dv[i=1:N], v[i+1] == v[i] + Δt * (dv[i] + dv[i+1])/2.0
    con_dq[i=1:N], q[i+1] == q[i] + Δt * (dq[i] + dq[i+1])/2.0
    con_dc[i=1:N], c[i+1] == c[i] + Δt * (dc[i] + dc[i+1])/2.0
    con_dg[i=1:N], g[i+1] == g[i] + Δt * (dg[i] + dg[i+1])/2.0
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

function solve!(sys)
    # Solves for the control and state
    println("Solving...")
    optimize!(sys)

    # Print optimization status
    if termination_status(sys) == MOI.OPTIMAL
        println("Solution is optimal")
    elseif termination_status(sys) == MOI.LOCALLY_SOLVED
        println("Local solution found")
    elseif termination_status(sys) == MOI.TIME_LIMIT && has_values(sys)
        println("Solution is suboptimal due to a time limit, but a primal solution is available")
    else
        println("The model was not solved correctly.")
    end
    println("Objective value = ", objective_value(sys), "\n")

    # retrieve state, costate, and control trajectories
    t = (1:(N+1))*Δt;
    t = (t[1:(end-1)] + t[2:end])/2.0;
    ss = value.(s);
    ss = (ss[1:(end-1)] + ss[2:end])/2.0;
    ee = value.(e);
    ee = (ee[1:(end-1)] + ee[2:end])/2.0;
    vv = value.(v);
    vv = (vv[1:(end-1)] + vv[2:end])/2.0;
    qq = value.(q);
    qq = (qq[1:(end-1)] + qq[2:end])/2.0;
    cc = value.(c);
    cc = (cc[1:(end-1)] + cc[2:end])/2.0;
    gg = value.(g);
    gg = (gg[1:(end-1)] + gg[2:end])/2.0;
    x = transpose(vcat(transpose.((ss, ee, vv, qq, cc, gg))...));
    duals = (dual.(con_ds), dual.(con_de), dual.(con_dv), dual.(con_dq), dual.(con_dc))
    λ = -transpose(vcat(transpose.(duals)...));
    aa = value.(α);
    aa = (aa[1:(end-1)] + aa[2:end])/2.0;
    dd = value.(d);
    dd = (dd[1:(end-1)] + dd[2:end])/2.0;
    u = transpose(vcat(transpose.((aa, dd))...));

    OCSolution(t, x, λ, u, x[end, end])
end

@timev sol = solve!(sys)

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
#mkpath("plots/")
fname(s) = "test/docs/jump_trapeze_" * string(N) * "_" * s * ".pdf"
function save_plots(sol)
    savefig(plot_state(sol), fname("state"))
    savefig(plot_costate(sol), fname("costate"))
    savefig(plot_control(sol), fname("control"))
    savefig(plot_objective(sol), fname("objective"))
end

save_plots(sol)

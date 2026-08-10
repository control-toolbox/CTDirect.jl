#!/usr/bin/env julia
# =============================================================================
# Shape probe for CTDirect  —  "1-D = scalar" boundary audit  (roadmap #610 §7)
# =============================================================================
#
# Purpose
# -------
# Measure, against the CURRENT working-tree CTDirect (no source modification), the
# *actual* shapes CTDirect hands to the user OCP functions stored in a CTModels model:
#
#     dynamics!(r, t, x, u, v)   lagrange(t, x, u, v)   mayer(x0, xf, v)
#     path!(r, t, x, u, v)       boundary!(r, x0, xf, v)
#
# The Handbook contract (philosophy/dimension-and-shape.md) says a 1-D quantity must
# reach the user function as a `Number`, never as a length-1 vector, with the coercion
# driven by the *declared* dimension (`only` when dim == 1, `identity` otherwise).
#
# The probe is a diagnostic, not a test suite: nothing asserts pass/fail. Every
# experiment is wrapped in try/catch, and every outcome is collected into capability
# matrices printed at the bottom.
#
# How it is run
# -------------
#     julia --project=probe/shape probe/shape/probe_shape.jl
#
# Blocks
# ------
#   A. Boundary shapes per scheme     — what x/u/v types each scheme actually passes
#   B. Scalar-written user functions  — does `r[1] = -x + u` work today?
#   C. Getter-level coercion blockers — who else consumes the getters and needs vectors
#   D. Solution side                  — does the returned Solution honour 1-D = scalar?
#   E. `only` primitive               — idempotence / views / Duals / allocations
# =============================================================================

using Pkg
Pkg.activate(@__DIR__)
Pkg.develop(; path=joinpath(@__DIR__, "..", ".."))   # the repo root (this CTDirect)
Pkg.instantiate()

using Printf: @printf, @sprintf
using ForwardDiff: ForwardDiff
import CommonSolve
import CTModels
import CTSolvers
import CTDirect
import NLPModelsIpopt

# =============================================================================
# Reporting helpers
# =============================================================================

const RESULTS = Vector{NTuple{4,String}}()   # (block, row, col, cell)

record!(block, row, col, cell) = push!(RESULTS, (block, row, col, string(cell)))

"Compact type label: Number / Vector{n} / <error>."
function shape_label(x)
    x isa Number && return "Number"
    x isa AbstractVector && return "Vec($(length(x)))"
    return string(typeof(x))
end

function section(title)
    println()
    println("=" ^ 78)
    println(title)
    println("=" ^ 78)
end

function try_show(label, f)
    try
        r = f()
        @printf("  %-52s → %s\n", label, r)
        return string(r)
    catch e
        msg = first(split(sprint(showerror, e), "\n"))
        @printf("  %-52s ✗ %s\n", label, msg)
        return "✗ " * string(typeof(e))
    end
end

# =============================================================================
# Type-recording OCP models (top level: no world-age surprises)
#
# Every user function records the *types it was handed* into SEEN, and computes with
# `[1]` indexing so it stays correct whether it receives a scalar or a 1-vector
# (the Handbook's "safe common denominator").
# =============================================================================

const SEEN = Dict{Symbol,Any}()
seen!(k, v) = (SEEN[k] = v; nothing)

# --- 1-D problem: n = m = q = 1 --------------------------------------------------
dyn_1d!(r, t, x, u, v) = (seen!(:dynamics, (x, u, v)); r[1] = -x[1] + u[1] + 0 * v[1]; nothing)
lag_1d(t, x, u, v) = (seen!(:lagrange, (x, u, v)); 0.5 * u[1]^2 + 0 * x[1] + 0 * v[1])
may_1d(x0, xf, v) = (seen!(:mayer, (x0, xf, v)); 0.0 * (x0[1] + xf[1] + v[1]))
path_1d!(r, t, x, u, v) = (seen!(:path, (x, u, v)); r[1] = x[1] + u[1]; nothing)
bnd_1d!(r, x0, xf, v) = (seen!(:boundary, (x0, xf, v)); r[1] = x0[1]; r[2] = xf[1]; nothing)

# --- n-D problem: n = 2, m = 2, q = 2 --------------------------------------------
dyn_2d!(r, t, x, u, v) = (seen!(:dynamics, (x, u, v)); r[1] = x[2]; r[2] = u[1] + 0 * v[1]; nothing)
lag_2d(t, x, u, v) = (seen!(:lagrange, (x, u, v)); 0.5 * u[1]^2 + 0 * x[1] + 0 * v[1])
may_2d(x0, xf, v) = (seen!(:mayer, (x0, xf, v)); 0.0 * (x0[1] + xf[1] + v[1]))
path_2d!(r, t, x, u, v) = (seen!(:path, (x, u, v)); r[1] = x[1] + u[1]; nothing)
bnd_2d!(r, x0, xf, v) = (seen!(:boundary, (x0, xf, v)); r[1] = x0[1]; r[2] = xf[1]; nothing)

"""
Build a recording OCP with declared dimensions (n, m, q), free final time when q ≥ 1
is *not* requested (kept fixed here so the probe stays about shapes, not times).
"""
function build_recording(n::Int, m::Int, q::Int)
    pre = CTModels.PreModel()
    CTModels.time!(pre; t0=0.0, tf=1.0)
    CTModels.state!(pre, n)
    CTModels.control!(pre, m)
    q > 0 && CTModels.variable!(pre, q)
    CTModels.time_dependence!(pre; autonomous=false)
    if n == 1
        CTModels.dynamics!(pre, dyn_1d!)
        CTModels.objective!(pre, :min; lagrange=lag_1d, mayer=may_1d)
        CTModels.constraint!(pre, :path; f=path_1d!, lb=[-Inf], ub=[Inf], label=:p)
        CTModels.constraint!(
            pre, :boundary; f=bnd_1d!, lb=[-1.0, 0.0], ub=[-1.0, 0.0], label=:b
        )
    else
        CTModels.dynamics!(pre, dyn_2d!)
        CTModels.objective!(pre, :min; lagrange=lag_2d, mayer=may_2d)
        CTModels.constraint!(pre, :path; f=path_2d!, lb=[-Inf], ub=[Inf], label=:p)
        CTModels.constraint!(
            pre, :boundary; f=bnd_2d!, lb=[-1.0, 0.0], ub=[-1.0, 0.0], label=:b
        )
    end
    CTModels.definition!(pre, quote end)
    return CTModels.build(pre)
end

"Discretize `ocp` with `scheme` and return the internal DOCP."
function make_docp(ocp, scheme::Symbol; grid_size=4)
    disc = CTDirect.Collocation(; grid_size=grid_size, scheme=scheme)
    dm = CTSolvers.discretize(ocp, disc)
    return dm.cache.docp
end

"""
Drive `__objective` and `__constraints!` once and return the recorded argument shapes
per user function.
"""
function record_shapes(ocp, scheme::Symbol)
    empty!(SEEN)
    docp = make_docp(ocp, scheme)
    xu = collect(0.1 * ones(docp.dim_NLP_variables))
    CTDirect.__objective(xu, docp)
    CTDirect.__constraints!(zeros(docp.dim_NLP_constraints), xu, docp)
    return copy(SEEN)
end

const SCHEMES = (
    :euler,
    :euler_implicit,
    :trapeze,
    :midpoint,
    :gauss_legendre_2_constant_control,
    :gauss_legendre_2,
    :gauss_legendre_3,
)

const FUNS = (:dynamics, :lagrange, :mayer, :path, :boundary)

function block_A()
    section("A. Boundary shapes today — what CTDirect passes to the user functions")
    for (tag, (n, m, q)) in (("1-D  (n=m=q=1)", (1, 1, 1)), ("n-D  (n=m=q=2)", (2, 2, 2)))
        ocp = build_recording(n, m, q)
        println("\n--- $tag ---")
        for scheme in SCHEMES
            seen = try
                record_shapes(ocp, scheme)
            catch e
                @printf("  %-36s ✗ %s\n", scheme, typeof(e))
                for f in FUNS
                    record!("A/$tag", string(scheme), string(f), "✗")
                end
                continue
            end
            for f in FUNS
                got = get(seen, f, nothing)
                if got === nothing
                    record!("A/$tag", string(scheme), string(f), "—")
                    continue
                end
                a, b, c = got
                # (x, u, v) for dynamics/lagrange/path ; (x0, xf, v) for mayer/boundary
                cell = "$(shape_label(a))/$(shape_label(b))/$(shape_label(c))"
                record!("A/$tag", string(scheme), string(f), cell)
            end
            cells = join(
                [
                    "$(f)=" * last(
                        RESULTS[findlast(
                            r -> r[2] == string(scheme) && r[3] == string(f), RESULTS
                        )],
                    ) for f in FUNS
                ],
                "  ",
            )
            @printf("  %-36s %s\n", scheme, cells)
        end
    end
end

# =============================================================================
# B. Scalar-written user functions — the code the convention is supposed to enable
# =============================================================================

# Pure scalar style: NO `[1]` indexing anywhere. Correct iff the caller passes scalars.
dyn_scalar!(r, t, x, u, v) = (r[1] = -x + u; nothing)
lag_scalar(t, x, u, v) = 0.5 * u^2
bnd_scalar!(r, x0, xf, v) = (r[1] = x0; r[2] = xf; nothing)

function build_scalar_style()
    pre = CTModels.PreModel()
    CTModels.time!(pre; t0=0.0, tf=1.0)
    CTModels.state!(pre, 1)
    CTModels.control!(pre, 1)
    CTModels.time_dependence!(pre; autonomous=false)
    CTModels.dynamics!(pre, dyn_scalar!)
    CTModels.objective!(pre, :min; lagrange=lag_scalar)
    CTModels.constraint!(
        pre, :boundary; f=bnd_scalar!, lb=[-1.0, 0.0], ub=[-1.0, 0.0], label=:b
    )
    CTModels.definition!(pre, quote end)
    return CTModels.build(pre)
end

function block_B()
    section("B. Scalar-written user functions (`r[1] = -x + u`, no `[1]` on inputs)")
    println("   The convention exists to make this style legal. Measured today:")
    ocp = build_scalar_style()
    for scheme in SCHEMES
        cell = try
            docp = make_docp(ocp, scheme)
            xu = collect(0.1 * ones(docp.dim_NLP_variables))
            CTDirect.__objective(xu, docp)
            CTDirect.__constraints!(zeros(docp.dim_NLP_constraints), xu, docp)
            "✓ works"
        catch e
            "✗ " * string(typeof(e))
        end
        @printf("  %-36s %s\n", scheme, cell)
        record!("B", string(scheme), "scalar-style", cell)
    end
end

# =============================================================================
# C. Getter-level coercion blockers — other consumers of the same getters
# =============================================================================

function build_free_tf()
    pre = CTModels.PreModel()
    CTModels.variable!(pre, 1)
    CTModels.time!(pre; t0=0.0, indf=1)
    CTModels.state!(pre, 1)
    CTModels.control!(pre, 1)
    CTModels.time_dependence!(pre; autonomous=false)
    CTModels.dynamics!(pre, dyn_1d!)
    CTModels.objective!(pre, :min; mayer=may_1d)
    CTModels.definition!(pre, quote end)
    return CTModels.build(pre)
end

"""
Option A ("coerce inside the getters") would hand a scalar to every *internal* consumer
too, not just to the user functions. This block measures whether each of those consumers
survives a scalar — i.e. the true cost of Option A.
"""
function block_C()
    section("C. Option A audit — other consumers of the same getters, fed a scalar")

    ocp = try
        build_free_tf()
    catch e
        println("  (free-tf model could not be built: $(typeof(e)))")
        nothing
    end

    if ocp !== nothing
        # get_time_grid: v = get_OCP_variable(xu, docp) → CTModels.final_time(ocp, v)
        record!(
            "C",
            "get_time_grid → CTModels.final_time(ocp, v)",
            "vector v",
            try_show(
                "CTModels.final_time(ocp, [0.5])", () -> CTModels.final_time(ocp, [0.5])
            ),
        )
        record!(
            "C",
            "get_time_grid → CTModels.final_time(ocp, v)",
            "scalar v",
            try_show("CTModels.final_time(ocp, 0.5)", () -> CTModels.final_time(ocp, 0.5)),
        )
    end

    # ode/common.jl getter(): V[:, i] .= get_OCP_state_at_time_step(...)
    record!(
        "C",
        "getter(): V[:, i] .= <state>",
        "scalar v",
        try_show("V = zeros(1,3); V[:,1] .= 2.0", function ()
            V = zeros(1, 3)
            V[:, 1] .= 2.0
            return V[1, 1]
        end),
    )

    # DOCP_data.jl:553  v[:] = getter(nlp_solution; val=:variable)
    record!(
        "C",
        "build_OCP_solution: v[:] = <variable>",
        "scalar v",
        try_show("v = zeros(1); v[:] = 2.0", function ()
            v = zeros(1)
            v[:] = 2.0
            return v[1]
        end),
    )
    record!(
        "C",
        "build_OCP_solution: v[:] = <variable>",
        "vector v",
        try_show("v = zeros(1); v[:] = [2.0]", function ()
            v = zeros(1)
            v[:] = [2.0]
            return v[1]
        end),
    )

    # midpoint.jl:55-59  xs = 0.5 * (x_i + x_{i+1})
    record!(
        "C",
        "midpoint: xs = 0.5*(xi + xip1)",
        "scalar v",
        try_show("0.5 * (2.0 + 3.0)", () -> 0.5 * (2.0 + 3.0)),
    )

    # trapeze.jl:132-136  x_next = xi ; x_next += h/2 * (work_i + work_{i+1})
    record!(
        "C",
        "trapeze: x_next = xi; x_next += h*(w_i+w_ip1)",
        "scalar v",
        try_show("2.0 + 0.5 * ([1.0] + [1.0])", function ()
            x_next = 2.0
            x_next += 0.5 * ([1.0] + [1.0])
            return x_next
        end),
    )
    record!(
        "C",
        "trapeze: x_next = xi; x_next += h*(w_i+w_ip1)",
        "vector v",
        try_show("[2.0] + 0.5 * ([1.0] + [1.0])", function ()
            x_next = [2.0]
            x_next += 0.5 * ([1.0] + [1.0])
            return x_next[1]
        end),
    )

    # irk.jl:204-208 / irk_stagewise.jl:368-374 — the stage state is a WORK BUFFER,
    # not a getter result: no getter-level coercion can ever reach it.
    record!(
        "C",
        "IRK stage state: @. work_xij = xi (+ h a_jl k_il)",
        "scalar v",
        try_show("w = zeros(1); @. w = 2.0 ; w", function ()
            w = zeros(1)
            @. w = 2.0
            return w[1]
        end),
    )
    record!(
        "C",
        "IRK stage state: passed to dynamics as",
        "scalar v",
        try_show("typeof(work_xij) reaching CTModels.dynamics", () -> "Vec(1) — always"),
    )
end

# =============================================================================
# D. Solution side — does the built Solution honour 1-D = scalar?
# =============================================================================

function block_D()
    section("D. Solution side (CTModels build_solution) for a 1-D problem")
    cell = try
        pre = CTModels.PreModel()
        CTModels.time!(pre; t0=0.0, tf=1.0)
        CTModels.state!(pre, 1)
        CTModels.control!(pre, 1)
        CTModels.time_dependence!(pre; autonomous=false)
        CTModels.dynamics!(pre, (r, t, x, u, v) -> (r[1] = -x[1] + u[1]; nothing))
        CTModels.objective!(pre, :min; lagrange=(t, x, u, v) -> 0.5 * u[1]^2)
        CTModels.constraint!(
            pre,
            :boundary;
            f=(r, x0, xf, v) -> (r[1] = x0[1]; r[2] = xf[1]; nothing),
            lb=[-1.0, 0.0],
            ub=[-1.0, 0.0],
            label=:b,
        )
        CTModels.definition!(pre, quote end)
        ocp = CTModels.build(pre)

        dm = CTSolvers.discretize(ocp, CTDirect.Collocation(; grid_size=20, scheme=:midpoint))
        init = CTModels.build_initial_guess(ocp, ())
        sol = CommonSolve.solve(
            dm,
            init,
            CTSolvers.ADNLP(; backend=:optimized),
            CTSolvers.Ipopt(; max_iter=200, tol=1e-8, print_level=0, sb="yes");
            display=false,
        )
        x = CTModels.state(sol)
        u = CTModels.control(sol)
        p = CTModels.costate(sol)
        @printf("  %-46s → %s\n", "state(sol)(0.5)", shape_label(x(0.5)))
        @printf("  %-46s → %s\n", "control(sol)(0.5)", shape_label(u(0.5)))
        @printf("  %-46s → %s\n", "costate(sol)(0.5)", shape_label(p(0.5)))
        record!("D", "state(sol)(t)", "n=1", shape_label(x(0.5)))
        record!("D", "control(sol)(t)", "m=1", shape_label(u(0.5)))
        record!("D", "costate(sol)(t)", "n=1", shape_label(p(0.5)))
        "ok"
    catch e
        msg = first(split(sprint(showerror, e), "\n"))
        @printf("  solve/solution inspection ✗ %s\n", msg)
        record!("D", "solve", "n=1", "✗ " * string(typeof(e)))
        "✗"
    end
    return cell
end

# =============================================================================
# E. `only` primitive — the coercion CTDirect would use
# =============================================================================

function block_E()
    section("E. `only` as the dimension-driven coercion")
    xu = collect(1.0:10.0)
    w = @view xu[3:3]
    try_show("only(view of length 1)", () -> only(w))
    try_show("only(2.0)  (idempotent on a scalar)", () -> only(2.0))
    try_show("only([2.0])", () -> only([2.0]))
    try_show("(2.0)[1]   (scalar indexing)", () -> (2.0)[1])
    try_show("only(view) under ForwardDiff", function ()
        g = ForwardDiff.derivative(t -> only(@view [t, 0.0][1:1])^2, 3.0)
        return g
    end)
    try_show("only(empty view) — must never be reached for dim 0", () -> only(@view xu[3:2]))

    # allocation contract
    f_only(a) = only(@view a[3:3])
    f_only(xu)
    n = @allocated f_only(xu)
    @printf("  %-46s → %d bytes\n", "@allocated only(@view xu[3:3])", n)
    record!("E", "only(view)", "allocations", string(n) * " B")
end

# =============================================================================
# Matrix printing
# =============================================================================

function print_matrix()
    section("CAPABILITY MATRIX")
    blocks = unique(first.(RESULTS))
    for b in blocks
        rows = filter(r -> r[1] == b, RESULTS)
        println("\n### $b")
        cols = unique(getindex.(rows, 3))
        rownames = unique(getindex.(rows, 2))
        w = maximum(length.(rownames)) + 2
        @printf("%-*s", w, "")
        for c in cols
            @printf("%-34s", c)
        end
        println()
        for rn in rownames
            @printf("%-*s", w, rn)
            for c in cols
                idx = findlast(r -> r[2] == rn && r[3] == c, rows)
                @printf("%-34s", idx === nothing ? "" : rows[idx][4])
            end
            println()
        end
    end
end

# =============================================================================
# main
# =============================================================================

block_A()
block_B()
block_C()
block_D()
block_E()
print_matrix()

println()
println("Legend for block A cells: <arg1>/<arg2>/<arg3>")
println("  dynamics / lagrange / path : x / u / v")
println("  mayer / boundary           : x0 / xf / v")
println("Contract target (Handbook): 1-D quantity → Number; n-D → Vec(n).")

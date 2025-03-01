include("deps.jl")
using Plots
using Printf

println("Test: refined grid options")
test4 = false # use time_grid=:refined to do this internally, ie in solve call a continuation_steps(...) ?
test5 = false # time_grid=:optimized ?
test6 = false # time_grid=:optimized ?

# 4. pseudo grid refinement with manual grid input
if test4
    # here the grids are uniform so we could just pass grid_size
    # but this illustrate the possibility of grid refinement
    include("../problems/goddard.jl")
    N_target = 250

    # basic solve for comparison
    N = N_target
    println("\nBasic solve")
    @time begin
        sol = solve(goddard, time_grid=LinRange(0, 1, N + 1), display=false)
        @printf("steps %4d, objective %9.6f, iterations %4d\n", N, sol.objective, sol.iterations)
    end

    # init
    println("\nStep continuation")
    @time begin
        N = 10
        sol = solve(goddard, time_grid=LinRange(0, 1, N + 1), display=false, max_iter=5)
        @printf(
            "steps %4d, objective %9.6f, iterations set to %4d\n",
            N,
            sol.objective,
            sol.iterations
        )
        # loop on steps
        while N < N_target
            global N = min(N * 2, N_target)
            global sol = solve(
                goddard,
                time_grid=LinRange(0, 1, N + 1),
                display=false,
                init=sol,
                max_iter=5,
            )
            @printf(
                "steps %4d, objective %9.6f, iterations set to %4d\n",
                N,
                sol.objective,
                sol.iterations
            )
        end
        # final solve
        N = N_target
        sol = solve(goddard, time_grid=LinRange(0, 1, N + 1), display=false, init=sol)
        @printf("steps %4d, objective %9.6f, iterations %4d\n", N, sol.objective, sol.iterations)
    end
end

# 5. time grid as optimization variables

# NB The function below is a dirty hack.
# How is this handled by AD ?
function dt(t, v)
    if t == 0
        dt = v[1]
    elseif t == N_vars
        dt = v[end]
    else
        dt = 0.5 * (v[Int(t)] + v[Int(t)+1])
    end
    return dt
end

if test5
    # number of time steps to be optimized
    N_vars = 10

    ocp = Model(variable=true, autonomous=false)
    state!(ocp, 2)
    control!(ocp, 1)
    variable!(ocp, N_vars)
    time!(ocp, t0=0, tf=N_vars)
    constraint!(ocp, :initial, lb=[0, 0], ub=[0, 0])
    constraint!(ocp, :final, lb=[1, 0], ub=[1, 0])
    constraint!(ocp, :control, lb=-1, ub=1)
    constraint!(ocp, :variable, lb=0.01 * ones(N_vars), ub=10 * ones(N_vars))
    dynamics!(ocp, (t, x, u, v) -> [x[2], u] * dt(t, v))
    # min tf
    objective!(ocp, :mayer, (x0, xf, v) -> sum(v))
    #=    # min energy fixed tf
    fixed_tf = 2.5
    constraint!(ocp, :boundary, f=(x0, xf, v)->sum(v), lb=fixed_tf, ub=fixed_tf)
    objective!(ocp, :lagrange, (t, x, u, v) -> u^2 * v[Int(t)])
    =#

    sol = solve(ocp, grid_size=N_vars)

    # actual time grid
    v = sol.variable
    T_opt = zeros(N_vars + 1)
    for i = 1:N_vars
        T_opt[i+1] = T_opt[i] + v[i]
    end
    println("Optimized time steps ", T_opt)
    println("And tf: ", sum(sol.variable))

    #plot(sol)
    U_opt = zeros(N_vars + 1)
    # ffs julia
    for i = 1:(N_vars+1)
        U_opt[i] = sol.control(sol.times[i])
    end
    p = plot(T_opt, U_opt, markershape=:circle, show=true)
end

# goddard test case: does not work very well
# it could be the optimization finds 'bad' values for the time steps so it can cheat the ODE and get a better objective...
# also maybe a derivatives problem for function dt
# this feature needs a proper implementation anyway
if test6
    N_vars = 30
    goddard = Model(variable=true, autonomous=false)
    Cd = 310
    Tmax = 3.5
    β = 500
    b = 2
    r0 = 1
    v0 = 0
    vmax = 0.1
    m0 = 1
    mf = 0.6
    x0 = [r0, v0, m0]
    state!(goddard, 3)
    control!(goddard, 1)
    variable!(goddard, N_vars)
    time!(goddard, t0=0, tf=N_vars)
    constraint!(goddard, :initial, lb=x0, ub=x0)
    constraint!(goddard, :final, rg=3, lb=mf, ub=mf)
    constraint!(goddard, :state, lb=[r0, v0, mf], ub=[r0 + 0.2, vmax, m0])
    constraint!(goddard, :control, lb=0, ub=1)
    constraint!(goddard, :variable, lb=0.05 / N_vars * ones(N_vars), ub=Inf * ones(N_vars))
    objective!(goddard, :mayer, (x0, xf, v) -> xf[1], :max)
    function F0(x)
        r, v, m = x
        D = Cd * v^2 * exp(-β * (r - 1))
        return [v, -D / m - 1 / r^2, 0]
    end
    function F1(x)
        r, v, m = x
        return [0, Tmax / m, -b * Tmax]
    end
    dynamics!(goddard, (t, x, u, v) -> (F0(x) + u * F1(x)) * dt(t, v))

    # start will small time steps ?
    sol = solve(goddard, grid_size=N_vars, tol=1e-12, init=(variable=zeros(N_vars),))

    # actual time grid
    v = sol.variable
    T_opt = zeros(N_vars + 1)
    for i = 1:N_vars
        T_opt[i+1] = T_opt[i] + v[i]
    end
    println("Optimized time steps ", T_opt)
    println("And tf: ", sum(sol.variable))

    U_opt = sol.control.(sol.times)
    p = plot(T_opt, U_opt, markershape=:circle, show=true)
end

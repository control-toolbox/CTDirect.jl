using CTDirect
using CTProblems
using CTBase # for plot
using Plots

Plots.default(show=true)

prob = Problem(:goddard, :all_constraints)
ocp = prob.model
#println(constraints(ocp)) +++ not valid anymore 

plot_solution = true
plot_state_constraints = false
plot_control_constraints = false
plot_mixed_constraints = false
plot_state_box = false
plot_control_box = false

check_state_constraints = false
check_control_constraints = false
check_mixed_constraints = false
check_state_box = false
check_control_box = false

# solve problem
println("Solving test problem...")
init = [1.01, 0.05, 0.8, 0.1]
sol = solve(ocp, grid_size=50, print_level=0, tol=1e-8, mu_strategy="adaptive", init=init)
t0 = sol.times[1]
tf = last(sol.times)
t = sol.times
if plot_solution
    plot(sol)
end

# state contraints
if plot_state_constraints
    ncx = sol.infos[:dim_state_constraints]
    if ncx > 0
        PCX = Array{Plots.Plot,1}(undef, ncx)
        for i = 1:ncx
            PCX[i] = plot(
                t -> sol.infos[:state_constraints](t)[i],
                t0,
                tf,
                label="state_constraint v <= 0.1",
                legend=:topleft,
            )
            PCX[i] = plot!(
                twinx(),
                t -> sol.infos[:mult_state_constraints](t)[i],
                t0,
                tf,
                color=:red,
                label="multiplier",
                xticks=:none,
            ) #2nd scale
        end
        plot(PCX..., layout=(ncx, 1))
    end
end

# control constraints
if plot_control_constraints
    ncu = sol.infos[:dim_control_constraints]
    if ncu > 0
        PCU = Array{Plots.Plot,1}(undef, ncu)
        for i = 1:ncu
            PCU[i] = plot(
                t -> sol.infos[:control_constraints](t)[i],
                t0,
                tf,
                label="control_constraint u <= 1",
                legend=:topleft,
            )
            PCU[i] = plot!(
                twinx(),
                t -> sol.infos[:mult_control_constraints](t)[i],
                t0,
                tf,
                color=:red,
                label="multiplier",
                xticks=:none,
            ) #2nd scale
        end
        plot(PCU..., layout=(ncu, 1))
    end
end

# mixed constraints
if plot_mixed_constraints
    ncxu = sol.infos[:dim_mixed_constraints]
    if ncxu > 0
        PCXU = Array{Plots.Plot,1}(undef, ncxu)
        for i = 1:ncxu
            PCXU[i] = plot(
                t -> sol.infos[:mixed_constraints](t)[i],
                t0,
                tf,
                label="mixed_constraint m >= 0.6",
                legend=:topleft,
            )
            PCXU[i] = plot!(
                twinx(),
                t -> sol.infos[:mult_mixed_constraints](t)[i],
                t0,
                tf,
                color=:red,
                label="multiplier",
                xticks=:none,
            ) #2nd scale
        end
        plot(PCXU..., layout=(ncxu, 1))
    end
end

# state box  +++ not generic !
if plot_state_box
    PX = Array{Plots.Plot,1}(undef, 2) # only boxes on r, v
    for i = 1:2
        PX[i] = plot(
            t -> sol.state(t)[i],
            t0,
            tf,
            label="state with box mult (green: LB, red: UB)",
            legend=:topleft,
        )
        PX[i] = plot!(
            twinx(),
            t -> sol.infos[:mult_state_box_lower](t)[i],
            t0,
            tf,
            color=:green,
            xticks=:none,
            label=:none,
            linestyle=:dash,
        )
        PX[i] = plot!(
            twinx(),
            t -> sol.infos[:mult_state_box_upper](t)[i],
            t0,
            tf,
            color=:red,
            xticks=:none,
            label=:none,
            linestyle=:dash,
        )
    end
    PPX = plot(PX..., layout=(2, 1))
end

# control box   +++ not generic !
if plot_control_box
    PU = plot(
        t -> sol.control(t)[1],
        t0,
        tf,
        label="control with box mult (green: LB u>=0, red: UB unused)",
        legend=:topleft,
    )
    PU = plot!(
        twinx(),
        t -> sol.infos[:mult_control_box_lower](t)[1],
        t0,
        tf,
        color=:green,
        xticks=:none,
        label=:none,
        linestyle=:dash,
    )
    PU = plot!(
        twinx(),
        t -> sol.infos[:mult_control_box_upper](t)[1],
        t0,
        tf,
        color=:red,
        xticks=:none,
        label=:none,
        linestyle=:dash,
    )
    plot(PU)
end

# +++ add automatic tests for multipliers, ie sign and complementarity with constraints ?

# state constraint v <= 0.1
if check_state_constraints
    a = sol.infos[:mult_state_constraints](t)[1]
    println(a)
    println(all(>=(0), a))
    b = sol.infos[:state_constraints](t)[1]
    println(b)
    println(a .* b)
    println(norm(a .* b))
end

# control constraint u <= 1

# mixed constraint m >= 0.6

# state box r >= 1 and v >= 0
# control box u >= 0 

using CTDirect
using CTProblems
using CTBase # for plot
using Plots

# goddard with state constraint - maximize altitude
prob = Problem(:goddard, :state_constraint)
ocp = prob.model

# explicitely set specific constraints for plot
# println(constraints(ocp))
remove_constraint!(ocp,:state_constraint_r)
remove_constraint!(ocp,:state_constraint_v)
remove_constraint!(ocp,:control_constraint)
# state constraint
constraint!(ocp, :state, x->x[2], -Inf, 0.1, :state_con_vmax)
# control constraint
constraint!(ocp, :control, u->u, -Inf, 1, :control_con_umax)
# mixed constraint
constraint!(ocp, :mixed, (x,u)->x[3], 0.6, Inf, :mixed_con_mmin)
# state box
constraint!(ocp, :state, Index(1), 1, Inf, :state_box_rmin)
constraint!(ocp, :state, Index(2), 0, Inf, :state_box_vmin)
# control box
constraint!(ocp, :control, Index(1), 0, Inf, :control_box_umin)
# println(constraints(ocp))

# initial guess (constant state and control functions)
init = [1.01, 0.05, 0.8, 0.1]

# solve problem
sol = solve(ocp, grid_size=100, print_level=5, tol=1e-12, mu_strategy="adaptive", init=init)

# plot
#plot(sol)

# test plots for constraint and multipliers
t0 = sol.times[1]
tf = last(sol.times)

#=
# state constraints
ncx = sol.infos[:dim_state_constraints]
if ncx > 0
    PCX = Array{Plots.Plot, 1}(undef, ncx);
    for i in 1:ncx
        PCX[i] = plot(t -> sol.infos[:state_constraints](t)[i], t0, tf, label="state_constraint v <= 0.1", legend=:topleft)
        PCX[i] = plot!(twinx(), t -> sol.infos[:mult_state_constraints](t)[i], t0, tf, color=:red, label="multiplier", xticks=:none) #2nd scale
    end
    P1 = plot(PCX..., layout = (ncx,1));
else
    P1 = nothing
end

# control constraints
ncu = sol.infos[:dim_control_constraints]
if ncu > 0
    PCU = Array{Plots.Plot, 1}(undef, ncu);
    for i in 1:ncu
        PCU[i] = plot(t -> sol.infos[:control_constraints](t)[i], t0, tf, label="control_constraint u <= 1", legend=:topleft)
        PCU[i] = plot!(twinx(), t -> sol.infos[:mult_control_constraints](t)[i], t0, tf, color=:red, label="multiplier", xticks=:none) #2nd scale
    end
    P2 = plot(PCU..., layout = (ncu,1));
else
    P2 = nothing
end

# mixed constraints
ncxu = sol.infos[:dim_mixed_constraints]
if ncxu > 0
    PCXU = Array{Plots.Plot, 1}(undef, ncxu);
    for i in 1:ncxu
        PCXU[i] = plot(t -> sol.infos[:mixed_constraints](t)[i], t0, tf, label="mixed_constraint m >= 0.6", legend=:topleft)
        PCXU[i] = plot!(twinx(), t -> sol.infos[:mult_mixed_constraints](t)[i], t0, tf, color=:red, label="multiplier", xticks=:none) #2nd scale
    end
    P3 = plot(PCXU..., layout = (ncxu,1));
else
    P3 = nothing
end

plot(P1, P2, P3, layout = (3,1))
=#

# state box with multipliers
PX = Array{Plots.Plot, 1}(undef, 2); # only boxes on r, v
for i in 1:2
    PX[i] = plot(t -> sol.state(t)[i], t0, tf, label="state with box mult (green: LB, red: UB)", legend=:topleft)
    PX[i] = plot!(twinx(),t -> sol.infos[:mult_state_box_lower](t)[i], t0, tf, color=:green, xticks=:none, label=:none, linestyle=:dash)
    PX[i] = plot!(twinx(),t -> sol.infos[:mult_state_box_upper](t)[i], t0, tf, color=:red, xticks=:none, label=:none, linestyle=:dash)
end
PPX = plot(PX..., layout = (2, 1))

#=
# control box with multipliers
PU = plot(t -> sol.control(t)[1], t0, tf, label="control with box mult (green: LB u>=0, red: UB unused)", legend=:topleft)
PU = plot!(twinx(),t -> sol.infos[:mult_control_box_lower](t)[1], t0, tf, color=:green, xticks=:none, label=:none, linestyle=:dash)
PU = plot!(twinx(),t -> sol.infos[:mult_control_box_upper](t)[1], t0, tf, color=:red, xticks=:none, label=:none, linestyle=:dash)
plot(PU)
=#

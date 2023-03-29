using CTDirect
using CTProblems
using CTBase # for plot
using Plots;

# goddard with state constraint - maximize altitude
prob = Problem(:goddard, :state_constraint)
ocp = prob.model

# explicitely add a state constraint for plot
#println(constraints(ocp))
remove_constraint!(ocp,:state_constraint_v)
constraint!(ocp, :state, x->x[2], 0, 0.1, :state_con_v)
#println(constraints(ocp))

# initial guess (constant state and control functions)
init = [1.01, 0.05, 0.8, 0.1]

# solve problem
sol = solve(ocp, grid_size=100, print_level=5, tol=1e-12, mu_strategy="adaptive", init=init)

# plot
plot(sol)

# additional plots for constraint and multipliers (use info field from sol)
# +++todo: retrieve individual labels for constraints
t0 = sol.times[1]
tf = last(sol.times)
ncx = sol.infos[:dim_state_constraints]
PCX = Array{Plots.Plot, 1}(undef, ncx);
for i in 1:ncx
    PCX[i] = plot(t -> sol.infos[:state_constraints](t)[i], t0, tf, label="constraint", legend=:topleft)
    PCX[i] = plot!(twinx(), t -> sol.infos[:mult_state_constraints](t)[i], t0, tf, color=:red, label="multiplier", xticks=:none) #2nd scale
    #PCX[i] = plot!(t -> sol.infos[:mult_state_constraints](t)[i], t0, tf, color=:red, label="multiplier")
end
plot(PCX..., layout = (ncx,1));

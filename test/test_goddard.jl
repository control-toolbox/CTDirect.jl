# goddard with state constraint - maximize altitude  +++ use newer one in CTProblems
println("Goddard test")
prob = Problem(:goddard, :state_constraint)
ocp = prob.model

# initial guess (constant state and control functions)
init = [1.01, 0.05, 0.8, 0.1]

# all constraints formulation
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
sol = solve(ocp, grid_size=10, print_level=0, init=init)
@testset verbose = true showtiming = true ":goddard :all_constraints" begin
    @test sol.objective ≈ prob.solution.objective atol=5e-3
end

# +++ add automatic tests for multipliers, ie sign and complementarity with constraints ?
t = sol.times
# state constraint v <= 0.1
a = sol.infos[:mult_state_constraints](t)[1]
println(a)
println(all(>=(0), a) )
b = sol.infos[:state_constraints](t)[1]
println(b)
println(a.*b)
println(norm(a.*b))

# control constraint u <= 1

# mixed constraint m >= 0.6

# state box r >= 1 and v >= 0
# control box u >= 0 

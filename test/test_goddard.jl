# goddard with state constraint - maximize altitude
println("Goddard test")
prob = Problem(:goddard, :state_constraint)
ocp = prob.model

# initial guess (constant state and control functions)
init = [1.01, 0.05, 0.8, 0.1]

# state constraint formulation
println("State constraint formulation")
#println(constraints(ocp))
sol = solve(ocp, grid_size=10, print_level=0, init=init)
@testset verbose = true showtiming = true ":goddard :state_constraint" begin
    @test sol.objective ≈ prob.solution.objective atol=5e-3
end

#= # mixed constraint formulation
println("Mixed constraint formulation")
#remove_constraint!(ocp,+++)
constraint!(ocp, :mixed, (x,u)->x[2], 0, vmax, :mixed_con1)
sol = solve(ocp, grid_size=10, print_level=0, init=init)
@testset verbose = true showtiming = true ":goddard :mixed_constraint" begin
    @test sol.objective ≈ prob.solution.objective atol=5e-3
end =#

# box constraint formulation
println("Box constraint formulation")
#println(constraints(ocp))
remove_constraint!(ocp, :state_constraint_r)
remove_constraint!(ocp, :state_constraint_v)
remove_constraint!(ocp, :control_constraint)
constraint!(ocp, :state, 1:3, [1,0,0.6], [Inf,0.1,1], :box_state)
constraint!(ocp, :control, Index(1), 0, 1, :box_control)
#println(constraints(ocp)) #display(constraints(ocp))
sol = solve(ocp, grid_size=10, print_level=0, init=init)
println(sol.objective)
@testset verbose = true showtiming = true ":goddard :box_constraint" begin
    @test sol.objective ≈ prob.solution.objective atol=5e-3
end

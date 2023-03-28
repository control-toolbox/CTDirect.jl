# goddard with state constraint - maximize altitude
println("Goddard test")
prob = Problem(:goddard, :state_constraint)
ocp = prob.model
println("Is solvable ? ", CTDirect.is_solvable(ocp))

init = [1.01, 0.25, 0.5, 0.4]

#= # state constraint formulation
println("State constraint formulation")
sol = solve(ocp, grid_size=10, print_level=0, init=init)
@testset verbose = true showtiming = true ":goddard :state_constraint" begin
    @test sol.objective ≈ prob.solution.objective atol=1e-1
end =#

#= # mixed constraint formulation
println("Mixed constraint formulation")
#remove_constraint!(ocp,+++)
constraint!(ocp, :mixed, (x,u)->x[2], 0, vmax, :mixed_con1)
sol = solve(ocp, grid_size=10, print_level=0, init=init)
@testset verbose = true showtiming = true ":goddard :mixed_constraint" begin
    @test sol.objective ≈ prob.solution.objective atol=1e-1
end =#

# box constraint formulation
println("Box constraint formulation")
println(ocp.constraints)
remove_constraint!(ocp, :state_constraint_1)
remove_constraint!(ocp, :state_constraint_2)
remove_constraint!(ocp, :state_constraint_3)
remove_constraint!(ocp, :control_constraint_1)
constraint!(ocp, :state, 1:3, [1,0,0.6], [Inf,0.1,1], :box_state)
constraint!(ocp, :control, Index(1), 0, 1, :box_control)
println(ocp.constraints)
sol = solve(ocp, grid_size=10, print_level=0, init=init)
@testset verbose = true showtiming = true ":goddard :box_constraint" begin
    @test sol.objective ≈ prob.solution.objective atol=1e-2
end


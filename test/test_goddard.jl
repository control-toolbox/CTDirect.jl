# goddard with state constraint - maximize altitude
println("Goddard test")
prob = Problem(:goddard, :state_constraint)
ocp = prob.model
println("Is solvable ? ", CTDirect.is_solvable(ocp))

init = [1.01, 0.25, 0.5, 0.4]

# state constraint formulation
println("State constraint formulation")
sol = solve(ocp, grid_size=10, print_level=0, init=init)
@testset verbose = true showtiming = true ":goddard :state_constraint" begin
    @test sol.objective ≈ -1.0 atol=1e-1
end

# box constraint formulation
println("Box constraint formulation")
remove_constraint!(ocp, :state_con2)
vmax = 0.1
constraint!(ocp, :state, Index(2), 0, vmax, :box_con1)
sol = solve(ocp, grid_size=10, print_level=0, init=init)
@testset verbose = true showtiming = true ":goddard :box_constraint" begin
    @test sol.objective ≈ -1.0 atol=1e-1
end

# mixed constraint formulation
println("Mixed constraint formulation")
remove_constraint!(ocp, :box_con1)
constraint!(ocp, :mixed, (x,u)->x[2], 0, vmax, :mixed_con1)
sol = solve(ocp, grid_size=10, print_level=0, init=init)
@testset verbose = true showtiming = true ":goddard :mixed_constraint" begin
    @test sol.objective ≈ -1.0 atol=1e-1
end



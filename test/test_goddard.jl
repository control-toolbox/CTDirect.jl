# goddard with state constraint - maximize altitude
println("Goddard test")
prob = Problem(:goddard, :all_constraint)
init = [1.01, 0.05, 0.8, 0.1]
sol = solve(prob.model, grid_size=10, print_level=0, init=init)
@testset verbose = true showtiming = true ":goddard :all_constraints" begin
    @test sol.objective ≈ prob.solution.objective rtol=1e-2
end

# +++ add automatic tests for multipliers, ie sign and complementarity with constraints ?
#=
t = sol.times
# state constraint v <= 0.1
a = sol.infos[:mult_state_constraints](t)[1]
println(a)
println(all(>=(0), a) )
b = sol.infos[:state_constraints](t)[1]
println(b)
println(a.*b)
println(norm(a.*b))
=#

# control constraint u <= 1

# mixed constraint m >= 0.6

# state box r >= 1 and v >= 0
# control box u >= 0 

# manual test, including ie plots
using Plots
Plots.default(show = true)
include("test_common.jl")

include("./problems/beam.jl")
include("./problems/goddard.jl")
include("./problems/double_integrator.jl")
include("./problems/truck_trailer.jl")

#sol = solve_problem(truck_trailer(); display=true)
#println("J ", objective(sol), " tf ", variable(sol), " iter ", iterations(sol))
sol1 = solve_problem(truck_trailer(); discretizer=:direct_shooting, display=true)
println("J ", objective(sol1), " tf ", variable(sol1), " iter ", iterations(sol1))
sol2 = solve_problem(truck_trailer(); discretizer=:direct_shooting, grid_size=50, control_steps=10, display=true)
println("J ", objective(sol2), " tf ", variable(sol2))
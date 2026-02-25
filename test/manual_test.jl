# manual test, including ie plots
using Plots
Plots.default(show = true)
include("test_common.jl")

include("./problems/beam.jl")
include("./problems/goddard.jl")
include("./problems/double_integrator.jl")
include("./problems/truck_trailer.jl")

sol = solve_problem(truck_trailer(); display=true)
#sol1 = solve_problem(double_integrator_minenergy2(); display=true, modeler=:exa, solver=:madnlp)


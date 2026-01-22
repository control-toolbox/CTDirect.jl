# manual test, including ie plots
using Plots
Plots.default(show = true)
include("test_common.jl")

include("./problems/beam.jl")
sol = solve_problem(beam(); display=true)
#plot(sol)

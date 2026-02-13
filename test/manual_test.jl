# manual test, including ie plots
using Plots
Plots.default(show = true)
include("test_common.jl")

include("./problems/beam.jl")
include("./problems/goddard.jl")
include("./problems/double_integrator.jl")

sol = solve_problem(goddard(); display=true)
sol1 = solve_problem(goddard2(); display=true, modeler=:exa, solver=:madnlp)


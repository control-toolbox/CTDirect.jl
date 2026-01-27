# manual test, including ie plots
using Plots
Plots.default(show = true)
include("test_common.jl")

include("./problems/beam.jl")
include("./problems/goddard.jl")
sol = solve_problem(beam2(); display=true, solver=:madnlp, modeler=:exa)

#plot(sol) even basic plot is broken for all julia versions -_- 

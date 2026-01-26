# manual test, including ie plots
using Plots
Plots.default(show = true)
include("test_common.jl")

include("./problems/goddard.jl")
sol = solve_problem(goddard(); display=true)
#plot(sol) even basic plot is broken for all julia versions -_- 

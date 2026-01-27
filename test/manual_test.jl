# manual test, including ie plots
using Plots
Plots.default(show = true)
include("test_common.jl")

include("./problems/beam.jl")
include("./problems/goddard.jl")


sol = solve_problem(goddard2(); display=true, solver=:ipopt, modeler=:adnlp)
sol1 = solve_problem(goddard2(); display=true, solver=:ipopt, modeler=:adnlp, init=sol)

#plot(sol) even basic plot is broken for all julia versions -_- 

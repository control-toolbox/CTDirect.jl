# manual test, including ie plots
using Plots
Plots.default(show = true)
include("test_common.jl")

include("./problems/beam.jl")
include("./problems/goddard.jl")
include("./problems/double_integrator.jl")

prob = goddard()

sol = solve_problem(prob; display=true)
sol = solve_problem(prob; display=true, adnlp_backend=:manual)
sol = solve_problem(prob; display=true)

#plot(sol) even basic plot is broken for all julia versions -_- 

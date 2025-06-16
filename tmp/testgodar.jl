using Pkg
Pkg.activate("tmp")

using CTDirect: CTDirect, solve, direct_transcription, set_initial_guess, build_OCP_solution
using CTModels: CTModels, objective, state, control, variable, costate, time_grid, iterations
using CTParser: CTParser, @def, prefix!
prefix!(:CTModels) # tell CTParser def macro to use CTModels instead of OptimalControl

using NLPModelsIpopt
using Plots


include("../test/problems/goddard.jl")
problem = goddard()


sol = solve(problem[:ocp], grid_size=500, disc_method= :trapeze, display = false)
#@btime direct_sol = solve(ocp; grid_size=500, disc_method = :trapeze, display = false)
display(plot(sol))

solmid = solve(problem[:ocp], grid_size=500, disc_method= :midpoint, display = false)
display(plot(solmid))



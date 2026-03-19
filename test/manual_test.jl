# manual test, including ie plots
using Plots
Plots.default(show = true)
include("test_common.jl")

include("./problems/beam.jl")
include("./problems/goddard.jl")
include("./problems/double_integrator.jl")
include("./problems/truck_trailer.jl")

#=  Uno
=#
sol = solve_problem(beam(); solver=:uno, display=true)

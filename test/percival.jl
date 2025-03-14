using CTDirect
using Plots
using ADNLPModels
using Percival

include("./problems/goddard.jl")
include("./problems/double_integrator.jl")
include("./problems/simple_integrator.jl")
prob = simple_integrator() #double_integrator_minenergy() #goddard_all()
docp, nlp = direct_transcription(prob.ocp, init=prob.init, adnlp_backend=:default)#, matrix_free=true)
ADNLPModels.get_adbackend(nlp)
dsol = percival(nlp, verbose = 5)
#sol = OptimalControlSolution(docp, dsol)
#plot(sol)

# test is ok
#nlp = ADNLPModel(
#  x -> (x[1] - 1)^2 + 100 * (x[2] - x[1]^2)^2,
#  [-1.2; 1.0],
#  x -> [x[1]^2 + x[2]^2],
#  [1.0],
#  [1.0],
#)
#output = percival(nlp, verbose = 1)
using CTDirect
using CTBase
using Plots
using ADNLPModels
using MadNLP
using NLPModelsIpopt
using Percival
using BenchmarkTools

#= basic test is ok
nlp = ADNLPModel(
  x -> (x[1] - 1)^2 + 100 * (x[2] - x[1]^2)^2,
  [-1.2; 1.0],
  x -> [x[1]^2 + x[2]^2],
  [1.0],
  [1.0],
)
#output = percival(nlp, verbose = 1)=#

include("./problems/double_integrator.jl")
prob = double_integrator_minenergy()
# All solvers perform worse with adnlp_backend=:default that uses ForwardDiff only
# Percival performs the same with matrix_free=true that disables Jacobian / Hessian matrices
docp, nlp = direct_transcription(prob.ocp, grid_size=25)
println(ADNLPModels.get_adbackend(nlp))
println("ipopt")
@btime ipopt(nlp, print_level=0)
println("madnlp")
@btime madnlp(nlp, print_level=MadNLP.ERROR)
println("percival")
# @btime dsol = percival(nlp)
# same allocs/time as calling percival(nlp)
solver = PercivalSolver(nlp)
@btime solve!(solver, nlp)

# check solution
dsol = percival(nlp, verbose=1)
sol = OptimalControlSolution(docp, dsol)
display(plot(sol))

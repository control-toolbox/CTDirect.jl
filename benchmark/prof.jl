# Profiling
using CTDirect
using CTBase
using NLPModelsIpopt
using HSL
using Printf

#using BenchmarkTools
#using Traceur
#using Profile
#using PProf
using JET

precompile = true
test_time = true
test_code_warntype =false
test_jet = false
test_ipopt = false

# define OCP
prob = include("../problems/swimmer.jl")
ocp = prob[:ocp]
println("Load problem ", prob[:name])

# compilation
if precompile
  println("Precompilation")
  sol = solve(ocp, grid_size=50, display=false, max_iter=2)
end

# full solve
if test_time
  println("Timed solve")
  @timev sol = solve(ocp, grid_size=50, print_level=0)
end

if test_code_warntype
  println("@code_warntype objective")
  @code_warntype CTDirect.DOCP_objective(x0, docp)
end

if test_jet
  println("@report_opt objective")
  @report_opt CTDirect.DOCP_objective(x0, docp)
  #═════ 48 possible errors found ═════
end

# ipopt statistics
if test_ipopt
  docp = directTranscription(ocp, grid_size=200, init=(state=[1,0.2,0.5], control=0.5))
  dsol = solve(docp, print_level=5, tol=1e-12, print_timing_statistics="yes")
  println("\n Discrete solution")
  println(dsol)
end
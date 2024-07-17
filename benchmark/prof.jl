# Profiling
include("../test/deps.jl")

#using BenchmarkTools
#using Traceur
#using Profile
#using PProf
using JET

objective = :lagrange
precompile = true
test_time = false
test_code_warntype = true
test_jet = false
test_ipopt = false

# define OCP
if objective == :lagrange
  prob = include("../problems/fuller.jl")
  ocp = prob[:ocp]
  docp = direct_transcription(ocp)
  println("Load problem (lagrange): ", prob[:name])
else
  prob = include("../problems/jackson.jl")
  ocp = prob[:ocp]
  docp = direct_transcription(ocp)
  println("Load problem (mayer): ", prob[:name])
end

# precompilation
if precompile
  println("Precompilation")
  solve(ocp, grid_size=50, display=false, max_iter=2)
end

# full solve
if test_time
  println("Timed solve")
  @timev sol = solve(ocp, grid_size=50, print_level=0)
end


if test_code_warntype
  println("@code_warntype objective")
  # x0, xf, obj: Any ...
  @code_warntype CTDirect.DOCP_objective(CTDirect.DOCP_initial_guess(docp), docp)
  # OK !
  @code_warntype CTDirect.DOCP_constraints!(zeros(docp.dim_NLP_constraints), CTDirect.DOCP_initial_guess(docp), docp)
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
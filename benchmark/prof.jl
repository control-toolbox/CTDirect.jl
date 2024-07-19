# Profiling
include("../test/deps.jl")

#using BenchmarkTools
#using Traceur
#using Profile
#using PProf
using JET

precompile = true

test_time = true
#test = :objective
test = :constraints
test_code_warntype = false
test_jet = false

# define OCP
prob = include("../problems/fuller.jl")
#prob = include("../problems/jackson.jl")
#prob = include("../problems/goddard.jl")
ocp = prob[:ocp]
grid_size = 200
docp, nlp = direct_transcription(ocp, grid_size=grid_size)
println("Load problem ", prob[:name])


# full solve
if test_time
  if precompile
    println("Precompilation")
    solve(ocp, grid_size=grid_size, display=false, max_iter=2)
  end
  println("Timed solve")
  @timev sol = solve(ocp, grid_size=grid_size, print_level=0)
end

if precompile
  println("Precompilation")
  if test == :objective
    CTDirect.DOCP_objective(CTDirect.DOCP_initial_guess(docp), docp)
  else
    CTDirect.DOCP_constraints!(zeros(docp.dim_NLP_constraints), CTDirect.DOCP_initial_guess(docp), docp)
  end
end

if test_code_warntype
  if test == :objective
    # NB. Pb with the mayer part: obj is type unstable (Any) because ocp.mayer is Union(Mayer,nothing), even for mayer problems (also, we should not even enter this code part for lagrange problems since has_mayer us defined as const in DOCP oO ...).
    @code_warntype CTDirect.DOCP_objective(CTDirect.DOCP_initial_guess(docp), docp)
  else
    # OK ! only index is Any but typing it seems worse...
    @code_warntype CTDirect.DOCP_constraints!(zeros(docp.dim_NLP_constraints), CTDirect.DOCP_initial_guess(docp), docp)
  end
end

if test_jet
  if test == :objective
    # 4 possible errors
    # due to the ocp.mayer type problem cf above
    @report_opt CTDirect.DOCP_objective(CTDirect.DOCP_initial_guess(docp), docp)
  else
    # 46 possible errors: some getindex (Integer vs Int...)
    # ArgsAtTimeStep, setPathConstraintsAtTimeStep...
    @report_opt CTDirect.DOCP_constraints!(zeros(docp.dim_NLP_constraints), CTDirect.DOCP_initial_guess(docp), docp)
  end
end

#=
# ipopt statistics
if test_ipopt
  docp = directTranscription(ocp, grid_size=200, init=(state=[1,0.2,0.5], control=0.5))
  dsol = solve(docp, print_level=5, tol=1e-12, print_timing_statistics="yes")
  println("\n Discrete solution")
  println(dsol)
end
=#
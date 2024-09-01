# Profiling
include("../test/deps.jl")

using BenchmarkTools
#using Traceur
#using Profile
#using PProf
using JET

precompile = true

test_solve = true
#test = :objective
test = :constraints
test_code_warntype = false
test_jet = false

# define OCP
include("../test/problems/goddard.jl")
prob = goddard_all()
ocp = prob[:ocp]
grid_size = 100
docp, nlp = direct_transcription(ocp, grid_size = grid_size)
println("Load problem ", prob[:name])

if precompile
    println("Precompilation")
    direct_solve(ocp, grid_size = grid_size, display = false, max_iter = 2)
    CTDirect.DOCP_objective(CTDirect.DOCP_initial_guess(docp), docp)
    CTDirect.DOCP_constraints!(
        zeros(docp.dim_NLP_constraints),
        CTDirect.DOCP_initial_guess(docp),
        docp,
    )
end

# full solve
if test_solve
    println("Timed solve")
    @timev sol = direct_solve(ocp, grid_size = grid_size, display = false)
    @btime sol = direct_solve(ocp, grid_size = grid_size, display = false)
end

if test_code_warntype
    if test == :objective
        # NB. Pb with the mayer part: obj is type unstable (Any) because ocp.mayer is Union(Mayer,nothing), even for mayer problems (also, we should not even enter this code part for lagrange problems since has_mayer us defined as const in DOCP oO ...).
        @code_warntype CTDirect.DOCP_objective(CTDirect.DOCP_initial_guess(docp), docp)
    else
        # OK !
        @code_warntype CTDirect.DOCP_constraints!(
            zeros(docp.dim_NLP_constraints),
            CTDirect.DOCP_initial_guess(docp),
            docp,
        )
    end
end

if test_jet
    if test == :objective
        # 4 possible errors
        # due to the ocp.mayer type problem cf above
        @report_opt CTDirect.DOCP_objective(CTDirect.DOCP_initial_guess(docp), docp)
    else
        # 50 possible errors: some getindex (Integer vs Int...)
        # all variables x,u,v
        @report_opt CTDirect.DOCP_constraints!(
            zeros(docp.dim_NLP_constraints),
            CTDirect.DOCP_initial_guess(docp),
            docp,
        )
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

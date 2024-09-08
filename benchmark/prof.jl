# +++ TODO: make function with bools as args ?
# Profiling
using CTDirect
using CTBase

using LinearAlgebra
using NLPModelsIpopt
using BenchmarkTools
#using Traceur
#using Profile
#using PProf
using JET

precompile = true
test_objective = true
test_constraints = true
test_transcription = true
test_solve = true

test_code_warntype = false
test_jet = false

# define OCP
include("../test/problems/goddard.jl")
prob = goddard_all()
ocp = prob[:ocp]
grid_size = 100
println("Load problem ", prob[:name])

if precompile
    println("Precompilation")
    if test_transcription
        docp, nlp = direct_transcription(ocp, grid_size = grid_size)
    end
    if test_solve
        direct_solve(ocp, grid_size = grid_size, display = false, max_iter = 2)
    end
    if test_objective
        CTDirect.DOCP_objective(CTDirect.DOCP_initial_guess(docp), docp)
    end
    if test_constraints
        CTDirect.DOCP_constraints!(
            zeros(docp.dim_NLP_constraints),
            CTDirect.DOCP_initial_guess(docp),
            docp,
        )
    end
end

# evaluation
if test_objective
    println("Timed objective")
    @btime CTDirect.DOCP_objective(CTDirect.DOCP_initial_guess(docp), docp)
end
if test_constraints
    println("Timed constraints")
    @btime CTDirect.DOCP_constraints!(
        zeros(docp.dim_NLP_constraints),
        CTDirect.DOCP_initial_guess(docp),
        docp,
    )
end

# transcription
if test_transcription
    println("Timed transcription")
    @btime docp, nlp = direct_transcription(ocp, grid_size = grid_size)
end

# full solve
if test_solve
    println("Timed full solve")
    @btime sol = direct_solve(ocp, grid_size = grid_size, display = false)
end

if test_code_warntype
    if test_objective
        # NB. Pb with the mayer part: obj is type unstable (Any) because ocp.mayer is Union(Mayer,nothing), even for mayer problems (also, we should not even enter this code part for lagrange problems since has_mayer us defined as const in DOCP oO ...).
        @code_warntype CTDirect.DOCP_objective(CTDirect.DOCP_initial_guess(docp), docp)
    end
    if test_constraints
        # OK !
        @code_warntype CTDirect.DOCP_constraints!(
            zeros(docp.dim_NLP_constraints),
            CTDirect.DOCP_initial_guess(docp),
            docp,
        )
    end
end

if test_jet
    if test_objective
        # 4 possible errors
        # due to the ocp.mayer type problem cf above
        @report_opt CTDirect.DOCP_objective(CTDirect.DOCP_initial_guess(docp), docp)
    end
    if test_constraints
        # 50 possible errors: some getindex (Integer vs Int...)
        # all variables x,u,v
        @report_opt CTDirect.DOCP_constraints!(
            zeros(docp.dim_NLP_constraints),
            CTDirect.DOCP_initial_guess(docp),
            docp,
        )
    end
end

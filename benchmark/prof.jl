# tests to check allocations in particular
using CTDirect
using CTBase

using LinearAlgebra
using NLPModelsIpopt
using BenchmarkTools
using Profile

include("../test/problems/goddard.jl")


function unit(test_obj=true)
    
    # define problem and variables
    prob = goddard_all()
    ocp = prob[:ocp]
    docp,_ = direct_transcription(ocp)
    xu = CTDirect.DOCP_initial_guess(docp)
    c = Vector{Float64}(undef, docp.dim_NLP_constraints)

    # DOCP_objective
    if test_obj
        CTDirect.DOCP_objective(xu, docp) # compile
        Profile.clear_malloc_data()
        a = @allocated begin
            obj = CTDirect.DOCP_objective(xu, docp)
        end
        println("DOCP_objective ", a)
    end

end

#=
# call to setPointConstraints
v = CTDirect.get_optim_variable(xu, docp)
a = @allocated begin
    CTDirect.setPointConstraints!($docp, c, $xu, $v)
end; a > 0 && @show a
# call to setPathConstraints 
# call to setConstraintsBlock
# call to DOCP_constraints (also check for remaining undef in c)
=#

# check vectorize too !
# tests to check allocations in particular
using CTDirect
using CTBase

using LinearAlgebra
using NLPModelsIpopt
using BenchmarkTools
using Profile

include("../test/problems/goddard.jl")


function test_unit(test_obj=false, test_cons=true, grid_size=10)
    
    # define problem and variables
    prob = goddard_all()
    ocp = prob[:ocp]
    docp,_ = direct_transcription(ocp, grid_size=grid_size) # nlp creates more allocs !
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

    # DOCP_constraints
    if test_cons
        CTDirect.DOCP_constraints!(c, xu, docp) # compile
        @btime CTDirect.DOCP_constraints!($c, $xu, $docp)
        Profile.clear_malloc_data()
        b = @allocated begin
            CTDirect.DOCP_constraints!(c, xu, docp)
        end
        println("DOCP_constraints! ", b)
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
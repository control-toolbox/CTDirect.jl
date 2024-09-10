# tests to check allocations in particular
using CTDirect
using CTBase

using LinearAlgebra
using NLPModelsIpopt
using BenchmarkTools
using Profile

include("../test/problems/goddard.jl")

function test_unit(test_obj=false, test_cons=false, test_dyn=true, grid_size=10)
    
    # define problem and variables
    prob = goddard_all()
    ocp = prob[:ocp]
    docp,_ = direct_transcription(ocp, grid_size=grid_size) # nlp creates more allocs !
    xu = CTDirect.DOCP_initial_guess(docp)
    #c = Vector{Float64}(undef, docp.dim_NLP_constraints)
    c = similar(xu, docp.dim_NLP_constraints)
    work = similar(xu, docp.dim_NLP_x)
    t = xu[end]
    v = CTDirect.get_optim_variable(xu, docp)
    x, u = CTDirect.get_variables_at_time_step(xu, docp, docp.dim_NLP_steps)

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
        #@btime CTDirect.DOCP_constraints!($c, $xu, $docp)
        Profile.clear_malloc_data()
        a = @allocated begin
            CTDirect.DOCP_constraints!(c, xu, docp)
        end
        println("DOCP_constraints! ", a)
    end

    # dynamics and work array
    if test_dyn
        docp.ocp.dynamics(t, x, u, v) # compile
        Profile.clear_malloc_data()
        a = @allocated begin
            #docp.ocp.dynamics(t, x, u, v) # 384 ie 6 float64 (args, or vect args + ret ?)
            #work[1:docp.dim_OCP_x] = docp.ocp.dynamics(t, x, u, v) # 416 +32
            #work[1:docp.dim_OCP_x] .= docp.ocp.dynamics(t, x, u, v) # 496 worse
            #@views work[1:docp.dim_OCP_x] = docp.ocp.dynamics(t, x, u, v) # 416 same
        end
        println("dynamics and work array ", a)
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
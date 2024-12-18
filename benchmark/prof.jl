# tests to check allocations in particular
using CTDirect
using CTBase

using LinearAlgebra
using NLPModelsIpopt
using BenchmarkTools
using Profile
using PProf
using JET

include("../test/problems/goddard.jl")
include("../test/problems/simple_integrator.jl")

# local version of mayer cost
function local_mayer(obj, x0, xf, v)
    obj[1] = xf[3]
    return
end

function init(prob ;grid_size, disc_method)
    ocp = prob[:ocp]
    docp = CTDirect.DOCP(ocp, grid_size=grid_size, time_grid=CTDirect.__time_grid(), disc_method=disc_method)
    xu = CTDirect.DOCP_initial_guess(docp)
    return prob, docp, xu
end


function test_unit(prob ;test_get=false, test_obj=true, test_cons=true, test_trans=true, test_solve=true, warntype=false, jet=false, profile=false, grid_size=250, disc_method=:trapeze)

    if profile
        Profile.Allocs.clear()
    end

    # define problem and variables
    prob, docp, xu = init(prob; grid_size=grid_size, disc_method=disc_method)
    disc = docp.discretization
    #= OK, same as calling the functions with docp
    NLP_objective = (xu) -> CTDirect.DOCP_objective(xu, docp)
    NLP_constraints! = (c, xu) -> CTDirect.DOCP_constraints!(c, xu, docp) =#
    c = fill(666.666, docp.dim_NLP_constraints)
    work = similar(xu, docp.dim_NLP_x)

    # getters
    if test_get
        println("Getters")
        print("t"); @btime CTDirect.get_final_time($xu, $docp)
        print("x"); @btime CTDirect.get_OCP_state_at_time_step($xu, $docp, 1)
        print("u"); @btime CTDirect.get_OCP_control_at_time_step($xu, $docp, 1)
        print("v"); @btime CTDirect.get_OCP_variable($xu, $docp)
        if warntype
            @code_warntype CTDirect.get_final_time(xu, docp)
            @code_warntype CTDirect.get_time_grid(xu, docp)
            @code_warntype CTDirect.get_OCP_state_at_time_step(xu, docp, 1)
            @code_warntype CTDirect.get_OCP_control_at_time_step(xu, docp, 1)
            @code_warntype CTDirect.get_OCP_variable(xu, docp)
        end
    end

    # DOCP_objective
    if test_obj
        print("Objective"); @btime CTDirect.DOCP_objective($xu, $docp)
        warntype && @code_warntype CTDirect.DOCP_objective(xu, docp)
        jet && display(@report_opt CTDirect.DOCP_objective(xu, docp))
        if profile
            Profile.Allocs.@profile sample_rate=1.0 CTDirect.DOCP_objective(xu, docp)
            results = Profile.Allocs.fetch()
            PProf.Allocs.pprof()
        end
    end

    # DOCP_constraints
    if test_cons
        print("Constraints"); @btime CTDirect.DOCP_constraints!($c, $xu, $docp)
        any(c.==666.666) && error("undefined values in constraints ",c)
        warntype && @code_warntype CTDirect.DOCP_constraints!(c, xu, docp)
        jet && display(@report_opt CTDirect.DOCP_constraints!(c, xu, docp))
        if profile
            Profile.Allocs.@profile sample_rate=1.0 CTDirect.DOCP_constraints!(c, xu, docp)
            results = Profile.Allocs.fetch()
            PProf.Allocs.pprof()
        end
    end

    # transcription
    if test_trans
        print("Transcription"); @btime direct_transcription($prob.ocp, grid_size=$grid_size, disc_method=$disc_method)
    end

    # solve
    if test_solve
        sol = direct_solve(prob.ocp, display=false, grid_size=grid_size, disc_method=disc_method)
        if !isapprox(sol.objective, prob.obj, rtol=1e-2)
            error("objective mismatch: ", sol.objective, " vs ", prob.obj)
        end
        print("Solve"); @btime direct_solve($prob.ocp, display=false, grid_size=$grid_size, disc_method=$disc_method)
    end

end


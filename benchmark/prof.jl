# tests to check allocations in particular
using CTDirect
import CTModels

using LinearAlgebra
using NLPModelsIpopt
using BenchmarkTools
using Profile
using PProf
using JET


function init(ocp; grid_size, disc_method)
    docp = CTDirect.DOCP(ocp, grid_size=grid_size, time_grid=CTDirect.__time_grid(), disc_method=disc_method)
    xu = CTDirect.DOCP_initial_guess(docp)
    return docp, xu
end


function test_unit(ocp; test_obj=true, test_cons=true, test_trans=true, test_solve=true, warntype=false, jet=false, profile=false, grid_size=100, disc_method=:trapeze)

    if profile
        Profile.Allocs.clear()
    end

    # define problem and variables
    docp, xu = init(ocp, grid_size=grid_size, disc_method=disc_method)
    disc = docp.discretization
    #= OK, same as calling the functions with docp
    NLP_objective = (xu) -> CTDirect.DOCP_objective(xu, docp)
    NLP_constraints! = (c, xu) -> CTDirect.DOCP_constraints!(c, xu, docp) =#
    c = fill(666.666, docp.dim_NLP_constraints)
    work = similar(xu, docp.dim_NLP_x)

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
        print("Transcription"); @btime direct_transcription($ocp, grid_size=$grid_size, disc_method=$disc_method)
    end

    # solve
    if test_solve
        sol = direct_solve(ocp, display=false, grid_size=grid_size, disc_method=disc_method)
        #if !isapprox(sol.objective, prob.obj, rtol=1e-2)
        #    error("objective mismatch: ", sol.objective, " vs ", prob.obj)
        #end
        print("Solve"); @btime direct_solve($ocp, display=false, grid_size=$grid_size, disc_method=$disc_method)
    end

end


using CTDirect
using NLPModelsIpopt
import CTModels

using BenchmarkTools
using JET
using Profile
using PProf

function test_model(ocp; warntype=false, jet=false, profile=false)

    if profile
        Profile.Allocs.clear()
    end

    # test objective and constraints
    docp,_ = direct_transcription(ocp)
    xu = CTDirect.DOCP_initial_guess(docp)
    c = fill(666.666, docp.dim_NLP_constraints)

    print("Objective"); @btime CTDirect.DOCP_objective($xu, $docp)
    warntype && @code_warntype CTDirect.DOCP_objective(xu, docp)
    jet && display(@report_opt CTDirect.DOCP_objective(xu, docp))
    if profile 
        Profile.Allocs.@profile sample_rate=1.0 CTDirect.DOCP_objective(xu, docp)
        results = Profile.Allocs.fetch()
        PProf.Allocs.pprof()
    end

    print("Constraints"); @btime CTDirect.DOCP_constraints!($c, $xu, $docp)
    any(c.==666.666) && error("undefined values in constraints ",c)
    warntype && @code_warntype CTDirect.DOCP_constraints!(c, xu, docp)
    jet && display(@report_opt CTDirect.DOCP_constraints!(c, xu, docp))
    if profile
        Profile.Allocs.@profile sample_rate=1.0 CTDirect.DOCP_constraints!(c, xu, docp)
        results = Profile.Allocs.fetch()
        PProf.Allocs.pprof()
    end

    #print("Transcription"); @btime direct_transcription($ocp)
    #warntype && @code_warntype direct_transcription(ocp)

    #print("Solve"); @btime direct_solve($ocp, display=false)

end
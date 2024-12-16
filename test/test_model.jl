using CTDirect
using NLPModelsIpopt
import CTModels

using BenchmarkTools
using JET
using Profile
using PProf


# define problem with new model: simple integrator
function simple_integrator_model()
    pre_ocp = CTModels.PreModel()
    CTModels.state!(pre_ocp, 1)
    CTModels.control!(pre_ocp, 2)
    CTModels.time!(pre_ocp, t0=0.0, tf=1.0)
    f!(r, t, x, u, v) = r .=  .- x[1] .- u[1] .+ u[2] 
    CTModels.dynamics!(pre_ocp, f!)
    l!(r, t, x, u, v) = r .= (u[1] .+ u[2]).^2
    CTModels.objective!(pre_ocp, :min, lagrange=l!)
    function bc!(r, x0, xf, v)
        r[1] = x0[1]
        r[2] = xf[1]
    end
    CTModels.constraint!(pre_ocp, :boundary, f=bc!, lb=[-1, 0], ub=[-1, 0], label=:boundary)
    CTModels.constraint!(pre_ocp, :control, rg=1:2, lb=[0, 0], ub=[Inf, Inf], label=:control_rg)
    CTModels.definition!(pre_ocp, Expr(:simple_integrator_min_energy))
    ocp = CTModels.build_model(pre_ocp)
    return ocp
end

# define problem with new model: double integrator
function double_integrator_mintf_model()
    pre_ocp = CTModels.PreModel()
    CTModels.state!(pre_ocp, 2)
    CTModels.control!(pre_ocp, 1)
    CTModels.variable!(pre_ocp, 1)
    CTModels.time!(pre_ocp, t0=0.0, indf=1)
    function f!(r, t, x, u, v)
        r[1] = x[2]
        r[2] = u[1]
    end 
    CTModels.dynamics!(pre_ocp, f!)
    function mayer!(r, x0, xf, v)
        r[1] = v[1]
    end 
    CTModels.objective!(pre_ocp, :min, mayer=mayer!)
    function bc!(r, x0, xf, v)
        r[1] = x0[1]
        r[2] = x0[2]
        r[3] = xf[1]
        r[4] = xf[2]
    end
    CTModels.constraint!(pre_ocp, :boundary, f=bc!, lb=[0, 0, 1, 0], ub=[0, 0, 1, 0], label=:boundary)
    CTModels.constraint!(pre_ocp, :control, rg=1:1, lb=[-1], ub=[1], label=:control_rg)
    CTModels.constraint!(pre_ocp, :variable, rg=1:1, lb=[0.05], ub=[Inf], label=:variable_rg)
    CTModels.definition!(pre_ocp, Expr(:double_integrator_min_tf))
    ocp = CTModels.build_model(pre_ocp)
    return ocp
end


# define problem with new model: fuller
function fuller_model()
    pre_ocp = CTModels.PreModel()
    CTModels.state!(pre_ocp, 2)
    CTModels.control!(pre_ocp, 1)
    CTModels.time!(pre_ocp, t0=0.0, tf=3.5)
    function f!(r, t, x, u, v)
        r[1] = x[2]
        r[2] = u[1]
    end 
    CTModels.dynamics!(pre_ocp, f!)
    function l!(r, t, x, u, v)
        r[1] = x[1]^2
    end
    CTModels.objective!(pre_ocp, :min, lagrange=l!)
    function bc!(r, x0, xf, v)
        r[1] = x0[1]
        r[2] = x0[2]
        r[3] = xf[1]
        r[4] = xf[2]
    end
    CTModels.constraint!(pre_ocp, :boundary, f=bc!, lb=[0, 1, 0, 0], ub=[0, 1, 0, 0], label=:boundary)
    CTModels.constraint!(pre_ocp, :control, rg=1:1, lb=[-1], ub=[1], label=:control_rg)
    CTModels.definition!(pre_ocp, Expr(:fuller_min_energy))
    ocp = CTModels.build_model(pre_ocp)
    return ocp
end


function test_model(ocp; warntype=false, jet=false, profile=false)

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

    print("Transcription"); @btime direct_transcription($ocp)
    # more allocs here, try to type the nlp constraints functions in docp, check problem size vs main ?

    print("Solve"); @btime direct_solve($ocp, display=false)

end
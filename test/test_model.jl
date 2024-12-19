using CTDirect
using NLPModelsIpopt
import CTModels

using BenchmarkTools
using JET
using Profile
using PProf

# beam

# jackson

# robbins

# vanderpol

# define problem with new model: simple integrator
function simple_integrator_model()
    pre_ocp = CTModels.PreModel()
    CTModels.state!(pre_ocp, 1)
    CTModels.control!(pre_ocp, 2)
    CTModels.time!(pre_ocp, t0=0.0, tf=1.0)
    f!(r, t, x, u, v) = r .=  .- x[1] .- u[1] .+ u[2] 
    CTModels.dynamics!(pre_ocp, f!)
    l(t, x, u, v) = (u[1] .+ u[2]).^2
    CTModels.objective!(pre_ocp, :min, lagrange=l)
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
    mayer(x0, xf, v) = v[1]
    CTModels.objective!(pre_ocp, :min, mayer=mayer)
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
    l(t, x, u, v) = x[1]^2
    CTModels.objective!(pre_ocp, :min, lagrange=l)
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

function algal_bacterial_model()

    s_in = 0.5
    β = 23e-3
    γ = 0.44
    dmax = 1.5
    ϕmax = 6.48; ks = 0.09;
    ρmax = 27.3e-3; kv = 0.57e-3;
    μmax = 1.0211; qmin = 2.7628e-3;
    x0 = [0.1629, 0.0487, 0.0003, 0.0177, 0.035, 0]
    ϕ(s) = ϕmax * s / (ks + s)
    ρ(v) = ρmax * v / (kv + v)
    μ(q) = μmax * (1 - qmin / q)

    pre_ocp = CTModels.PreModel()
    CTModels.state!(pre_ocp, 6)
    CTModels.control!(pre_ocp, 2)
    CTModels.time!(pre_ocp, t0=0.0, tf=20.0)
    function bc!(r, x0, xf, v)
        r .= x0
    end
    CTModels.constraint!(pre_ocp, :boundary, f=bc!, lb=x0, ub=x0, label=:boundary)
    CTModels.constraint!(pre_ocp, :state, rg=1:6, lb=[0, 0, 0, qmin, 0, 0], ub=[Inf,Inf,Inf,Inf,Inf,Inf], label=:state_rg)
    CTModels.constraint!(pre_ocp, :control, rg=1:2, lb=[0, 0], ub=[1, dmax], label=:control_rg)
    function f!(r, t, x, u, v)
        r[1] = u[2]*(s_in - x[1]) - ϕ(x[1])*x[2]/γ
        r[2] = ((1 - u[1])*ϕ(x[1]) - u[2])*x[2]
        r[3] = u[1]*β*ϕ(x[1])*x[2] - ρ(x[3])*x[5] - u[2]*x[3]
        r[4] = ρ(x[3]) - μ(x[4])*x[4]
        r[5] = (μ(x[4]) - u[2])*x[5]
        r[6] = u[2]*x[5]
    end
    CTModels.dynamics!(pre_ocp, f!)
    mayer(x0, xf, v) = xf[6]
    CTModels.objective!(pre_ocp, :max, mayer=mayer)
    CTModels.definition!(pre_ocp, Expr(:algal_bacterial))
    ocp = CTModels.build_model(pre_ocp)
    return ocp
end

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
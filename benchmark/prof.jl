# tests to check allocations in particular
using CTDirect
using CTBase

using LinearAlgebra
using NLPModelsIpopt
using BenchmarkTools
using Profile
using PProf
using JET

#include("../test/problems/goddard.jl")
include("../test/problems/simple_integrator.jl")

# local version of mayer cost
function local_mayer(obj, x0, xf, v)
    obj[1] = xf[3]
    return
end

function init(; grid_size, disc_method)
    #prob = goddard_all()
    #prob = goddard()
    prob = simple_integrator()
    ocp = prob[:ocp]
    docp = CTDirect.DOCP(ocp, grid_size=grid_size, time_grid=CTDirect.__time_grid(), disc_method=string(disc_method))
    xu = CTDirect.DOCP_initial_guess(docp)
    return prob, docp, xu
end


function test_unit(; test_get=false, test_dyn=false, test_unit_cons=false, test_mayer=false, test_obj=true, test_block=false, test_cons=true, test_trans=true, test_solve=true, warntype=false, jet=false, profile=false, grid_size=100, disc_method=:trapeze)

    # define problem and variables
    prob, docp, xu = init(grid_size=grid_size, disc_method=disc_method)
    disc = docp.discretization
    #= OK, same as calling the functions with docp
    NLP_objective = (xu) -> CTDirect.DOCP_objective(xu, docp)
    NLP_constraints! = (c, xu) -> CTDirect.DOCP_constraints!(c, xu, docp) =#
    c = fill(666.666, docp.dim_NLP_constraints)
    work = similar(xu, docp.dim_NLP_x)

    # getters
    if test_get
        println("Getters")
        print("t")
        @btime CTDirect.get_final_time($xu, $docp)
        print("x")
        @btime CTDirect.get_OCP_state_at_time_step($xu, $docp, 1)
        print("u")
        @btime CTDirect.get_OCP_control_at_time_step($xu, $docp, 1)
        print("v")
        @btime CTDirect.get_OCP_variable($xu, $docp)
        if warntype
            @code_warntype CTDirect.get_final_time(xu, docp)
            @code_warntype CTDirect.get_time_grid(xu, docp)
            @code_warntype CTDirect.get_OCP_state_at_time_step(xu, docp, 1)
            @code_warntype CTDirect.get_OCP_control_at_time_step(xu, docp, 1)
            @code_warntype CTDirect.get_OCP_variable(xu, docp)
        end
    end

    f = similar(xu, docp.dim_NLP_x)

    if test_dyn
        print("dynamics_ext")
        @btime $docp.dynamics_ext($f, $t, $x, $u, $v)
        warntype && @code_warntype docp.dynamics_ext(f, t, x, u, v)
        Profile.clear_malloc_data()
        docp.dynamics_ext(f, t, x, u, v)
    end

    if test_unit_cons
        print("u cons")
        @btime $docp.control_constraints[2]($c, $t, $u, $v)
        print("x cons")
        @btime $docp.state_constraints[2]($c, $t, $x, $v)
        print("xu cons")
        @btime $docp.mixed_constraints[2]($c, $t, $x, $u, $v)
        if warntype
            @code_warntype docp.control_constraints[2](c, t, u, v)
            @code_warntype docp.state_constraints[2](c, t, x, v)
            @code_warntype docp.mixed_constraints[2](c, t, x, u, v)
        end
    end

    # objective
    if test_mayer
        n = docp.dim_OCP_x
        nx = docp.dim_NLP_x
        m = docp.dim_NLP_u
        N = docp.dim_NLP_steps
        x0 = CTDirect.get_OCP_state_at_time_step(xu, docp, 1)
        xf = CTDirect.get_OCP_state_at_time_step(xu, docp, N + 1)
        v = CTDirect.get_OCP_variable(xu, docp)
        obj = similar(xu, 1)

        # local mayer
        println("")
        print("Local Mayer: views for x0/xf and scalar v")
        @btime local_mayer($obj, (@view $xu[1:$n]), (@view $xu[($nx+$m)*$N+1:($nx+$m)*$N+$n]), $xu[end]) # OK
        print("Local Mayer: param scal/vec getters")
        @btime local_mayer($obj, $x0, $xf, $v) # OK
        print("OCP Mayer: param scal/vec getters")
        @btime $docp.ocp.mayer($obj, $x0, $xf, $v) # 3 allocs (112)

        warntype && @code_warntype docp.ocp.mayer(obj, x0, xf, v)
        jet && display(@report_opt docp.ocp.mayer(obj, x0, xf, v))
        if profile
            Profile.Allocs.@profile sample_rate = 1.0 docp.ocp.mayer(obj, x0, xf, v)
            results = Profile.Allocs.fetch()
            PProf.Allocs.pprof()
        end
    end

    if test_obj
        print("Objective")
        @btime CTDirect.DOCP_objective($xu, $docp)
        warntype && @code_warntype CTDirect.DOCP_objective(xu, docp)
        jet && display(@report_opt CTDirect.DOCP_objective(xu, docp))
        if profile
            Profile.Allocs.@profile sample_rate = 1.0 CTDirect.DOCP_objective(xu, docp)
            results = Profile.Allocs.fetch()
            PProf.Allocs.pprof()
        end
    end

    if test_block
        times = CTDirect.get_time_grid(xu, docp) # type OK
        i = 1
        v = CTDirect.get_OCP_variable(xu, docp)
        work = CTDirect.setWorkArray(docp, xu, times, v)
        print("Constraints block")
        @btime CTDirect.setConstraintBlock!($docp, $c, $xu, $v, $times, $i, $work)
        warntype && @code_warntype CTDirect.setConstraintBlock!(docp, c, xu, v, times, i, work)
        jet && display(@report_opt CTDirect.setConstraintBlock!(docp, c, xu, v, times, i, work))
        if profile
            Profile.Allocs.@profile sample_rate = 1.0 CTDirect.setConstraintBlock!(docp, c, xu, v, times, i, work)
            results = Profile.Allocs.fetch()
            PProf.Allocs.pprof()
        end
    end

    # DOCP_constraints
    if test_cons
        print("Constraints")
        @btime CTDirect.DOCP_constraints!($c, $xu, $docp)
        any(c .== 666.666) && error("undefined values in constraints ", c)
        warntype && @code_warntype CTDirect.DOCP_constraints!(c, xu, docp)
        jet && display(@report_opt CTDirect.DOCP_constraints!(c, xu, docp))
        if profile
            Profile.Allocs.@profile sample_rate = 1.0 CTDirect.DOCP_constraints!(c, xu, docp)
            results = Profile.Allocs.fetch()
            PProf.Allocs.pprof()
        end
    end

    # transcription
    if test_trans
        print("Transcription")
        @btime direct_transcription($prob.ocp, grid_size=$grid_size)
    end

    # solve
    if test_solve
        sol = direct_solve(prob.ocp, display=false, grid_size=grid_size)
        if !isapprox(sol.objective, prob.obj, rtol=1e-2)
            error("objective mismatch: ", sol.objective, " vs ", prob.obj)
        end
        print("Solve")
        @btime direct_solve($prob.ocp, display=false, grid_size=$grid_size)
    end

end

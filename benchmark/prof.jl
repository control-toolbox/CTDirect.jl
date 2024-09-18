# tests to check allocations in particular
using CTDirect
using CTBase

using LinearAlgebra
using NLPModelsIpopt
using BenchmarkTools
using Profile
using PProf

include("../test/problems/goddard.jl")
#include("../test/problems/double_integrator.jl")

# local version of dynamics
Cd = 310
beta = 500
b = 2
Tmax = 3.5
# compact function is not better...
function compact_dynamics(t, x, u, vv)
    r, v, m = x
    D = Cd * v^2 * exp(-beta * (r - 1))
    return [v,
            -D / m - 1 / r^2 + u * Tmax / m ,
            - b * u * Tmax]
end
function F0(x, Cd, beta)
    r, v, m = x
    D = Cd * v^2 * exp(-beta * (r - 1))
    return [v, -D / m - 1 / r^2, 0]
end
function F1(x, Tmax, b)
    r, v, m = x
    return [0, Tmax / m, -b * Tmax]
end
function local_dynamics(t, x, u, v)
    return F0(x, Cd, beta) + u * F1(x, Tmax, b)
end


function dummy_dynamics(t, x, u, vv)
    return [x[2], x[1], u]
end

function init(;in_place, grid_size, discretization)
    if in_place
        prob = goddard_all_inplace()
        #prob = double_integrator_a()
    else
        prob = goddard_all()
        #prob = double_integrator_mintf()
    end
    ocp = prob[:ocp]
    discretization = string(discretization)
    if discretization == "midpoint"
        disc_method = CTDirect.Midpoint()
    elseif discretization == "trapeze"
        disc_method = CTDirect.Trapeze()
    else
        error("Unknown discretization method:", discretization)
    end
    docp = CTDirect.DOCP(ocp, grid_size=grid_size, time_grid=CTDirect.__time_grid(), discretization=disc_method)
    xu = CTDirect.DOCP_initial_guess(docp)
    return prob, docp, xu
end

function test_basic()

    a = @allocated begin x = 1. end
    println("x = 1. ALLOC ", a)
    a = @allocated begin x = [1.] end
    println("x = [1.] ALLOC ", a)
    a = @allocated begin x = [1.,2] end
    println("x = [1.,2] ALLOC ", a)
    a = @allocated begin x = [1.,2,3] end
    println("x = [1.,2,3] ALLOC ", a)
    a = @allocated begin x = [1.,2,3,4] end
    println("x = [1.,2,3,4] ALLOC ", a)
    a = @allocated begin x = [1.,2,3,4,5] end
    println("x = [1.,2,3,4,5] ALLOC ", a)
    a = @allocated begin x = [1,2,3,4,5] end
    println("x = [1,2,3,4,5] ALLOC ", a)
    a = @allocated begin x = [1.,2.,3.,4.,5.] end
    println("x = [1.,2.,3.,4.,5.] ALLOC ", a)

end


function test_getters(; warntype=false, grid_size=100, discretization=:trapeze, in_place=false)

    # harcdoded arguments
    a = @allocated begin t_1 = 0. end
    b = @allocated begin x_1 = [1.,0.,1.] end
    #b = @allocated begin x_1 = [0.,0.] end
    c = @allocated begin u_1 = [1.] end
    d = @allocated begin v_1 = .1 end
    println("Allocation for hardcoded t,x,u,v: ",a, " ", b, " ", c, " ", d)

    # getters for arguments
    docp, xu = init(in_place=in_place, grid_size=grid_size, discretization=discretization)
    time_grid = CTDirect.get_time_grid(xu, docp)
    a = @allocated begin t_2 = time_grid[1] end
    b = @allocated begin x_2 = CTDirect.get_state_at_time_step(xu, docp, 1) end
    c = @allocated begin u_2 = CTDirect.get_control_at_time_step(xu, docp, 1) end
    d = @allocated begin v_2 = CTDirect.get_optim_variable(xu, docp) end
    e = @allocated begin xx_2 = docp._x(x_2) end
    f = @allocated begin uu_2 = docp._u(u_2) end
    println("Allocation for getters t,x,u,v: ",a, " ", b, " ", c, " ", d, " and vectorized x, u: ", e, " ", f)

    # dynamics
    a = @allocated begin docp.ocp.dynamics(t_2, xx_2, uu_2, v_2) end; println("Allocation for ocp dynamics (getter vectorized args)+vectorization: ",a+e+f)
    a = @allocated begin docp.dynamics_ext(t_2, x_2, u_2, v_2) end; println("Allocation for docp dynamics_ext (includes vectorization): ",a)

end


function test_unit(;test_get=false, test_dyn=false, test_unit_cons=false, test_obj=true, test_cons=true, test_trans=true, test_solve=true, warntype=false, profile=false, grid_size=100, discretization=:trapeze, in_place=false)
    
    # define problem and variables
    prob, docp, xu = init(in_place=in_place, grid_size=grid_size, discretization=discretization)
    c = fill(666.666, docp.dim_NLP_constraints)
    work = similar(xu, docp.dim_NLP_x)

    # getters
    if test_get
        print("t "); @btime $docp.get_final_time($xu)
        print("t bis"); @btime $docp.get_optim_variable($xu)[$docp.ocp.final_time]
        print("v "); @btime $docp.get_optim_variable($xu)
        print("x "); @btime CTDirect.get_state_at_time_step($xu, $docp, $docp.dim_NLP_steps)
        print("u "); @btime CTDirect.get_control_at_time_step($xu, $docp, $docp.dim_NLP_steps) 
        if warntype
            @code_warntype docp.get_final_time(xu)
            @code_warntype docp.get_optim_variable(xu)
            @code_warntype CTDirect.get_state_at_time_step(xu, docp, docp.dim_NLP_steps)
            @code_warntype CTDirect.get_control_at_time_step(xu, docp, docp.dim_NLP_steps)
        end
    end

    t = docp.get_final_time(xu)
    v = docp.get_optim_variable(xu)
    x = CTDirect.get_state_at_time_step(xu, docp, docp.dim_NLP_steps)
    u = CTDirect.get_control_at_time_step(xu, docp, docp.dim_NLP_steps)
    f = similar(xu, docp.dim_NLP_x)

    if test_dyn
        if in_place
            print("dynamics_ext"); @btime $docp.dynamics_ext($f, $t, $x, $u, $v)
            warntype && @code_warntype docp.dynamics_ext(f, t, x, u, v)
            Profile.clear_malloc_data()
            docp.dynamics_ext(f, t, x, u, v)
        else
            print("dynamics_ext"); @btime $docp.dynamics_ext($t, $x, $u, $v)
            warntype && @code_warntype docp.dynamics_ext(t, x, u, v)
            Profile.clear_malloc_data()
            docp.dynamics_ext(t, x, u, v)
        end
    end

    if test_unit_cons
        if in_place
            print("u cons"); @btime $docp.control_constraints[2]($c, $t, $u, $v)
            print("x cons"); @btime $docp.state_constraints[2]($c, $t, $x, $v)
            print("xu cons"); @btime $docp.mixed_constraints[2]($c, $t, $x, $u, $v)
            if warntype
                @code_warntype docp.control_constraints[2](c, t, u, v)
                @code_warntype docp.state_constraints[2](c, t, x, v)
                @code_warntype docp.mixed_constraints[2](c, t, x, u, v)
            end
        else
            print("u cons"); @btime $docp.control_constraints[2]($t, $u, $v)
            print("x cons"); @btime $docp.state_constraints[2]($t, $x, $v)
            print("xu cons"); @btime $docp.mixed_constraints[2]($t, $x, $u, $v)
            if warntype
                @code_warntype docp.control_constraints[2](t, u, v)
                @code_warntype docp.state_constraints[2](t, x, v)
                @code_warntype docp.mixed_constraints[2](t, x, u, v)
            end
        end
    end

    # DOCP_objective
    if test_obj
        print("Objective"); @btime CTDirect.DOCP_objective($xu, $docp)
        warntype && @code_warntype CTDirect.DOCP_objective(xu, docp)
        #Profile.clear_malloc_data()
        if profile
            Profile.Allocs.@profile sample_rate=1.0 CTDirect.DOCP_objective(xu, docp)
            results = Profile.Allocs.fetch()
            PProf.Allocs.pprof()
        end
    end

    # DOCP_constraints
    if test_cons
        print("Constraints"); @btime CTDirect.DOCP_constraints!($c, $xu, $docp)
        if any(c.==666.666)
            error("undefined values in constraints ",c)
        end
        warntype && @code_warntype CTDirect.DOCP_constraints!(c, xu, docp)
        #Profile.clear_malloc_data()
        if profile
            Profile.Allocs.@profile sample_rate=1.0 CTDirect.DOCP_constraints!(c, xu, docp)
            results = Profile.Allocs.fetch()
            PProf.Allocs.pprof()
        end
    end

    # transcription
    if test_trans
        print("Transcription"); @btime direct_transcription($prob.ocp, grid_size=$grid_size)
    end

    # solve
    if test_solve
        sol = direct_solve(prob.ocp, display=false, grid_size=grid_size)
        if !isapprox(sol.objective, prob.obj, rtol=1e-2)
            error("objective mismatch: ", sol.objective, " vs ", prob.obj)
        end
        print("Solve"); @btime direct_solve($prob.ocp, display=false, grid_size=$grid_size)
    end

end

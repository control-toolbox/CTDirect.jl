# tests to check allocations in particular
using CTDirect
using CTBase

using LinearAlgebra
using NLPModelsIpopt
using BenchmarkTools
using Profile
using PProf

include("../test/problems/goddard.jl")

#= local version of dynamics
function F0(x, Cd, beta)
    r, v, m = x
    D = Cd * v^2 * exp(-beta * (r - 1))
    return [v, -D / m - 1 / r^2, 0]
end
function F1(x, Tmax, b)
    r, v, m = x
    return [0, Tmax / m, -b * Tmax]
end
Cd = 310
beta = 500
b = 2
Tmax = 3.5
function local_dynamics(t, x, u, v)
    return F0(x, Cd, beta) + u * F1(x, Tmax, b)
end=#
#= compact function is not better...
function compact_dynamics(t, x, u, vv)
    r, v, m = x
    D = Cd * v^2 * exp(-beta * (r - 1))
    return [v,
            -D / m - 1 / r^2 + u * Tmax / m ,
            - b * u * Tmax]
end=#



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


function test_unit(;test_get=false, test_dyn=false, test_unit_cons=false, test_obj=false, test_cons=true, test_trans=false, test_solve=false, warntype=false, grid_size=100, discretization=:trapeze, in_place=false)
    
    # define problem and variables
    if in_place
        prob = goddard_all_inplace()
    else
        prob = goddard_all()
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
    c = fill(666.666, docp.dim_NLP_constraints)
    work = similar(xu, docp.dim_NLP_x)

    # getters
    if test_get
        print("t "); @btime CTDirect.get_final_time($xu, $docp)
        print("t bis"); @btime CTDirect.get_optim_variable($xu, $docp)[$docp.ocp.final_time]
        print("v "); @btime CTDirect.get_optim_variable($xu, $docp)
        print("x "); @btime CTDirect.get_state_at_time_step($xu, $docp, $docp.dim_NLP_steps)
        print("u "); @btime CTDirect.get_control_at_time_step($xu, $docp, $docp.dim_NLP_steps) 
        if warntype
            @code_warntype CTDirect.get_final_time(xu, docp)
            @code_warntype CTDirect.get_optim_variable(xu, docp)
            @code_warntype CTDirect.get_state_at_time_step(xu, docp, docp.dim_NLP_steps)
            @code_warntype CTDirect.get_control_at_time_step(xu, docp, docp.dim_NLP_steps)
        end
    end

    t = CTDirect.get_final_time(xu, docp)
    v = CTDirect.get_optim_variable(xu, docp)
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
        Profile.Allocs.@profile sample_rate=1.0 CTDirect.DOCP_objective(xu, docp)
        results = Profile.Allocs.fetch()
        PProf.Allocs.pprof()
    end

    # DOCP_constraints
    if test_cons
        print("Constraints"); @btime CTDirect.DOCP_constraints!($c, $xu, $docp)
        if any(c.==666.666)
            error("undefined values in constraints ",c)
        end
        warntype && @code_warntype CTDirect.DOCP_constraints!(c, xu, docp)
        #Profile.clear_malloc_data()
        Profile.Allocs.@profile sample_rate=1.0 CTDirect.DOCP_constraints!(c, xu, docp)
        results = Profile.Allocs.fetch()
        PProf.Allocs.pprof()
    end

    # transcription
    if test_trans
        print("Transcription"); @btime direct_transcription($ocp, grid_size=$grid_size)
    end

    # solve
    if test_solve
        sol = direct_solve(ocp, display=false, grid_size=grid_size)
        if !isapprox(sol.objective, prob.obj, rtol=1e-2)
            error("objective mismatch: ", sol.objective, " vs ", prob.obj)
        end
        print("Solve"); @btime direct_solve($ocp, display=false, grid_size=$grid_size)
    end

end

#= OUTPLACE
julia> test_unit()
dynamics_ext  761.696 ns (14 allocations: 608 bytes)
u cons  1.186 μs (17 allocations: 576 bytes)
x cons  681.263 ns (17 allocations: 512 bytes)
xu cons  693.541 ns (20 allocations: 656 bytes)
Objective  160.520 ns (8 allocations: 368 bytes)
Constraints  381.958 μs (8475 allocations: 320.47 KiB)
Transcription  14.954 ms (165217 allocations: 20.43 MiB)
Solve  158.425 ms (2192013 allocations: 121.40 MiB)
julia> test_unit(in_place=true)
dynamics_ext  220.777 ns (16 allocations: 704 bytes)
u cons  415.065 ns (14 allocations: 448 bytes)
x cons  396.473 ns (14 allocations: 448 bytes)
xu cons  423.005 ns (17 allocations: 592 bytes)
Objective  135.592 ns (7 allocations: 352 bytes)
Constraints  210.285 μs (8222 allocations: 335.66 KiB)
Transcription  14.177 ms (165179 allocations: 20.50 MiB)
Solve  114.000 ms (2147182 allocations: 123.22 MiB)
=#
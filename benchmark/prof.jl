# tests to check allocations in particular
using CTDirect
using CTBase

using LinearAlgebra
using NLPModelsIpopt
using BenchmarkTools
using Profile

include("../test/problems/goddard.jl")

# local version of dynamics
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
end
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


function test_unit(;test_get=false, test_dyn=true, test_unit_cons=true, test_obj=true, test_cons=true, test_trans=true, test_solve=true, warntype=false, grid_size=100, discretization=:trapeze, in_place=false)
    
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
        print("vx "); @btime CTDirect.vget_state_at_time_step($xu, $docp, $docp.dim_NLP_steps)
        print("vu "); @btime CTDirect.vget_control_at_time_step($xu, $docp, $docp.dim_NLP_steps) 
        if warntype
            @code_warntype CTDirect.get_final_time(xu, docp)
            @code_warntype CTDirect.get_optim_variable(xu, docp)
            @code_warntype CTDirect.vget_state_at_time_step(xu, docp, docp.dim_NLP_steps)
            @code_warntype CTDirect.vget_control_at_time_step(xu, docp, docp.dim_NLP_steps)
        end
    end

    t = CTDirect.get_final_time(xu, docp)
    v = CTDirect.get_optim_variable(xu, docp)
    vx = CTDirect.vget_state_at_time_step(xu, docp, docp.dim_NLP_steps)
    vu = CTDirect.vget_control_at_time_step(xu, docp, docp.dim_NLP_steps)
    f = similar(xu, docp.dim_NLP_x)

    if test_dyn
        if in_place
            print("dynamics_ext"); @btime $docp.dynamics_ext($f, $t, $vx, $vu, $v)
        else
            print("dynamics_ext"); @btime $docp.dynamics_ext($t, $vx, $vu, $v)
        end
    end

    if test_unit_cons
        if in_place
            print("u cons"); @btime $docp.control_constraints[2]($c, $t, $vu, $v)
            print("x cons"); @btime $docp.state_constraints[2]($c, $t, $vx, $v)
            print("xu cons"); @btime $docp.mixed_constraints[2]($c, $t, $vx, $vu, $v)            
        else
            print("u cons"); @btime $docp.control_constraints[2]($t, $vu, $v)
            print("x cons"); @btime $docp.state_constraints[2]($t, $vx, $v)
            print("xu cons"); @btime $docp.mixed_constraints[2]($t, $vx, $vu, $v)
        end
    end

    # DOCP_objective
    if test_obj
        print("Objective"); @btime CTDirect.DOCP_objective($xu, $docp) 
    end

    # DOCP_constraints
    if test_cons
        print("Constraints"); @btime CTDirect.DOCP_constraints!($c, $xu, $docp)
        if any(c.==666.666)
            error("undefined values in constraints ",c)
        end
    end

    # transcription
    if test_trans
        print("Transcription"); @btime direct_transcription($ocp, grid_size=$grid_size)
    end

    # solve
    if test_solve
        print("Solve"); @btime direct_solve($ocp, display=false, grid_size=$grid_size)
    end

end

#= OUTPLACE
dynamics_ext  736.197 ns (14 allocations: 608 bytes)
u cons  1.311 μs (17 allocations: 576 bytes)
x cons  699.278 ns (17 allocations: 512 bytes)
xu cons  717.809 ns (20 allocations: 656 bytes)
Objective  163.728 ns (8 allocations: 368 bytes)
Constraints  376.813 μs (8878 allocations: 337.77 KiB)
Transcription  17.200 ms (186428 allocations: 21.28 MiB)
Solve  168.787 ms (2272314 allocations: 121.02 MiB)
=#

#= INPLACE
dynamics_ext  223.439 ns (16 allocations: 704 bytes)
u cons  411.290 ns (14 allocations: 448 bytes)
x cons  393.302 ns (14 allocations: 448 bytes)
xu cons  416.201 ns (17 allocations: 592 bytes)
Objective  132.614 ns (7 allocations: 352 bytes)
Constraints  216.509 μs (8322 allocations: 341.91 KiB)
Transcription  16.147 ms (185178 allocations: 21.30 MiB)
Solve  127.944 ms (2156581 allocations: 120.32 MiB)
=#
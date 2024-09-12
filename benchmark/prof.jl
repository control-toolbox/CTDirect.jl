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


function test_unit(;test_get=false, test_dyn=true, test_unit_cons=true, test_obj=true, test_cons=true, test_trans=true, test_solve=true, grid_size=100, discretization=:trapeze, in_place=false)
    
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
    #docp,_ = direct_transcription(ocp, grid_size=grid_size, discretization=discretization) # nlp creates more allocs !
    docp = CTDirect.DOCP(ocp, grid_size=grid_size, time_grid=CTDirect.__time_grid(), discretization=disc_method)
    xu = CTDirect.DOCP_initial_guess(docp)
    cons = fill(666.666, docp.dim_NLP_constraints)
    work = similar(xu, docp.dim_NLP_x)

    # getters
    #NB same numbers with @allocated
    #t: 16 with getter vs 0 with get_optim_variable[index]... type problem ? note that times are actually type unstable between fixed times / free times OCPs...
    #x: 80 with scal/vec*, 80 with always vec, 0 with view
    #u:  0 with scal*/vec, 0 with view
    #v:  0 with scal*/vec
    if test_get
        print("t "); @btime CTDirect.get_final_time($xu, $docp)
        print("t bis"); @btime CTDirect.get_optim_variable($xu, $docp)[$docp.ocp.final_time]
        #print("x "); @btime CTDirect.get_state_at_time_step($xu, $docp, $docp.dim_NLP_steps)
        print("vx "); @btime CTDirect.vget_state_at_time_step($xu, $docp, $docp.dim_NLP_steps)
        #print("u "); @btime CTDirect.get_control_at_time_step($xu, $docp, $docp.dim_NLP_steps)
        print("vu "); @btime CTDirect.vget_control_at_time_step($xu, $docp, $docp.dim_NLP_steps) 
        print("v "); @btime CTDirect.get_optim_variable($xu, $docp)
    end

    # dynamics (nb ocp.dynamics idem)
    # 384 with standard scal/vec getters (@allocated)
    # 400 with @btime -_- ffs
    # dynamics x u  649.630 ns (9 allocations: 400 bytes)
    # dynamics raw  659.401 ns (7 allocations: 432 bytes)
    # vdynamics vx vu  636.708 ns (7 allocations: 432 bytes)
    # vdynamics vraw  635.243 ns (7 allocations: 432 bytes)
    t = CTDirect.get_final_time(xu, docp)
    #x = CTDirect.get_state_at_time_step(xu, docp, docp.dim_NLP_steps)
    #u = CTDirect.get_control_at_time_step(xu, docp, docp.dim_NLP_steps)
    v = CTDirect.get_optim_variable(xu, docp)
    vx = CTDirect.vget_state_at_time_step(xu, docp, docp.dim_NLP_steps)
    vu = CTDirect.vget_control_at_time_step(xu, docp, docp.dim_NLP_steps)
    if test_dyn
        #print("dynamics x u"); @btime $docp.dynamics($t, $x, $u, $v)
        #print("dynamics raw"); @btime $docp.dynamics(1., [1.,1.,1.], 1., 1.)
        print("dynamics vx vu"); @btime $docp.dynamics($t, $vx, $vu, $v)
        print("dynamics_ext"); @btime $docp.dynamics_ext($t, $vx, $vu, $v)
    end

    if test_unit_cons
        println(typeof(docp.control_constraints[2](t, vu, v)))
        println(typeof(docp.state_constraints[2](t, vx, v)))
        println(typeof(docp.mixed_constraints[2](t, vx, vu, v)))
        print("u cons"); @btime $docp.control_constraints[2]($t, $vu, $v)
        print("x cons"); @btime $docp.state_constraints[2]($t, $vx, $v)
        print("xu cons"); @btime $docp.mixed_constraints[2]($t, $vx, $vu, $v)
    end

    # DOCP_objective
    if test_obj
        print("Objective"); @btime CTDirect.DOCP_objective($xu, $docp) 
    end

    # DOCP_constraints
    if test_cons
        print("Constraints"); @btime CTDirect.DOCP_constraints!($cons, $xu, $docp)
        if any(cons.==666.666)
            error("undefined values in constraints ",cons)
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

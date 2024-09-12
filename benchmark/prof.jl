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


function test_unit(;test_obj=true, test_cons=true, test_dyn=true, test_get=true, grid_size=10, discretization=:trapeze, in_place=false)
    
    # define problem and variables
    if in_place
        prob = goddard_all_inplace()
    else
        prob = goddard_all()
    end
    ocp = prob[:ocp]
    docp,_ = direct_transcription(ocp, grid_size=grid_size, discretization=discretization) # nlp creates more allocs !
    xu = CTDirect.DOCP_initial_guess(docp)
    cons = fill(666.666, docp.dim_NLP_constraints)
    work = similar(xu, docp.dim_NLP_x)
    vdynamics = CTDirect.vectorize(docp.dynamics, docp.dim_OCP_x, docp.dim_NLP_u)

    # getters
    #NB same numbers with @allocated
    #t: 16 with getter vs 0 with get_optim_variable[index]... type problem ? note that times are actually type unstable between fixed times / free times OCPs...
    #x: 80 with scal/vec*, 80 with always vec, 0 with view
    #u:  0 with scal*/vec, 0 with view
    #v:  0 with scal*/vec
    if test_get
        print("t "); @btime CTDirect.get_final_time($xu, $docp)
        print("t bis"); @btime CTDirect.get_optim_variable($xu, $docp)[$docp.ocp.final_time]
        print("x "); @btime CTDirect.get_state_at_time_step($xu, $docp, $docp.dim_NLP_steps)
        print("vx "); @btime CTDirect.vget_state_at_time_step($xu, $docp, $docp.dim_NLP_steps)
        print("u "); @btime CTDirect.get_control_at_time_step($xu, $docp, $docp.dim_NLP_steps)
        print("vu "); @btime CTDirect.vget_control_at_time_step($xu, $docp, $docp.dim_NLP_steps) 
        print("v "); @btime CTDirect.get_optim_variable($xu, $docp)
    end

    # DOCP_objective
    if test_obj
        CTDirect.DOCP_objective(xu, docp) # compile
        Profile.clear_malloc_data()
        a = @allocated begin
            obj = CTDirect.DOCP_objective(xu, docp)
        end
        println("DOCP_objective ", a)
    end

    # DOCP_constraints
    if test_cons
        CTDirect.DOCP_constraints!(cons, xu, docp) # compile
        if any(cons.==666.666)
            error("undefined values in constraints ",cons)
        end
        #@btime CTDirect.DOCP_constraints!($c, $xu, $docp)
        Profile.clear_malloc_data()
        a = @allocated begin
            CTDirect.DOCP_constraints!(cons, xu, docp)
        end
        println("DOCP_constraints! ", a)
    end

    # dynamics (nb ocp.dynamics idem)
    # 384 with standard scal/vec getters (@allocated)
    # 400 with @btime -_- ffs
    # dynamics x u  649.630 ns (9 allocations: 400 bytes)
    # dynamics raw  659.401 ns (7 allocations: 432 bytes)
    # vdynamics vx vu  636.708 ns (7 allocations: 432 bytes)
    # vdynamics vraw  635.243 ns (7 allocations: 432 bytes)
    t = CTDirect.get_final_time(xu, docp)
    x = CTDirect.get_state_at_time_step(xu, docp, docp.dim_NLP_steps)
    u = CTDirect.get_control_at_time_step(xu, docp, docp.dim_NLP_steps)
    v = CTDirect.get_optim_variable(xu, docp)
    vx = CTDirect.vget_state_at_time_step(xu, docp, docp.dim_NLP_steps)
    vu = CTDirect.vget_control_at_time_step(xu, docp, docp.dim_NLP_steps)
    if test_dyn
        print("dynamics x u"); @btime $docp.dynamics($t, $x, $u, $v)
        print("dynamics raw"); @btime $docp.dynamics(1., [1.,1.,1.], 1., 1.)
        print("vdynamics vx vu"); @btime $vdynamics($t, $vx, $vu, $v)
        print("vdynamics vraw"); @btime $vdynamics($xu[end], (@view $xu[1:3]), (@view $xu[4:4]), $xu[end])
    end
end

# check vectorize too !
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
    cons = fill(-666.666, docp.dim_NLP_constraints)
    work = similar(xu, docp.dim_NLP_x)

    # getters
    a = @allocated begin
        t = CTDirect.get_final_time(xu, docp)
    end
    b = @allocated begin
        x = CTDirect.get_state_at_time_step(xu, docp, docp.dim_NLP_steps)
    end
    c = @allocated begin
        u = CTDirect.get_control_at_time_step(xu, docp, docp.dim_NLP_steps)
    end
    d = @allocated begin
        v = CTDirect.get_optim_variable(xu, docp)
    end
    println("getters t ", a, " x ", b, " u ", c, " v ", d)


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
        if any(cons.==-666.666)
            error("undefined values in constraints ",cons)
        end
        #@btime CTDirect.DOCP_constraints!($c, $xu, $docp)
        Profile.clear_malloc_data()
        a = @allocated begin
            CTDirect.DOCP_constraints!(cons, xu, docp)
        end
        println("DOCP_constraints! ", a)
    end

    # dynamics
    if test_dyn
        docp.ocp.dynamics(t, x, u, v) # compile
        Profile.clear_malloc_data()
        println(typeof(t))
        a = @allocated begin
            docp.dynamics(t, x, u, v) # 384 ie 6 float64 (args, or vect args + ret ?); 368 sans le t ie 16 de moins, encore le 16 !
            #nb idem ocp.dynamics
        end
        #b = @allocated begin
        #    local_dynamics(t, x, u, v) # 368 comme sans le t...
        #end
        println("dynamics ", a)
    end

end

# check vectorize too !
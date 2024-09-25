# tests to check allocations in particular
using CTDirect
using CTBase

using LinearAlgebra
using NLPModelsIpopt
using BenchmarkTools
using Profile
using PProf
using JET

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

function init(;in_place, grid_size, disc_method)
    if in_place
        prob = goddard_all_inplace()
        #prob = double_integrator_a()
    else
        prob = goddard_all()
        #prob = double_integrator_mintf()
    end
    ocp = prob[:ocp]
    docp = CTDirect.DOCP(ocp, grid_size=grid_size, time_grid=CTDirect.__time_grid(), disc_method=string(disc_method))
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


function test_unit(;test_get=true, test_dyn=false, test_unit_cons=false, test_obj=false, test_block=true, test_cons=true, test_trans=false, test_solve=false, warntype=true, jet=false, profile=false, grid_size=100, disc_method=:trapeze, in_place=false)
    
    # define problem and variables
    prob, docp, xu = init(in_place=in_place, grid_size=grid_size, disc_method=disc_method)
    disc = docp.discretization
    c = fill(666.666, docp.dim_NLP_constraints)
    work = similar(xu, docp.dim_NLP_x)

    # getters
    if test_get
        print("t "); @btime CTDirect.get_final_time($xu, $docp)
        print("v "); @btime CTDirect.get_optim_variable($xu, $docp)
        print("x "); @btime CTDirect.get_state_at_time_step($xu, $docp, $docp.dim_NLP_steps)
        print("u "); @btime CTDirect.get_control_at_time_step($xu, $docp, $docp.dim_NLP_steps) 
        if warntype
            @code_warntype CTDirect.get_final_time(xu, docp)
            @code_warntype CTDirect.get_optim_variable(xu, docp)
            @code_warntype CTDirect.get_state_at_time_step(xu, docp, docp.dim_NLP_steps)
            @code_warntype CTDirect.get_control_at_time_step(xu, docp, docp.dim_NLP_steps)
            @code_warntype CTDirect.get_time_grid!(xu, docp)
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

    # objective
    if test_obj
        print("Objective"); @btime CTDirect.DOCP_objective($xu, $docp)
        warntype && @code_warntype CTDirect.DOCP_objective(xu, docp)
        if profile
            Profile.Allocs.@profile sample_rate=1.0 CTDirect.DOCP_objective(xu, docp)
            results = Profile.Allocs.fetch()
            PProf.Allocs.pprof()
        end
    end

    if test_block
        print("Constraints block")
        CTDirect.get_time_grid!(xu, docp)
        v = CTDirect.get_optim_variable(xu, docp)
        #@code_warntype work = CTDirect.setWorkArray(docp, xu, docp.NLP_time_grid, v)
        i = 1
        #@code_warntype CTDirect.setConstraintBlock!(docp, c, xu, v, docp.NLP_time_grid, i, work)

    end

    # DOCP_constraints
    if test_cons
        print("Constraints"); @btime CTDirect.DOCP_constraints!($c, $xu, $docp)
        if any(c.==666.666)
            error("undefined values in constraints ",c)
        end
        if warntype 
            println("code warntype")
            @code_warntype CTDirect.DOCP_constraints!(c, xu, docp)
            println("code warntype end")
        end
        if jet
            println("JET")
            display(@report_opt CTDirect.DOCP_constraints!(c, xu, docp))
            println("JET end")
        end
        if profile
            println("Profile")
            Profile.Allocs.@profile sample_rate=1.0 CTDirect.DOCP_constraints!(c, xu, docp)
            results = Profile.Allocs.fetch()
            PProf.Allocs.pprof()
            println("Profile end")
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



#= constraints profile
  Total:           0       6998 (flat, cum) 89.08%
     95            .          .               ui = get_control_at_time_step(xu, docp, i) 
     96            .          .            
     97            .          .               #1. state equation 
     98            .          .               if i <= docp.dim_NLP_steps 
     99            .          .                   # more variables 
    100            .        100                   fi = copy(work) # create new copy, not just a reference 
    101            .          .                   tip1 = time_grid[i+1] 
    102            .          .                   xip1 = get_state_at_time_step(xu, docp, i+1) 
    103            .          .                   uip1 = get_control_at_time_step(xu, docp, i+1) 
    104            .          .                   if docp.has_inplace 
    105            .          .                       docp.dynamics_ext(work, tip1, xip1, uip1, v) 
    106            .          .                   else 
    107            .          .                       # copy, do not create a new variable ! 
    108            .       1000                       work[:] = docp.dynamics_ext(tip1, xip1, uip1, v) #+++ jet runtime dispatch here 
    109            .          .                   end 
    110            .          .            
    111            .          .                   # trapeze rule with 'smart' update for dynamics (similar with @.) 
    112            .        800                   c[offset+1:offset+docp.dim_NLP_x] = xip1 - (xi + 0.5 * (tip1 - ti) * (fi + work)) #+++ jet runtime dispatch here even with explicit index ranges 
    113            .          .                   offset += docp.dim_NLP_x 
    114            .          .               end 
    115            .          .            
    116            .          .               # 2. path constraints 
    117            .          .               # Notes on allocations:.= seems similar 
    118            .          .               if docp.dim_u_cons > 0 
    119            .          .                   if docp.has_inplace 
    120            .          .                       docp.control_constraints[2]((@view c[offset+1:offset+docp.dim_u_cons]),ti, docp._u(ui), v) 
    121            .          .                   else 
    122            .       1632                       c[offset+1:offset+docp.dim_u_cons] = docp.control_constraints[2](ti, docp._u(ui), v) 
    123            .          .                   end 
    124            .          .               end 
    125            .          .               if docp.dim_x_cons > 0  
    126            .          .                   if docp.has_inplace 
    127            .          .                       docp.state_constraints[2]((@view c[offset+docp.dim_u_cons+1:offset+docp.dim_u_cons+docp.dim_x_cons]),ti, docp._x(xi), v) 
    128            .          .                   else 
    129            .       1531                       c[offset+docp.dim_u_cons+1:offset+docp.dim_u_cons+docp.dim_x_cons] = docp.state_constraints[2](ti, docp._x(xi), v) 
    130            .          .                   end 
    131            .          .               end 
    132            .          .               if docp.dim_mixed_cons > 0  
    133            .          .                   if docp.has_inplace 
    134            .          .                       docp.mixed_constraints[2]((@view c[offset+docp.dim_u_cons+docp.dim_x_cons+1:offset+docp.dim_u_cons+docp.dim_x_cons+docp.dim_mixed_cons]), ti, docp._x(xi), docp._u(ui), v) 
    135            .          .                   else 
    136            .       1935                       c[offset+docp.dim_u_cons+docp.dim_x_cons+1:offset+docp.dim_u_cons+docp.dim_x_cons+docp.dim_mixed_cons] = docp.mixed_constraints[2](ti, docp._x(xi), docp._u(ui), v) 
    137            .          .                   end 
    138            .          .               end 
    139            .          .            
    140            .          .           end 
    =#
using CTDirect
using CTBase
using Printf
using Statistics
using Plots; pyplot()

test1 = true
test2 = true
test3 = false

#################################################
# goddard max final altitude
println("Test: discrete continuation (goddard)")

Cd = 310
Tmax = 3.5
β = 500
b = 2
function F0(x)
    r, v, m = x
    D = Cd * v^2 * exp(-β*(r - 1))
    return [ v, -D/m - 1/r^2, 0 ]
end
function F1(x)
    r, v, m = x
    return [ 0, Tmax/m, -b*Tmax ]
end

ocp = Model(variable=true)
state!(ocp, 3)
control!(ocp, 1)
variable!(ocp, 1)
time!(ocp, 0, Index(1))
constraint!(ocp, :initial, [1,0,1], :initial_constraint)
constraint!(ocp, :final, Index(3), 0.6, :final_constraint)
constraint!(ocp, :state, 1:2:3, [1,0.6], [1.2,1], :state_box)
constraint!(ocp, :control, Index(1), 0, 1, :control_box)
constraint!(ocp, :variable, Index(1), 0.01, Inf, :variable_box)
constraint!(ocp, :state, Index(2), 0, Inf, :speed_limit)
objective!(ocp, :mayer, (x0, xf, v) -> xf[1], :max)
dynamics!(ocp, (x, u, v) -> F0(x) + u*F1(x) )

# solve unconstrained problem
sol0 = solve(ocp, print_level=0)
@printf("Objective for reference solution %.6f\n", sol0.objective)
sol = sol0

if test1
    # default init
    print("\nSolve for different speed limits, default initial guess\nvmax ")
    iter_list = []
    for vmax=0.15:-0.01:0.05
        print(vmax," ")
        remove_constraint!(ocp, :speed_limit)
        constraint!(ocp, :state, Index(2), 0, vmax, :speed_limit)
        global sol = solve(ocp, print_level=0)   
        #@printf("vmax %.2f objective %.6f iterations %d\n",vmax,sol.objective, sol.iterations)
        push!(iter_list, sol.iterations)
    end
    @printf("\nAverage iterations %d\n", mean(iter_list))

    # warm start
    print("\nDiscrete continuation on speed limit, with warm start\nvmax ")
    vmax_list = []
    obj_list = []
    iter_list = []
    for vmax=0.15:-0.01:0.05
        print(vmax," ")
        remove_constraint!(ocp, :speed_limit)
        constraint!(ocp, :state, Index(2), 0, vmax, :speed_limit)
        global sol = solve(ocp, print_level=0, init=OptimalControlInit(sol))   
        #@printf("vmax %.2f objective %.6f iterations %d\n",vmax,sol.objective, sol.iterations)
        push!(vmax_list, vmax)
        push!(obj_list, sol.objective)
        push!(iter_list, sol.iterations)
    end
    @printf("\nAverage iterations %d\n", mean(iter_list))

    # plot obj(vmax)
    pobj = plot(vmax_list, obj_list, label="r(tf)", xlabel="Speed limit (vmax)", ylabel="Maximal altitude r(tf)",seriestype=:scatter)

    # plot multiple solutions
    plot(sol0)
    p = plot!(sol)

    display(plot(pobj, p, layout=2))
end

if test2
    print("\nDiscrete continuation on maximal thrust\nTmax ")
    # reset speed constraint
    remove_constraint!(ocp, :speed_limit)
    vmax = 0.1
    constraint!(ocp, :state, Index(2), 0, vmax, :speed_limit)
    sol = solve(ocp, print_level=0)
    sol0 = sol

    # continuation on Tmax (using default init is slower)
    Tmax_list = []
    obj_list = []
    for Tmax_local=3.5:-0.5:1
        global Tmax = Tmax_local 
        print(Tmax," ")   
        global sol = solve(ocp, print_level=0, init=OptimalControlInit(sol))
        push!(Tmax_list, Tmax)
        push!(obj_list, sol.objective)
        #print(Tmax, " ", sol.objective, " ", sol.iterations)
    end
    
    # plot obj(vmax)
    pobj = plot(Tmax_list, obj_list, label="r(tf)", xlabel="Maximal thrust (Tmax)", ylabel="Maximal altitude r(tf)",seriestype=:scatter)

    # plot multiple solutions
    plot(sol0)
    p = plot!(sol)
    
    display(plot(pobj, p, layout=2, reuse=false))

end


# currently only one call to time is allowed...
if test3
    time!(ocp, 0, 0.1)
    sol = solve(ocp, print_level=0)
    println(sol.objective)
end
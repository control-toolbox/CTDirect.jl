using CTDirect
using CTBase
using Printf
using Statistics
using Plots; pyplot()

test1 = true
test2 = true
test3 = true
test4 = true

#################################################
# goddard max final altitude
println("Test: discrete continuation")

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
if (test1 || test2)
    sol0 = solve(ocp, print_level=0)
    @printf("Objective for goddard reference solution %.6f",  sol0.objective)
    sol = sol0
end

if test1
    # default init
    print("\nSolve goddard for different speed limits, default initial guess\nvmax ")
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
    sol2 = solve(ocp, print_level=0)
    sol0 = sol2

    # continuation on Tmax (using default init is slower)
    Tmax_list = []
    obj_list = []
    for Tmax_local=3.5:-0.5:1
        global Tmax = Tmax_local 
        print(Tmax," ")   
        global sol2 = solve(ocp, print_level=0, init=OptimalControlInit(sol2))
        push!(Tmax_list, Tmax)
        push!(obj_list, sol2.objective)
        #print(Tmax, " ", sol2.objective, " ", sol2.iterations)
    end
    
    # plot obj(vmax)
    pobj = plot(Tmax_list, obj_list, label="r(tf)", xlabel="Maximal thrust (Tmax)", ylabel="Maximal altitude r(tf)",seriestype=:scatter)

    # plot multiple solutions
    plot(sol0)
    p = plot!(sol2)
    
    display(plot(pobj, p, layout=2, reuse=false))

end


# continuation with parametric definition of the ocp
if test3
    println("\nDiscrete contiuation on parametric ocp")
    # definition of the parametric OCP
    relu(x) = max(0, x)
    μ = 10
    p_relu(x) = log(abs(1 + exp(μ*x)))/μ
    f(x) = 1-x
    m(x) = (p_relu∘f)(x)
    T = 2

    function myocp(ρ)
        @def ocp begin
            τ ∈ R, variable
            s ∈ [ 0, 1 ], time
            x ∈ R², state
            u ∈ R², control
            x₁(0) == 0
            x₂(0) == 1
            x₁(1) == 1
            ẋ(s) == [τ*(u₁(s)+2), (T-τ)*u₂(s)]
            -1 ≤ u₁(s) ≤ 1
            -1 ≤ u₂(s) ≤ 1
            0 ≤ τ ≤ T
            -(x₂(1)-2)^3 - ∫( ρ * ( τ*m(x₁(s))^2 + (T-τ)*m(x₂(s))^2 ) ) → min
        end
        return ocp
    end

    # continuation on rho
    init3 = OptimalControlInit()
    ρs = [0.1, 5, 10, 30, 100]
    for ρ in ρs
        ocp3 = myocp(ρ)
        sol3 = solve(ocp3, print_level=0, init=init3)
        global init3 = OptimalControlInit(sol3)
        @printf("Rho %.2f objective %.6f iterations %d\n", ρ, sol3.objective, sol3.iterations)
    end
end

# continuation on fixed final time
# NB. time! can be called only once, so we redefine the ocp
if (test4)
    println("\nDiscrete continuation on final time (double integrator)")
    # parametric ocp
    function ocp_T(T)
        @def ocp begin
            t ∈ [ 0, T ], time
            x ∈ R², state
            u ∈ R, control
            q = x₁
            v = x₂
            q(0) == 0
            v(0) == 0
            q(T) == 1
            v(T) == 0
            ẋ(t) == [ v(t), u(t) ]
            ∫(u(t)^2) → min
        end
        return ocp
    end

    # contiuation on final time
    init4 = OptimalControlInit()
    for T=1:5
        ocp4 = ocp_T(T) 
        sol4 = solve(ocp4, print_level=0, init=init4)
        global init4 = OptimalControlInit(sol4)
        @printf("T %.2f objective %.6f iterations %d\n", T, sol4.objective, sol4.iterations)
    end

end
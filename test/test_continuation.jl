using CTDirect
using CTBase
using Printf
using Plots

test1 = true
test2 = true
test3 = true

# continuation on fixed final time
# NB. time! can be called only once, so we redefine the ocp
if test1
    println("\nDiscrete continuation on final time")
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

    # continuation on final time
    init1 = OptimalControlInit()
    iter_list = []
    for T=1:5
        local ocp1 = ocp_T(T) 
        local sol1 = solve(ocp1, print_level=0, init=init1)
        global init1 = sol1
        @printf("T %.2f objective %.6f iterations %d\n", T, sol1.objective, sol1.iterations)
        push!(iter_list, sol1.iterations)
    end
    @printf("Average iterations %.2f\n", sum(iter_list)/length(iter_list))

    # recheck solution (T=2) with explicit / non-uniform grid
    ocpT2 = ocp_T(2)
    solT2 = solve(ocpT2, print_level=0)
    solT2_exp = solve(ocpT2, time_grid=LinRange(0,1,101),print_level=0)
    println("T=2 Check explicit grid ", (solT2.objective==solT2_exp.objective) && (solT2.iterations==solT2_exp.iterations))
    solT2_nonunif = solve(ocpT2, time_grid=[0,0.3,1,1.9,2],print_level=0)
    println("T=2 with non-uniform grid ", solT2_nonunif.objective)
    plot(solT2_nonunif, show=true)

end

# continuation with parametric definition of the ocp
if test2
    println("\nDiscrete continuation on parametric ocp")
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
    init2 = OptimalControlInit()
    iter_list = []
    ρs = [0.1, 5, 10, 30, 100]
    for ρ in ρs
        local ocp2 = myocp(ρ)
        local sol2 = solve(ocp2, print_level=0, init=init2)
        global init2 = sol2
        @printf("Rho %.2f objective %.6f iterations %d\n", ρ, sol2.objective, sol2.iterations)
        push!(iter_list, sol2.iterations)
    end
    @printf("Average iterations %.2f\n", sum(iter_list)/length(iter_list))
end


# goddard max final altitude
if (test3)
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
    r0 = 1
    v0 = 0
    m0 = 1
    mf = 0.6
    x0=[r0,v0,m0]
    vmax = 0.1
    state!(ocp, 3)
    control!(ocp, 1)
    variable!(ocp, 1)
    time!(ocp, t0=0, indf=1)
    constraint!(ocp, :initial, lb=x0, ub=x0)
    constraint!(ocp, :final, rg=3, lb=mf, ub=Inf)
    constraint!(ocp, :state, lb=[r0,v0,mf], ub=[r0+0.2,vmax,m0])
    constraint!(ocp, :control, lb=0, ub=1)
    constraint!(ocp, :variable, lb=0.01, ub=Inf)
    objective!(ocp, :mayer, (x0, xf, v) -> xf[1], :max)
    dynamics!(ocp, (x, u, v) -> F0(x) + u*F1(x) )

    # solve unconstrained problem
    sol0 = solve(ocp, print_level=0)
    @printf("\nObjective for goddard reference solution %.6f",  sol0.objective)

    # using a global variable in ocp definition
    print("\nDiscrete continuation on maximal thrust\nTmax ")
    # continuation on Tmax (using default init is slower)
    Tmax_list = []
    obj_list = []
    iter_list = []
    sol3 = sol0
    for Tmax_local=3.5:-0.5:1
        global Tmax = Tmax_local 
        print(Tmax," ")   
        global sol3 = solve(ocp, print_level=0, init=sol3)
        push!(Tmax_list, Tmax)
        push!(obj_list, sol3.objective)
        push!(iter_list, sol3.iterations)
    end
    @printf("\nAverage iterations %.2f\n", sum(iter_list)/length(iter_list))

    # plot obj(vmax)
    pobj = plot(Tmax_list, obj_list, label="r(tf)", xlabel="Maximal thrust (Tmax)", ylabel="Maximal altitude r(tf)",seriestype=:scatter)
    # plot multiple solutions
    plot(sol0)
    p = plot!(sol3)
    display(plot(pobj, p, layout=2, reuse=false))
end

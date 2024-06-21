println("Test: discrete continuation")

# continuation on final time with parametric ocp definition 
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
@testset verbose = true showtiming = true ":parametric_ocp :warm_start" begin
    init = OptimalControlInit()
    obj_list = []
    iter_list = []
    for T=1:5
        ocp = ocp_T(T) 
        sol = solve(ocp, print_level=0, init=init)
        init = sol
        push!(obj_list, sol.objective)
        push!(iter_list, sol.iterations)
    end
    @test last(obj_list) ≈ 0.096038 rtol=1e-2
    @test sum(iter_list)/length(iter_list) ≈ 2.8 rtol=1e-2
end

# recheck solution (T=2) with explicit / non-uniform grid
@testset verbose = true showtiming = true ":explicit_grid :non_uniform_grid" begin
    ocpT2 = ocp_T(2)
    solT2 = solve(ocpT2, print_level=0)
    solT2_exp = solve(ocpT2, time_grid=LinRange(0,1,101),print_level=0)
    @test (solT2.objective==solT2_exp.objective) && (solT2.iterations==solT2_exp.iterations)
    solT2_nonunif = solve(ocpT2, time_grid=[0,0.3,1,1.9,2],print_level=0)
    @test solT2_nonunif.objective ≈ 2.43 rtol=1e-2
end

# continuation with parametric definition of the ocp
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
@testset verbose = true showtiming = true ":parametric_ocp :warm_start" begin
    init = OptimalControlInit()
    obj_list = []
    iter_list = []
    for ρ in [0.1, 5, 10, 30, 100]
        ocp = myocp(ρ)
        sol = solve(ocp, print_level=0, init=init)
        init = sol
        push!(obj_list, sol.objective)
        push!(iter_list, sol.iterations)
    end
    @test last(obj_list) ≈ -148.022367 rtol=1e-2
    @test sum(iter_list)/length(iter_list) ≈ 43.8 rtol=1e-2
end

# goddard
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
sol0 = solve(ocp, print_level=0)

@testset verbose = true showtiming = true ":global_variable :warm_start" begin
    sol = sol0
    obj_list = []
    iter_list = []
    for Tmax_local=3.5:-0.5:1
        global Tmax = Tmax_local 
        sol = solve(ocp, print_level=0, init=sol)
        push!(obj_list, sol.objective)
        push!(iter_list, sol.iterations)
    end
    @test last(obj_list) ≈ 1.00359 rtol=1e-2
    @test sum(iter_list)/length(iter_list) ≈ 16.8 rtol=1e-2
end

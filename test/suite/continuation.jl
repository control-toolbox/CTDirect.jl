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
    for T=1:5
        ocp = ocp_T(T) 
        sol = solve(ocp, print_level=0, init=init)
        init = sol
        push!(obj_list, sol.objective)
    end
    @test obj_list ≈ [12, 1.5, 0.44, 0.19, 0.096] rtol=1e-2
end

#=
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
    for ρ in [0.1, 5, 10, 30, 100]
        ocp = myocp(ρ)
        sol = solve(ocp, print_level=0, init=init)
        init = sol
        push!(obj_list, sol.objective)
    end
    @test obj_list ≈ [-0.034, -1.7, -6.2, -35, -148] rtol=1e-2
end
=#

# goddard
sol0 = solve(goddard, print_level=0)

@testset verbose = true showtiming = true ":global_variable :warm_start" begin
    sol = sol0
    obj_list = []
    for Tmax_local=3.5:-0.5:1
        global Tmax = Tmax_local 
        sol = solve(goddard, print_level=0, init=sol)
        push!(obj_list, sol.objective)
    end
    @test obj_list ≈ [1.0125, 1.0124, 1.0120, 1.0112, 1.0092, 1.0036] rtol=1e-2
end

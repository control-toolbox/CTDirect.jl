using LinearAlgebra

"""
Compute the infty norm : Max_{t_i}(||x(t_i)||_2)

"""
function distance_infty(x::Function, x_sol::Function, T::Vector{Real})
    N = length(T)
    return maximum([norm(x(T[i]) - x_sol(T[i])) for i = 1:N])
end

"""
Compute the L2 norm : Int_{t_0}^{t_f}||u(t)|| dt

"""
function distance_L2(u::Function, u_sol::Function, T::Vector{Real})
    N = length(T)
    dT = T[2:end] - T[1:(end-1)]
    return sum(dT .* [abs(u(T[i]) - u_sol(T[i])) for i ∈ 1:(N-1)])
end

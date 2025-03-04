using CTDirect
using CTBase
using NLPModelsIpopt

T = 50  # Time horizon

# Define the vector field
f(u, v) = [u - u^3 - 10*u*v^2,  -(1 - u^2)*v]
f(x) = f(x...)

eps = 1e-1
asqrt(x) = sqrt(sqrt(x^2 + eps^2))

function lag(x,u)
    fx = f(x)
    unorm2 = u[1]^2 + u[2]^2
    fnorm2 = fx[1]^2 + fx[2]^2
    dotuf = u[1]*fx[1] + u[2]*fx[2]
    return asqrt(unorm2 * fnorm2) - dotuf
end

@def action begin
    t ∈ [0, T], time
    x ∈ R², state
    u ∈ R², control
    x(0) == [-1, 0]    # Starting point (left well)
    x(T) == [1, 0]     # End point (right well)
    ẋ(t) == u(t)       # Path dynamics
    #∫(asqrt(dot(u(t),u(t))*dot(f(x(t)),f(x(t)))) - dot(u(t),f(x(t)))) → min
    ∫(lag(x(t),u(t))) → min
end
# Linear interpolation for x₁
x1(t) = -(1 - t/T) + t/T
# Parabolic guess for x₂
x2(t) = 0.3(-x1(t)^2 + 1)
x(t) = [x1(t), x2(t)]
u(t) = f(x(t))
init = (state=x, control=u)

sol = direct_solve(action; grid_size=1000, disc_method=:gauss_legendre_3, init=init)
plot(sol)
for logeps in 2:5
  global eps = 10.0^-logeps
  global sol = direct_solve(action; grid_size=1000*logeps, disc_method=:gauss_legendre_3, init=sol)
  plot!(sol)
end
current()

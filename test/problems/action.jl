# minimum action problem

function action()

    T = 50

    # Define the vector field
    f(u, v) = [u - u^3 - 10*u*v^2, -(1 - u^2)*v]
    f(x) = f(x...)

    eps = 1e-1
    asqrt(x) = sqrt(sqrt(x^2 + eps^2))

    function lag(x, u)
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
        ∫(lag(x(t), u(t))) → min
    end


    # Init: Linear interpolation for x₁, Parabolic guess for x₂
    x1(t) = -(1 - t/T) + t/T
    x2(t) = 0.3(-x1(t)^2 + 1)
    x(t) = [x1(t), x2(t)]
    u(t) = f(x(t))
    init = (state=x, control=u)

    return ((ocp=action, obj=nothing, name="action", init=init))
end

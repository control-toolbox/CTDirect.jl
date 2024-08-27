# simple intergator

# min enery, dual control
function simple_integrator()
    ocp = Model()
    state!(ocp, 1)
    control!(ocp, 2)
    time!(ocp, t0 = 0, tf = 1)
    constraint!(ocp, :initial, val = -1)
    constraint!(ocp, :final, val = 0)
    constraint!(ocp, :control, lb = [0, 0], ub = [Inf, Inf])
    dynamics!(ocp, (x, u) -> -x - u[1] + u[2])
    objective!(ocp, :lagrange, (x, u) -> (u[1] + u[2])^2)

    return ((ocp = ocp, obj = nothing, name = "simple_integrator", init = nothing))
end

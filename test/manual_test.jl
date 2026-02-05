# manual test, including ie plots
using Plots
Plots.default(show = true)
include("test_common.jl")

include("./problems/beam.jl")
include("./problems/goddard.jl")
include("./problems/double_integrator.jl")

#sol = solve_problem(goddard(); display=true)
#sol1 = solve_problem(goddard2(); display=true, modeler=:exa)

prob = double_integrator_mintf()
v_const = 0.15
t_vec = [0, 0.1, v_const]
x_vec = [[0, 0], [1, 2], [5, -1]]
u_func = t -> (cos(10 * t) + 1) * 0.5
init = @init prob.ocp begin
    x(t_vec) := x_vec
    u(t) := u_func(t)
    tf := v_const
end
sol2 = solve_problem(prob; display=true, max_iter=0, init=init)


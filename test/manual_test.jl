# manual test, including ie plots
using Plots
Plots.default(show = true)
include("test_common.jl")

include("./problems/beam.jl")
include("./problems/goddard.jl")
include("./problems/double_integrator.jl")

prob = double_integrator_mintf()

x_const = [0.5, 0.2]
u_const = 0.5
v_const = 0.15
x_func = t -> [t^2, sqrt(t)]
u_func = t -> (cos(10 * t) + 1) * 0.5
t_vec = [0, 0.1, v_const]
x_vec = [[0, 0], [1, 2], [5, -1]]
x_matrix = [0 0; 1 2; 5 -1]
u_vec = [0, 0.3, 0.1]

init = (state=(t_vec, x_vec), control=(t_vec, u_vec), variable=v_const)

sol = solve_problem(prob; display=true, init=init, max_iter=0)


#plot(sol) even basic plot is broken for all julia versions -_- 

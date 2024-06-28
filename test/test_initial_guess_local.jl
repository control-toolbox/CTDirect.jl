using CTDirect
using CTBase
using NLPModelsIpopt
using HSL
using JLD2
using JSON3
using Printf
using Plots

# double integrator min tf
@def ocp begin
    tf ∈ R, variable
    t ∈ [ 0, tf ], time
    x ∈ R², state
    u ∈ R, control
    -1 ≤ u(t) ≤ 1
    x(0) == [ 0, 0 ]
    x(tf) == [ 1, 0 ]
    0.1 ≤ tf ≤ Inf 
    ẋ(t) == [ x₂(t), u(t) ] 
    tf → min
end

# use cases
t_grid = [0, .1, 1] # Normalised
x_fun(t) = [t, t^2]
x = [ .5, 1 ]
x_matrix = [0 0; 1 2; 5 -1]
x_grid = [[0, 0], [1, 2], [5, -1]]
u_fun(t) = 2t - 1
u = 2
u_grid = [0, 3, .1]
v = 1
N = 400

# UC1: time/matrix for x, function for u, value for v
sol = solve(ocp, init=(time=t_grid, state=x_matrix, control=u_fun, variable=v), print_level=0)
@printf("%-56s %.3f at %d iterations\n", "Init xmatrix, ufunc, vvec", sol.objective, sol.iterations)

# UC2: time/vecvec for x, function for u, value for v
sol = solve(ocp, init=(time=t_grid, state=x_grid, control=u_fun, variable=v), print_level=0)
@printf("%-56s %.3f at %d iterations\n", "Init xfvecvec, ufunc, vvec", sol.objective, sol.iterations)

# ? UC2: function for x, pb with dims u vs N ??
#sol = solve(ocp, init=(grid_size=N, state=x_fun, control=u_grid), print_level=5)

# ? UC3: constant for x, matrix for u (default time grid) SAME PROBLEM
#sol = solve(ocp, init=(state=x, control=u_grid), print_level=5) # Default uniform time grid (N = 100)

# UC4: function for x, constant for u, default for v
sol = solve(ocp, init=(state=x_fun, control=u), print_level=0)
@printf("%-56s %.3f at %d iterations\n", "Init xfunc, uconst, vdefault", sol.objective, sol.iterations)

# UC5: warm start (different grid size from previous solution)
sol = solve(ocp, init=sol, grid_size=N, print_level=0)
@printf("%-56s %.3f at %d iterations\n", "Init warm start", sol.objective, sol.iterations)

#=
# matrix
tf=0.1
steps=11
time = LinRange(0,tf,steps)
x_matrix = stack(x_func.(time),dims=1)
u_matrix = stack(u_func.(time),dims=1)
sol = solve(ocp, print_level=0, init=(time=time, state=x_matrix, control=u_matrix), max_iter=0)
plot(sol, show=true)
sol = solve(ocp, print_level=0, init=(time=time, state=x_matrix, control=u_matrix), max_iter=maxiter)
@printf("%-56s %.3f at %d iterations\n", "Matrix for x, u", sol.objective, sol.iterations)
=#
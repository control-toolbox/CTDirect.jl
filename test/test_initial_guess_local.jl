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
t_grid = [0, .1, 1]
t_matrix = [0 .1 1]
x_fun(t) = [t, t^2]
#x = [ .5, 1 ]
x_matrix = [0 0; 1 2; 5 -1]
x_grid = [[0, 0], [1, 2], [5, -1]]
u_fun(t) = 2t - 1
u = 0.7 # will be projected inside [-1,1] bounds
#u_grid = [0, 3, .1]
v = 1
N = 400

# OK UC1.0: vec for t, vec for x, function for u, value for v
sol = solve(ocp, init=(time=t_grid, state=x_grid, control=u_fun, variable=v), print_level=0)
@printf("%-56s %.3f at %d iterations\n", "Init tvec, xvec, ufunc, vconst", sol.objective, sol.iterations)

# OK UC1.0: mat for t, vec for x, function for u, value for v
sol = solve(ocp, init=(time=t_matrix, state=x_grid, control=u_fun, variable=v), print_level=0)
@printf("%-56s %.3f at %d iterations\n", "Init tmat, xvec, ufunc, vconst", sol.objective, sol.iterations)

# OK UC1.2: vec for t, mat for x, function for u, value for v
sol = solve(ocp, init=(time=t_grid, state=x_matrix, control=u_fun, variable=v), print_level=0)
@printf("%-56s %.3f at %d iterations\n", "Init tvec, xmat, ufunc, vconst", sol.objective, sol.iterations)

# OK UC1.3: matrix for t and x, function for u, value for v
sol = solve(ocp, init=(time=t_matrix, state=x_matrix, control=u_fun, variable=v), print_level=0)
@printf("%-56s %.3f at %d iterations\n", "Init tmat, xmat, ufunc, vconst", sol.objective, sol.iterations)

# ? UC2: function for x, pb with dims u vs N ??
#sol = solve(ocp, init=(grid_size=N, state=x_fun, control=u_grid), print_level=5)

# ? UC3: constant for x, matrix for u (default time grid) SAME PROBLEM
#sol = solve(ocp, init=(state=x, control=u_grid), print_level=5) # Default uniform time grid (N = 100)

# OK UC4: function for x, constant for u, default for v
sol = solve(ocp, init=(state=x_fun, control=u), print_level=0)
@printf("%-56s %.3f at %d iterations\n", "Init xfunc, uconst, vdefault", sol.objective, sol.iterations)

# OK UC5: warm start (different grid size from previous solution)
sol = solve(ocp, init=sol, grid_size=N, print_level=0, max_iter=0)
@printf("%-56s %.3f at %d iterations\n", "Init warm start", sol.objective, sol.iterations)
plot(sol, show=true)
error("stop")

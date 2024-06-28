using CTDirect
using CTBase
using NLPModelsIpopt
using HSL
using JLD2
using JSON3
using Printf
using Plots

# +++ check init with 0 max iter and test value x(tf) u(t)


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

# NB. in ipopt variables will be projected inside bounds before optimization starts
v = 1 #tf
t_grid = [0, .1, v]
t_matrix = [0 .1 v]
x_fun(t) = [t, t^2]
x_const = [ .5, 1 ]
x_matrix = [0 0; 1 2; 5 -1]
x_grid = [[0, 0], [1, 2], [5, -1]]
u_fun(t) = 2t - 1
u_const = 0.7
u_grid = [0, 0.3, .1]

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

# OK UC2: constant for x, vec for u, v const (ie tf=1)
sol = solve(ocp, init=(state=x_const, time=t_grid, control=u_grid, variable=v), print_level=0)
@printf("%-56s %.3f at %d iterations\n", "Init xconst, uvec, vconst", sol.objective, sol.iterations)

# OK UC3: function for x, vec for u, v default (ie tf=0.1 !)
sol = solve(ocp, init=(state=x_fun, time=t_grid, control=u_grid), print_level=0)
@printf("%-56s %.3f at %d iterations\n", "Init xfun, uvec, vdefault", sol.objective, sol.iterations)

# OK UC4: function for x, constant for u
sol = solve(ocp, init=(state=x_fun, control=u_const, variable = v), print_level=0)
@printf("%-56s %.3f at %d iterations\n", "Init xfunc, uconst, vdefault", sol.objective, sol.iterations)

# OK UC5: warm start (different grid size from previous solution)
sol = solve(ocp, init=sol, grid_size=400, print_level=0)
@printf("%-56s %.3f at %d iterations\n", "Init warm start", sol.objective, sol.iterations)


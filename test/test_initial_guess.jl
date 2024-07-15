include("deps.jl")
using Printf
using Plots

println("Test: initial guess options\n")

# use 0 iterations to check initial guess, >0 to check cv
maxiter = 0

# test functions
function check_xf(sol, xf)
    return xf == sol.state(sol.times[end])
end
function check_uf(sol, uf)
    return uf == sol.control(sol.times[end])
end
function check_v(sol, v)
    return v == sol.variable
end

# reference solution
@def ocp begin
    tf ∈ R, variable
    t ∈ [ 0, tf ], time
    x ∈ R², state
    u ∈ R, control
    -1 ≤ u(t) ≤ 1
    x(0) == [ 0, 0 ]
    x(tf) == [ 1, 0 ]
    0.05 ≤ tf ≤ Inf 
    ẋ(t) == [ x₂(t), u(t) ] 
    tf → min
end
sol0 = solve(ocp, print_level=0)

# constant initial guess
x_const = [0.5, 0.2]
u_const = 0.5
v_const = 0.15

# functional initial guess
x_func = t->[t^2, sqrt(t)]
u_func = t->(cos(10*t)+1)*0.5

# interpolated initial gues
x_vec = [[0, 0], [1, 2], [5, -1]]
x_matrix = [0 0; 1 2; 5 -1]
u_vec = [0, 0.3, .1]


#################################################
# 1 Pass initial guess to all-in-one solve call
println("1. Passing the initial guess at the main solve level")

# 1.a default initial guess
sol = solve(ocp, print_level=0, max_iter=maxiter)
if maxiter > 0
    @printf("%-56s %.3f at %d iterations\n", "Default initial guess", sol.objective, sol.iterations)
else
    println(check_xf(sol, [0.1, 0.1]) && check_uf(sol, 0.1) && check_v(sol, 0.1))
end
sol = solve(ocp, print_level=0, init=(), max_iter=maxiter)
if maxiter > 0
    @printf("%-56s %.3f at %d iterations\n", "Default initial guess", sol.objective, sol.iterations)
else
    println(check_xf(sol, [0.1, 0.1]) && check_uf(sol, 0.1) && check_v(sol, 0.1))
end
sol = solve(ocp, print_level=0, init=nothing, max_iter=maxiter)
if maxiter > 0
    @printf("%-56s %.3f at %d iterations\n", "Default initial guess", sol.objective, sol.iterations)
else
    println(check_xf(sol, [0.1, 0.1]) && check_uf(sol, 0.1) && check_v(sol, 0.1))
end

# 1.b constant initial guess
sol = solve(ocp, print_level=0, init=(state=x_const,), max_iter=maxiter)
if maxiter > 0
    @printf("%-56s %.3f at %d iterations\n", "Constant x; default u,v", sol.objective, sol.iterations)
else
    println(check_xf(sol, x_const))
end

sol = solve(ocp, print_level=0, init=(control=u_const,), max_iter=maxiter)
if maxiter > 0
    @printf("%-56s %.3f at %d iterations\n", "Constant u; default x,v", sol.objective, sol.iterations)
else
    println(check_uf(sol, u_const))
end

sol = solve(ocp, print_level=0, init=(variable=v_const,), max_iter=maxiter)
if maxiter > 0
    @printf("%-56s %.3f at %d iterations\n", "Constant v; default x,u", sol.objective, sol.iterations)
else
    println(check_v(sol, v_const))
end

sol = solve(ocp, print_level=0, init=(state=x_const, control=u_const), max_iter=maxiter)
if maxiter > 0
    @printf("%-56s %.3f at %d iterations\n", "Constant x,u; default v", sol.objective, sol.iterations)
else
    println(check_xf(sol, x_const) && check_uf(sol, u_const))
end

sol = solve(ocp, print_level=0, init=(state=x_const, variable=v_const), max_iter=maxiter)
if maxiter > 0
    @printf("%-56s %.3f at %d iterations\n", "Constant x,v; default u", sol.objective, sol.iterations)
else
    println(check_xf(sol, x_const) && check_v(sol, v_const))
end

sol = solve(ocp, print_level=0, init=(control=u_const, variable=v_const), max_iter=maxiter)
if maxiter > 0
    @printf("%-56s %.3f at %d iterations\n", "Constant u,v; default x", sol.objective, sol.iterations)
else
    println(check_uf(sol, u_const) && check_v(sol, v_const))
end

sol = solve(ocp, print_level=0, init=(state=x_const, control=u_const, variable=v_const), max_iter=maxiter)
if maxiter > 0
    @printf("%-56s %.3f at %d iterations\n", "Constant x,u,v", sol.objective, sol.iterations)
else
    println(check_xf(sol, x_const) && check_uf(sol, u_const) && check_v(sol, v_const))
end    

# 1.c functional initial guess
sol = solve(ocp, print_level=0, init=(state=x_func,), max_iter=maxiter)
if maxiter > 0
    @printf("%-56s %.3f at %d iterations\n", "Functional x; default u,v", sol.objective, sol.iterations)
else
    println(check_xf(sol, x_func(sol.times[end])))
end

sol = solve(ocp, print_level=0, init=(control=u_func,), max_iter=maxiter)
if maxiter > 0
    @printf("%-56s %.3f at %d iterations\n", "Functional u; default x,v", sol.objective, sol.iterations)
else
    println(check_uf(sol, u_func(sol.times[end])))
end

sol = solve(ocp, print_level=0, init=(state=x_func, control=u_func), max_iter=maxiter)
if maxiter > 0
    @printf("%-56s %.3f at %d iterations\n", "Functional x,u; default v", sol.objective, sol.iterations)
else
    println(check_xf(sol, x_func(sol.times[end])) && check_uf(sol, u_func(sol.times[end])))
end

# 1.d interpolated initial guess
t_vec = [0, .1, v_const]
sol = solve(ocp, print_level=0, init=(time=t_vec, state=x_vec, control=u_vec, variable=v_const), max_iter=maxiter)
if maxiter > 0
    @printf("%-56s %.3f at %d iterations\n", "Vector t,x,u; constant v", sol.objective, sol.iterations)
else
    println(check_xf(sol, x_vec[end]) && check_uf(sol, u_vec[end]) && check_v(sol, v_const))
end

t_matrix = [0 .1 v_const]
sol = solve(ocp, print_level=0, init=(time=t_matrix, state=x_vec, control=u_vec, variable=v_const), max_iter=maxiter)
if maxiter > 0
    @printf("%-56s %.3f at %d iterations\n", "Matrix t; vector x,u; constant v", sol.objective, sol.iterations)
else
    println(check_xf(sol, x_vec[end]) && check_uf(sol, u_vec[end]) && check_v(sol, v_const))
end

sol = solve(ocp, print_level=0, init=(time=t_vec, state=x_matrix, control=u_vec, variable=v_const), max_iter=maxiter)
if maxiter > 0
    @printf("%-56s %.3f at %d iterations\n", "Matrix x; vector t,u; constant v", sol.objective, sol.iterations)
else
    println(check_xf(sol, x_vec[end]) && check_uf(sol, u_vec[end]) && check_v(sol, v_const))
end

# 1.e mixed initial guess
sol = solve(ocp, print_level=0, init=(time=t_vec, state=x_vec, control=u_func, variable=v_const), max_iter=maxiter)
if maxiter > 0
    @printf("%-56s %.3f at %d iterations\n", "Mixed: vector t,x; functional u; constant v", sol.objective, sol.iterations)
else
    println(check_xf(sol, x_vec[end]) && check_uf(sol, u_func(sol.times[end])) && check_v(sol, v_const))
end

# 1.f warm start
sol = solve(ocp, print_level=0, init=sol0, max_iter=maxiter)
if maxiter > 0
    @printf("%-56s %.3f at %d iterations\n", "Warm start from reference solution", sol.objective, sol.iterations)
else
    println(check_xf(sol, sol.state(sol.times[end])) && check_uf(sol, sol.control(sol.times[end])) && check_v(sol, sol.variable))
end

#################################################
# 2 Setting the initial guess at the DOCP level
println("\n2. Setting the initial guess at the DOCP level")
docp = direct_transcription(ocp)
# mixed init
set_initial_guess(docp, (time=t_vec, state=x_vec, control=u_func, variable=v_const))
dsol = CTDirect.solve_docp(docp, print_level=0, max_iter=maxiter)
sol = build_solution(docp, dsol)
if maxiter > 0
    @printf("%-56s %.3f at %d iterations\n", "Mixed initial guess set in DOCP", sol.objective, sol.iterations)
else
    println(check_xf(sol, x_vec[end]) && check_uf(sol, u_func(sol.times[end])) && check_v(sol, v_const))
end

# warm start
set_initial_guess(docp, sol0)
dsol = CTDirect.solve_docp(docp, print_level=0, max_iter=maxiter)
sol = build_solution(docp, dsol)
if maxiter > 0
    @printf("%-56s %.3f at %d iterations\n", "Warm start set in DOCP", sol.objective, sol.iterations)
else
    println(check_xf(sol, sol.state(sol.times[end])) && check_uf(sol, sol.control(sol.times[end])) && check_v(sol, sol.variable))
end

#################################################
# 3 Passing the initial guess to solve call
println("\n3. Passing the initial guess to solve call")
set_initial_guess(docp, ()) # reset init in docp
# mixed init
dsol = CTDirect.solve_docp(docp, init=(time=t_vec, state=x_vec, control=u_func, variable=v_const), print_level=0, max_iter=maxiter)
sol = build_solution(docp, dsol)
if maxiter > 0
    @printf("%-56s %.3f at %d iterations\n", "Mixed initial guess passed to solve", sol.objective, sol.iterations)
else
    println(check_xf(sol, x_vec[end]) && check_uf(sol, u_func(sol.times[end])) && check_v(sol, v_const))
end

# warm start
dsol = CTDirect.solve_docp(docp, init=sol0, print_level=0, max_iter=maxiter)
sol = build_solution(docp, dsol)
if maxiter > 0
    @printf("%-56s %.3f at %d iterations\n", "Warm start passed to solve", sol.objective, sol.iterations)
else
    println(check_xf(sol, sol.state(sol.times[end])) && check_uf(sol, sol.control(sol.times[end])) && check_v(sol, sol.variable))
end
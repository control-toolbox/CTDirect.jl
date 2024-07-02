include("common_deps.jl")
using Printf
using Plots

println("Test: initial guess options\n")

# use 0 iterations to check initial guess, >0 to check cv
maxiter = 0

# test functions
function check_xf(sol, xf)
    println(xf)
    println(sol.state(sol.times[end]))
    return xf == sol.state(sol.times[end])
end
function check_uf(sol, uf)
    println(uf)
    println(sol.control(sol.times[end]))
    return uf == sol.control(sol.times[end])
end
function check_v(sol, v)
    println(v)
    println(sol.variable)
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
v = 1 #tf
t_vec = [0, .1, v]
t_matrix = [0 .1 v]
x_vec = [[0, 0], [1, 2], [5, -1]]
x_matrix = [0 0; 1 2; 5 -1]
u_vec = [0, 0.3, .1]


#################################################
# 1 Pass initial guess to all-in-one solve call
println("1. Passing the initial guess at the main solve level")
# default init
sol = solve(ocp, print_level=0, max_iter=maxiter)
if maxiter > 0
    @printf("%-56s %.3f at %d iterations\n", "Default initial guess", sol.objective, sol.iterations)
else
    println(check_xf(sol, [0.1, 0.1]))
    println(check_uf(sol, 0.1))
    println(check_v(sol, 0.1))
end

# Constant initial guess (vector for x; default for u,v)
sol = solve(ocp, print_level=0, init=(state=x_const,), max_iter=maxiter)
if maxiter > 0
    @printf("%-56s %.3f at %d iterations\n", "Constant x; default for u,v)", sol.objective, sol.iterations)
else
    println(check_xf(sol, x_const))
end

# Constant initial guess (vector for u; default for x,v)
sol = solve(ocp, print_level=0, init=(control=u_const,), max_iter=maxiter)
@printf("%-56s %.3f at %d iterations\n", "Constant initial guess (vector for u; default for x,v)", sol.objective, sol.iterations)

# Constant initial guess (vector for v; default for x,u)
sol = solve(ocp, print_level=0, init=(variable=v_const,), max_iter=maxiter)
@printf("%-56s %.3f at %d iterations\n", "Constant initial guess (vector for v; default for x,u)", sol.objective, sol.iterations)

# Constant initial guess (vector for x,u; default for v)
sol = solve(ocp, print_level=0, init=(state=x_const, control=u_const), max_iter=maxiter)
@printf("%-56s %.3f at %d iterations\n", "Constant initial guess (vector for x,u; default for v)", sol.objective, sol.iterations)

# Constant initial guess (vector for x,v; default for u)
sol = solve(ocp, print_level=0, init=(state=x_const, variable=v_const), max_iter=maxiter)
@printf("%-56s %.3f at %d iterations\n", "Constant initial guess (vector for x,v; default for u)", sol.objective, sol.iterations)

# Constant initial guess (vector for u,v; default for x)
sol = solve(ocp, print_level=0, init=(control=u_const, variable=v_const), max_iter=maxiter)
@printf("%-56s %.3f at %d iterations\n", "Constant initial guess (vector for u,v; default for x)", sol.objective, sol.iterations)

# Constant initial guess (vector for x,u,v)
sol = solve(ocp, print_level=0, init=(state=x_const, control=u_const, variable=v_const), max_iter=maxiter)
@printf("%-56s %.3f at %d iterations\n", "Constant initial guess x,u,v (compact call)", sol.objective, sol.iterations)


# Functional initial guess for x; default for u,v)
sol = solve(ocp, print_level=0, init=(state=x_func,), max_iter=maxiter)
@printf("%-56s %.3f at %d iterations\n", "Functional initial guess for x; default for u,v", sol.objective, sol.iterations)

# Functional initial guess for u; default for x,v)
sol = solve(ocp, print_level=0, init=(control=u_func,), max_iter=maxiter)
@printf("%-56s %.3f at %d iterations\n", "Functional initial guess for u; default for x,v", sol.objective, sol.iterations)

# Functional initial guess for x,u; default for v)
sol = solve(ocp, print_level=0, init=(state=x_func, control=u_func), max_iter=maxiter)
@printf("%-56s %.3f at %d iterations\n", "Functional x,u; default v (compact call)", sol.objective, sol.iterations)

# Functional for x; constant for u; default for v)
sol = solve(ocp, print_level=0, init=(state=x_func, control=u_const), max_iter=maxiter)
@printf("%-56s %.3f at %d iterations\n", "Mixed functional/constant/default (compact call)", sol.objective, sol.iterations)

# warm start
sol = solve(ocp, print_level=0, init=sol0, max_iter=maxiter)
@printf("%-56s %.3f at %d iterations\n", "Warm start from reference solution", sol.objective, sol.iterations)


#################################################
# 2 Setting the initial guess at the DOCP level
println("\n2. Setting the initial guess at the DOCP level")
docp = directTranscription(ocp)
# mixed init
setInitialGuess(docp, (state=x_func, control=u_const))
dsol = solve(docp, print_level=0, max_iter=maxiter)
sol = OCPSolutionFromDOCP(docp, dsol)
@printf("%-56s %.3f at %d iterations\n", "Mixed initial guess set in DOCP (compact call)", sol.objective, sol.iterations)

# warm start
setInitialGuess(docp, sol0)
dsol = solve(docp, print_level=0, max_iter=maxiter)
sol = OCPSolutionFromDOCP(docp, dsol)
@printf("%-56s %.3f at %d iterations\n", "Warm start set in DOCP (compact call)", sol.objective, sol.iterations)

#################################################
# 3 Passing the initial guess to solve call
println("\n3. Passing the initial guess to solve call")
setInitialGuess(docp, ()) # reset init in docp
# mixed init
dsol = solve(docp, init=(state=x_func, control=u_const), print_level=0, max_iter=maxiter)
sol = OCPSolutionFromDOCP(docp, dsol)
@printf("%-56s %.3f at %d iterations\n", "Mixed initial guess passed to solve (compact call)", sol.objective, sol.iterations)

# warm start
dsol = solve(docp, init=sol0, print_level=0, max_iter=maxiter)
sol = OCPSolutionFromDOCP(docp, dsol)
@printf("%-56s %.3f at %d iterations\n", "Warm start passed to solve (compact call)", sol.objective, sol.iterations)
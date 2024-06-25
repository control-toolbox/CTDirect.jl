using CTDirect
using CTBase
using Printf

#################################################
# goddard max final altitude (all constraint types formulation)
println("Test goddard: initial guess options\n")
include("suite/goddard.jl")

# reference solution
sol0 = solve(ocp, print_level=0)

# use 0 iterations to retrieve initial guess as solution
maxiter = 1000

#################################################
# 1 Pass initial guess to all-in-one solve call
println("1. Passing the initial guess at the main solve level")
# default init
sol = solve(ocp, print_level=0, max_iter=maxiter)
@printf("%-56s %.3f at %d iterations\n", "Default initial guess", sol.objective, sol.iterations)

# constant initial guess
x_const = [1.05, 0.2, 0.8]
u_const = 0.5
v_const = 0.15

# Constant initial guess (vector for x; default for u,v)
sol = solve(ocp, print_level=0, init=(state=x_const,), max_iter=maxiter)
@printf("%-56s %.3f at %d iterations\n", "Constant initial guess (vector for x; default for u,v)", sol.objective, sol.iterations)

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

# functional initial guess
x_func = t->[1+t^2, sqrt(t), 1-t]
u_func = t->(cos(t)+1)*0.5

# Functional initial guess for x; default for u,v)
sol = solve(ocp, print_level=0, init=(state=x_func,), max_iter=maxiter)
@printf("%-56s %.3f at %d iterations\n", "Functional initial guess for x; default for u,v", sol.objective, sol.iterations)

# Functional initial guess for u; default for x,v)
sol = solve(ocp, print_level=0, init=(control=u_func,), max_iter=maxiter)
@printf("%-56s %.3f at %d iterations\n", "Functional initial guess for u; default for x,v", sol.objective, sol.iterations)

# Functional initial guess for x,u; default for v)
sol = solve(ocp, print_level=0, init=(state=x_func, control=u_func), max_iter=maxiter)
@printf("%-56s %.3f at %d iterations\n", "Functional x,u; default v (compact call)", sol.objective, sol.iterations)


# Functional initial guess for x; constant for u; default for v)
sol = solve(ocp, print_level=0, init=(state=x_func, control=u_const), max_iter=maxiter)
@printf("%-56s %.3f at %d iterations\n", "Mixed functional/constant/default (compact call)", sol.objective, sol.iterations)

# warm start
sol = solve(ocp, print_level=0, init=sol0, max_iter=maxiter)
@printf("%-56s %.3f at %d iterations\n", "Warm start from reference solution (compact call)", sol.objective, sol.iterations)

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
setInitialGuess(docp, OptimalControlInit()) # reset init in docp
# mixed init
dsol = solve(docp, init=(state=x_func, control=u_const), print_level=0, max_iter=maxiter)
sol = OCPSolutionFromDOCP(docp, dsol)
@printf("%-56s %.3f at %d iterations\n", "Mixed initial guess passed to solve (compact call)", sol.objective, sol.iterations)

# warm start
dsol = solve(docp, init=sol0, print_level=0, max_iter=maxiter)
sol = OCPSolutionFromDOCP(docp, dsol)
@printf("%-56s %.3f at %d iterations\n", "Warm start passed to solve (compact call)", sol.objective, sol.iterations)
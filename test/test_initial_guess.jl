using CTDirect
using Printf

#################################################
# goddard max final altitude (all constraint types formulation)
println("Test goddard (all constraints): initial guess options\n")
ocp = Model(variable=true)
Cd = 310
Tmax = 3.5
β = 500
b = 2
r0 = 1
v0 = 0
vmax = 0.1
m0 = 1
mf = 0.6
x0 = [ r0, v0, m0 ]
state!(ocp, 3)
control!(ocp, 1)
variable!(ocp, 1)
time!(ocp, 0, Index(1))
# use all possible types of constraints
# initial condition
constraint!(ocp, :initial, x0, :initial_constraint)
# final condition
constraint!(ocp, :final, Index(3), mf, :final_constraint)
# state constraint
constraint!(ocp, :state, (x,v)->x[2], -Inf, vmax, :state_con_v_ub)
# control constraint
constraint!(ocp, :control, (u,v)->u, -Inf, 1, :control_con_u_ub)
# mixed constraint
constraint!(ocp, :mixed, (x,u,v)->x[3], mf, Inf, :mixed_con_m_lb)
# variable constraint
constraint!(ocp, :variable, v->v, -Inf, 10, :variable_con_tf_ubx)
# state box
constraint!(ocp, :state, 1:2, [r0,v0], [r0+0.2, Inf], :state_box_rv)
# control box
constraint!(ocp, :control, Index(1), 0, Inf, :control_box_lb)
# variable box
constraint!(ocp, :variable, Index(1), 0.01, Inf, :variable_box_tfmin)
objective!(ocp, :mayer,  (x0, xf, v) -> xf[1], :max)
function F0(x)
    r, v, m = x
    D = Cd * v^2 * exp(-β*(r - 1))
    return [ v, -D/m - 1/r^2, 0 ]
end
function F1(x)
    r, v, m = x
    return [ 0, Tmax/m, -b*Tmax ]
end
dynamics!(ocp, (x, u, v) -> F0(x) + u*F1(x) )

# reference solution
sol0 = solve(ocp, print_level=0)
#= check ok
tf0 = sol0.variable
println("tf ", tf0, " obj ", sol0.objective)
ini0 = OCPInit(sol0)
println("sol v ", sol0.variable, " init v ", ini0.variable_init)
p1 = plot([0:0.01:tf0],t->sol0.state(t)[1])
p1 = plot!([0:0.01:tf0],t->ini0.state_init(t)[1])
p2 = plot([0:0.01:tf0],t->sol0.state(t)[2])
p2 = plot!([0:0.01:tf0],t->ini0.state_init(t)[2])
p3 = plot([0:0.01:tf0],t->sol0.state(t)[3])
p3 = plot!([0:0.01:tf0],t->ini0.state_init(t)[3])
p4 = plot([0:0.01:tf0],t->sol0.control(t))
p4 = plot!([0:0.01:tf0],t->ini0.control_init(t))
display(plot(p1,p2,p3,p4,layout=4))
=#

# use 0 iterations to retrieve initial guess as solution
maxiter = 1000

#################################################
# 1 Pass initial guess to all-in-one solve call
println("Passing the initial guess at the main solve level")
# default init
sol = solve(ocp, print_level=0, max_iter=maxiter)
@printf("%-56s %.3f at %d iterations\n", "Default initial guess (constant 0.1):", sol.objective, sol.iterations)

# constant initial guess
x_const = [1.05, 0.2, 0.8]
u_const = 0.5
v_const = 0.15

# Constant initial guess (vector for x; default for u,v)
sol = solve(ocp, print_level=0, init=OCPInit(state=x_const), max_iter=maxiter)
@printf("%-56s %.3f at %d iterations\n", "Constant initial guess (vector for x; default for u,v):", sol.objective, sol.iterations)

# Constant initial guess (vector for u; default for x,v)
sol = solve(ocp, print_level=0, init=OCPInit(control=u_const), max_iter=maxiter)
@printf("%-56s %.3f at %d iterations\n", "Constant initial guess (vector for u; default for x,v):", sol.objective, sol.iterations)

# Constant initial guess (vector for v; default for x,u)
sol = solve(ocp, print_level=0, init=OCPInit(variable=v_const), max_iter=maxiter)
@printf("%-56s %.3f at %d iterations\n", "Constant initial guess (vector for v; default for x,u):", sol.objective, sol.iterations)

# Constant initial guess (vector for x,u; default for v)
sol = solve(ocp, print_level=0, init=OCPInit(state=x_const, control=u_const), max_iter=maxiter)
@printf("%-56s %.3f at %d iterations\n", "Constant initial guess (vector for x,u; default for v):", sol.objective, sol.iterations)

# Constant initial guess (vector for x,v; default for u)
sol = solve(ocp, print_level=0, init=OCPInit(state=x_const, variable=v_const), max_iter=maxiter)
@printf("%-56s %.3f at %d iterations\n", "Constant initial guess (vector for x,v; default for u):", sol.objective, sol.iterations)

# Constant initial guess (vector for u,v; default for x)
sol = solve(ocp, print_level=0, init=OCPInit(control=u_const, variable=v_const), max_iter=maxiter)
@printf("%-56s %.3f at %d iterations\n", "Constant initial guess (vector for u,v; default for x):", sol.objective, sol.iterations)

# Constant initial guess (vector for x,u,v)
sol = solve(ocp, print_level=0, init=OCPInit(state=x_const, control=u_const, variable=v_const), max_iter=maxiter)
@printf("%-56s %.3f at %d iterations\n", "Constant initial guess (vector for x,u,v):", sol.objective, sol.iterations)

# functional initial guess
x_func = t->[1+t^2, sqrt(t), 1-t]
u_func = t->(cos(t)+1)*0.5

# Functional initial guess for x; default for u,v)
sol = solve(ocp, print_level=0, init=OCPInit(state=x_func), max_iter=maxiter)
@printf("%-56s %.3f at %d iterations\n", "Functional initial guess for x; default for u,v):", sol.objective, sol.iterations)

# Functional initial guess for u; default for x,v)
sol = solve(ocp, print_level=0, init=OCPInit(control=u_func), max_iter=maxiter)
@printf("%-56s %.3f at %d iterations\n", "Functional initial guess for u; default for x,v):", sol.objective, sol.iterations)

# Functional initial guess for x,u; default for v)
sol = solve(ocp, print_level=0, init=OCPInit(state=x_func, control=u_func), max_iter=maxiter)
@printf("%-56s %.3f at %d iterations\n", "Functional initial guess for x,u; default for v):", sol.objective, sol.iterations)

# Functional initial guess for x; constant for u; default for v)
sol = solve(ocp, print_level=0, init=OCPInit(state=x_func, control=u_const), max_iter=maxiter)
@printf("%-56s %.3f at %d iterations\n", "Functional for x; constant for u; default for v):", sol.objective, sol.iterations)

# warm start
sol = solve(ocp, print_level=0, init=init_function_u = OCPInit(sol0), max_iter=maxiter)
@printf("%-56s %.3f at %d iterations\n", "Warm start from reference solution", sol.objective, sol.iterations)

#################################################
# 2 Setting the initial guess at the DOCP level
println("\nSetting the initial guess at the DOCP level")
docp = directTranscription(ocp)
# constant vector init
setDOCPInit(docp, OCPInit(state=x_const, control=u_const, variable=v_const))
sol = solve(docp, print_level=0, max_iter=maxiter)
@printf("%-56s %.3f at %d iterations\n", "Constant initial guess set in DOCP", sol.objective, sol.iterations)
# mixed init
setDOCPInit(docp, OCPInit(state=x_func, control=u_const))
sol = solve(docp, print_level=0, max_iter=maxiter)
@printf("%-56s %.3f at %d iterations\n", "Func/const/default initial guess set in DOCP", sol.objective, sol.iterations)
# warm start
setDOCPInit(docp, OCPInit(sol0))
sol = solve(docp, print_level=0, max_iter=maxiter)
@printf("%-56s %.3f at %d iterations\n", "Warm start from reference solution set in DOCP", sol.objective, sol.iterations)

#################################################
# 3 Passing the initial guess to solve call
println("\nPassing the initial guess to solve call")
setDOCPInit(docp, OCPInit()) # reset init in docp
# constant vector init
sol = solve(docp, init=OCPInit(state=x_const, control=u_const, variable=v_const), print_level=0, max_iter=maxiter)
@printf("%-56s %.3f at %d iterations\n", "constant initial guess passed to solve", sol.objective, sol.iterations)
# mixed init
sol = solve(docp, init=OCPInit(state=x_func, control=u_const), print_level=0, max_iter=maxiter)
@printf("%-56s %.3f at %d iterations\n", "Func/const/default initial guess passed to solve", sol.objective, sol.iterations)
# warm start
sol = solve(docp, init=OCPInit(sol0), print_level=0, max_iter=maxiter)
@printf("%-56s %.3f at %d iterations\n", "Warm start from reference solution passed to solve", sol.objective, sol.iterations)

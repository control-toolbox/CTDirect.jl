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

#################################################
# 1 Pass initial guess to all-in-one solve call
println("Passing the initial guess at the main solve level")
# default init
sol = solveDirect(ocp, print_level=0)
@printf("%-56s %.3f at %d iterations\n", "Default initial guess (constant 0.1):", sol.objective, sol.iterations)

# constant initial guess (vector / function formats)
x_init = [1.05, 0.2, 0.8]
u_init = 0.5
v_init = 0.15

# Constant initial guess (vector for x; default for u,v)
sol = solveDirect(ocp, print_level=0, init=OptimalControlInit(x_init=x_init))
@printf("%-56s %.3f at %d iterations\n", "Constant initial guess (vector for x; default for u,v):", sol.objective, sol.iterations)

# Constant initial guess (vector for u; default for x,v)
sol = solveDirect(ocp, print_level=0, init=OptimalControlInit(u_init=u_init))
@printf("%-56s %.3f at %d iterations\n", "Constant initial guess (vector for u; default for x,v):", sol.objective, sol.iterations)

# Constant initial guess (vector for v; default for x,u)
sol = solveDirect(ocp, print_level=0, init=OptimalControlInit(v_init=v_init))
@printf("%-56s %.3f at %d iterations\n", "Constant initial guess (vector for v; default for x,u):", sol.objective, sol.iterations)

# Constant initial guess (vector for x,u; default for v)
sol = solveDirect(ocp, print_level=0, init=OptimalControlInit(x_init=x_init, u_init=u_init))
@printf("%-56s %.3f at %d iterations\n", "Constant initial guess (vector for x,u; default for v):", sol.objective, sol.iterations)

# Constant initial guess (vector for x,v; default for u)
sol = solveDirect(ocp, print_level=0, init=OptimalControlInit(x_init=x_init, v_init=v_init))
@printf("%-56s %.3f at %d iterations\n", "Constant initial guess (vector for x,v; default for u):", sol.objective, sol.iterations)

# Constant initial guess (vector for u,v; default for x)
sol = solveDirect(ocp, print_level=0, init=OptimalControlInit(u_init=u_init, v_init=v_init))
@printf("%-56s %.3f at %d iterations\n", "Constant initial guess (vector for u,v; default for x):", sol.objective, sol.iterations)

# Constant initial guess (vector for x,u,v)
sol = solveDirect(ocp, print_level=0, init=OptimalControlInit(x_init=x_init, u_init=u_init, v_init=v_init))
@printf("%-56s %.3f at %d iterations\n", "Constant initial guess (vector for x,u,v):", sol.objective, sol.iterations)


#= redo with non constant functions
# Constant initial guess (function for x; vector for u,v)
sol = solveDirect(ocp, print_level=0, init=OptimalControlInit(x_init=t->x_init, u_init=u_init, v_init=v_init))
@printf("%-56s %.3f at %d iterations\n", "Constant initial guess (function for x; vector for u,v):", sol.objective, sol.iterations)

# Constant initial guess (function for u; vector for x,v)
sol = solveDirect(ocp, print_level=0, init=init_function_u = OptimalControlInit(x_init=x_init, u_init=t->u_init, v_init=v_init))
@printf("%-56s %.3f at %d iterations\n", "Constant initial guess (function for u; vector for x,v):", sol.objective, sol.iterations)
=#

# warm start
sol = solveDirect(ocp, print_level=0, init=init_function_u = OptimalControlInit(sol))
@printf("%-56s %.3f at %d iterations\n", "Warm start from previous solution", sol.objective, sol.iterations)

#################################################
# 2 Setting the initial guess at the DOCP level
println("\nSetting the initial guess at the DOCP level")
docp = directTranscription(ocp)
# constant vector init
setDOCPInit(docp, OptimalControlInit(x_init=x_init, u_init=u_init, v_init=v_init))
sol = solveDOCP(docp, print_level=0)
@printf("%-56s %.3f at %d iterations\n", "Constant initial guess set in DOCP", sol.objective, sol.iterations)
# warm start
setDOCPInit(docp, OptimalControlInit(sol))
sol = solveDOCP(docp, print_level=0)
@printf("%-56s %.3f at %d iterations\n", "Warm start set in DOCP", sol.objective, sol.iterations)

#################################################
# 3 Passing the initial guess to solveDOCP call
println("\nPassing the initial guess to solveDOCP call")
setDOCPInit(docp, OptimalControlInit()) # reset init in docp
# constant vector init
sol = solveDOCP(docp, init=OptimalControlInit(x_init=x_init, u_init=u_init, v_init=v_init), print_level=0)
@printf("%-56s %.3f at %d iterations\n", "constant initial guess passed to solveDOCP", sol.objective, sol.iterations)
# warm start
sol = solveDOCP(docp, init=OptimalControlInit(sol), print_level=0)
@printf("%-56s %.3f at %d iterations\n", "Warm start passed to solveDOCP", sol.objective, sol.iterations)

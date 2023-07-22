#using BenchmarkTools
#using Traceur
using Profile
#using PProf
#using ProfileView
#using JET 

using CTDirect
using CTBase

println("Test: profiling")

# define OCP
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
constraint!(ocp, :initial, x0, :initial_constraint)
constraint!(ocp, :final, Index(3), mf, :final_constraint)
constraint!(ocp, :state, (x,v)->x[2], -Inf, vmax, :state_con_v_ub)
constraint!(ocp, :control, (u,v)->u, -Inf, 1, :control_con_u_ub)
constraint!(ocp, :mixed, (x,u,v)->x[3], mf, Inf, :mixed_con_m_lb)
constraint!(ocp, :variable, v->v, -Inf, 10, :variable_con_tf_ubx)
constraint!(ocp, :state, 1:2, [r0,v0], [r0+0.2, Inf], :state_box_rv)
constraint!(ocp, :control, Index(1), 0, Inf, :control_box_lb)
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

# create NLP functions
grid_size = 1000
init = OptimalControlInit()
ctd = CTDirect.CTDirect_data(ocp, grid_size, init)
xu = CTDirect.initial_guess(ctd)
l_var, u_var = CTDirect.variables_bounds(ctd)
lb, ub = CTDirect.constraints_bounds(ctd)

#=
Notes
- les return xu[i] et xu[i:j] generent effectivement des allocations (cf fichiers .mem avec julia --track-allocation=user)
- utiliser la macro @views (@view ne compile pas) devant chaque appel de getter ne change rien
- dans getter, return view(xu,i:i+2) alloue en fait davantage que return xu[i:i+2] (160 vs 64) ! 
- la version inplace des getters semble faire autant d'allocations... wtf
- attention a la tres petite taille des vecteurs, qui provoque des comportements parfois non lineaires (registres etc ?)
=#

println("Constraints")
println("compilation")
@time CTDirect.ipopt_constraint(xu, ctd)
println("basic cpu / allocs")
@time CTDirect.ipopt_constraint(xu, ctd)
println("code warntype")
@code_warntype CTDirect.ipopt_constraint(xu, ctd)
println("trace allocs")
Profile.clear_malloc_data()
CTDirect.ipopt_constraint(xu, ctd);
println("stfu julia")

#=
+++ todo: encapsuler fonctions ocp pour forcer le vectoriel (et le type xu ?)
=#



# Solver
println("First run for compilation")
@time sol = solve(ocp, grid_size=50, print_level=0, tol=1e-12)
println("Second run for benchmark")
@timev sol = solve(ocp, grid_size=50, print_level=5, tol=1e-12)

#= basic benchmark: goddard with 50 steps (second run, 30% compilation time on first one)
NLP stats: 
205var, 154const eq, 154 const ineq 
31570 nnz jac eq/ineq, 21115 nnz hess (non sparse matrices)
954 nnz jac eq, 154 nnz jac ineq, 561 nnz hess (sparse)

base:                               50iter / 89GB / 64s
inplace constraints:                50iter / 77GB / 58s
AD optimized backend*:              50iter / 387MB / 2.5s
(*slightly worse without inplace)
=#

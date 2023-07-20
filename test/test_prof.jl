#using BenchmarkTools
#using Traceur
using Profile
#using PProf
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
# profile objective function
println("Objective")
println("compilation")
@time CTDirect.ipopt_objective(xu, ctd)
println("basic cpu / allocs")
@time CTDirect.ipopt_objective(xu, ctd)
Profile.clear_malloc_data()
CTDirect.ipopt_objective(xu, ctd);
=#
#println("basic cpu / allocs")
#@time CTDirect.ipopt_objective(xu, ctd)
#println("benchmarktools")
#@btime CTDirect.ipopt_objective(xu, ctd)
#println("traceur") not very useful
#@trace CTDirect.ipopt_objective(xu, ctd)

#= Slice version
Objective
compilation
  0.004582 seconds (691 allocations: 41.600 KiB, 99.13% compilation time)
basic cpu / allocs
  0.000019 seconds (10 allocations: 320 bytes)
benchmarktools
  647.589 ns (10 allocations: 320 bytes)
=#
#= Slice version Julia 1.9
Objective
compilation
  0.062225 seconds (30.95 k allocations: 1.701 MiB, 99.80% compilation time)
basic cpu / allocs
  0.000027 seconds (10 allocations: 320 bytes)
benchmarktools
  668.711 ns (10 allocations: 320 bytes)


Notes
- les return xu[i] et xu[i:j] generent effectivement des allocations (cf fichiers .mem avec julia --track-allocation=user)
- utiliser la macro @views (@view ne compile pas) devant chaque appel de getter ne change rien
- dans getter, return view(xu,i:i+2) alloue en fait davantage que return xu[i:i+2] (160 vs 64) ! 
- la version inplace des getters semble faire autant d'allocations... wtf
- attention a la tres petite taille des vecteurs, qui provoque des comportements parfois non lineaires (registres etc ?)
- dans les constraintes, on a une allocation (32) pour le pas de temps h et pour chaque tip1, pourquoi ? A cause du free tf ?
=#

println("Constraints")
println("compilation")
@time CTDirect.ipopt_constraint(xu, ctd)
println("basic cpu / allocs")
@timev CTDirect.ipopt_constraint(xu, ctd)
#println("code warntype")
#@code_warntype CTDirect.ipopt_constraint(xu, ctd)
#println("trace allocs")
#Profile.clear_malloc_data()
#CTDirect.ipopt_constraint(xu, ctd);

# ProfileView. Pas vraiment lisible non plus -_-
#@profview for i in 1:100 
#  CTDirect.ipopt_constraint(xu, ctd)
#end

# PProf. bof, pas lisible
#Profile.clear()
#@profile CTDirect.ipopt_constraint(xu, ctd)
#pprof()

# Jet
#@report_opt CTDirect.ipopt_constraint(xu, ctd)

#println("stfu julia")


#= SANS TYPAGE SUR LES GETTERS
Constraints
compilation
  0.469292 seconds (422.76 k allocations: 23.478 MiB, 98.11% compilation time)
basic cpu / allocs
  0.007670 seconds (104.15 k allocations: 2.689 MiB)
=#

#= AVEC ELTYPE XU SUR LES XI, UI: idem -_-

REGLER EN PRIORITE LES INSTABILITES DE TYPES ET DYNAMIC DISPATCH 

+++ todo: encapsuler fonctions ocp pour forcer le vectoriel (et le type xu ?)

+++ le getter de x et les appels a dynamics refont des allocs (le controle est ok mais il est ici scalaire, c'est la difference ?)
+++ rq: declarer c en champ de ctd pour ne pas reallouer a chaque eval des contraintes ? idem pour xi, xip1, ui, uip1, fi, fip1 ???
PROBLEME DE TYPE ? ou bien retrouver le type utilise par adnlpmodels !
+++ a l'inverse declarer ctd en non mutable, quitte a splitter ??

+++ nb: on voit des allocs dans les fonctions ocp ! (cf test_prof.jl.mem)


=#



# Solver
@time sol = solve(ocp, grid_size=50, print_level=5, tol=1e-12)
@timev sol = solve(ocp, grid_size=50, print_level=5, tol=1e-12)
# solve a 50 steps: 192GB allocs :D
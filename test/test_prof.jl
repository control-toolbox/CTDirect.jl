#using BenchmarkTools
#using Traceur
#using Profile
#using PProf
using JET

using CTDirect
using CTBase

code_warntype = true
jet = true

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

#=
# Solver
println("First run for compilation @time")
@time sol = solve(ocp, grid_size=50, print_level=0, tol=1e-12)
println("Second run for benchmark @timev")
@timev sol = solve(ocp, grid_size=50, print_level=0, tol=1e-12)
#=
  0.835355 seconds (4.25 M allocations: 386.129 MiB, 3.20% gc time)
elapsed time (ns):  835354582
gc time (ns):       26690593
bytes allocated:    404886024
pool allocs:        4228969
non-pool GC allocs: 21080
malloc() calls:     38
free() calls:       38
minor collections:  4
full collections:   0
=#

# Unit tests for objective and constraints
grid_size=50
init=OptimalControlInit()
ctd = CTDirect_data(ocp, grid_size, init)
xu0 = initial_guess(ctd)

# Types
if code_warntype
    println("@code_warntype objective")
    @code_warntype ipopt_objective(xu0, ctd)
    #=
    Locals
    xf::Any
    x0::Any
    v::Any
    obj::Any
    N::Int64
    tf::Any
    t0::Any
    =#

    println("@code_warntype constraints")
    c = zeros(ctd.dim_NLP_constraints)
    @code_warntype ipopt_constraint!(c, xu0, ctd)
    #=
    @_5::Union{Nothing, Tuple{Int64, Int64}}
    index::Any
    x0::Any
    uf::Any
    xf::Any
    li::Any
    xli::Float64
    fi::Any
    ui::Any
    xi::Any
    ti::Any
    v::Any
    h::Any
    N::Int64
    tf::Any
    t0::Any
    i::Int64
    lip1::Any
    xlip1::Float64
    fip1::Any
    uip1::Any
    xip1::Any
    tip1::Any
    =#

    #=
    l_var, u_var = variables_bounds(ctd)
    lb, ub = constraints_bounds(ctd)
    nlp = ADNLPModel!(xu -> ipopt_objective(xu, ctd), 
                    xu0, 
                    l_var, u_var, 
                    (c, xu) -> ipopt_constraint!(c, xu, ctd), 
                    lb, ub, 
                    backend = :optimized)
    =#
end


# JET
if jet
    println("@report_opobjective")
    @report_opt ipopt_objective(xu0, ctd)
    #=
    ═════ 56 possible errors found ═════
    les int semblent ok (int64)
    =#
end
=#

# new formulation
@time docp = DirectTranscription(ocp, grid_size=50)
@time sol = solveDOCP(docp, print_level=0, tol=1e-12)
@timev docp = DirectTranscription(ocp, grid_size=50)
@timev sol = solveDOCP(docp, print_level=0, tol=1e-12)
#=
0.114651 seconds (573.19 k allocations: 44.455 MiB, 13.90% gc time)
elapsed time (ns):  114651167
gc time (ns):       15931072
bytes allocated:    46614240
pool allocs:        571773
non-pool GC allocs: 1382
malloc() calls:     38
free() calls:       0
minor collections:  1
full collections:   0

0.829493 seconds (3.71 M allocations: 342.543 MiB, 4.38% gc time)
elapsed time (ns):  829493252
gc time (ns):       36322533
bytes allocated:    359182552
pool allocs:        3685640
non-pool GC allocs: 19699
free() calls:       38
minor collections:  4
full collections:   0
=#

xu = initial_guess(docp)
println("@code_warntype objective")
@code_warntype ipopt_objective(xu, docp)
#=
Locals
  xf::Any
  x0::Any
  v::Any
  obj::Any
  N::Any
  tf::Any
  t0::Any
=#
println("@report_opt objective")
@report_opt ipopt_objective(xu, docp)
#=
═════ 50 possible errors found ═════
=#
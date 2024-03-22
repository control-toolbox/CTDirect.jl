using CTDirect
using ADNLPModels, NLPModels
#using BenchmarkTools
#using Traceur
#using Profile
#using PProf
using JET

code_warntype = true
jet = false

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

# full solve
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

nlp = getNLP(docp)
x0 = initial_guess(docp) #println(x0 == nlp.meta.x0) true ok
if (code_warntype == true)
  println("@code_warntype ipopt_objective")
  @code_warntype ipopt_objective(x0, docp)
#=
MethodInstance for CTDirect.ipopt_objective(::Vector{Float64}, ::CTDirect.DOCP)
  from ipopt_objective(xu, docp) @ CTDirect ~/CTDirect.jl/src/problem.jl:273
Arguments
  #self#::Core.Const(CTDirect.ipopt_objective)
  xu::Vector{Float64}
  docp::CTDirect.DOCP
Locals
  xf::Any
  x0::Any
  v::Union{Float64, Vector{Float64}}
  obj::Any
  N::Int64
  tf::Any
  t0::Any
Body::Any
1 ─       Core.NewvarNode(:(xf))
│         Core.NewvarNode(:(x0))
│         Core.NewvarNode(:(v))
│         (t0 = CTDirect.get_initial_time(xu, docp))
│         (tf = CTDirect.get_final_time(xu, docp))
│         (N = Base.getproperty(docp, :dim_NLP_steps))
│         (obj = 0)
│   %8  = Base.getproperty(docp, :has_mayer_cost)::Bool
└──       goto #3 if not %8
2 ─       (v = CTDirect.get_variable(xu, docp))
│         (x0 = CTDirect.get_state_at_time_step(xu, docp, 0))
│         (xf = CTDirect.get_state_at_time_step(xu, docp, N))
│   %13 = obj::Core.Const(0)
│   %14 = Base.getproperty(docp, :ocp)::CTBase.OptimalControlModel
│   %15 = Base.getproperty(%14, :mayer)::Union{Nothing, CTBase.Mayer}
│   %16 = x0::Any
│   %17 = xf::Any
│   %18 = (%15)(%16, %17, v)::Real
└──       (obj = %13 + %18)
3 ┄ %20 = Base.getproperty(docp, :has_lagrange_cost)::Bool
└──       goto #5 if not %20
4 ─ %22 = obj::Any
│   %23 = (N + 1)::Int64
│   %24 = Base.getproperty(docp, :dim_NLP_state)::Int64
│   %25 = (%23 * %24)::Int64
│   %26 = Base.getindex(xu, %25)::Float64
└──       (obj = %22 + %26)
5 ┄ %28 = Base.getproperty(docp, :ocp)::CTBase.OptimalControlModel
│   %29 = CTDirect.is_min(%28)::Bool
└──       goto #7 if not %29
6 ─       return obj
7 ─ %32 = -obj::Any
└──       return %32

=#

  println("@code_warntype obj")
  @code_warntype obj(nlp, x0)
#=
MethodInstance for NLPModels.obj(::ADNLPModel{Float64, Vector{Float64}, Vector{Int64}}, ::Vector{Float64})
  from obj(nlp::ADNLPModel, x::AbstractVector) @ ADNLPModels ~/.julia/packages/ADNLPModels/Q4sHr/src/nlp.jl:533
Arguments
  #self#::Core.Const(NLPModels.obj)
  nlp::ADNLPModel{Float64, Vector{Float64}, Vector{Int64}}
  x::Vector{Float64}
Body::Any
1 ─ %1  = NLPModels.length(x)::Int64
│   %2  = Base.getproperty(nlp, :meta)::NLPModelMeta{Float64, Vector{Float64}}
│   %3  = Base.getproperty(%2, :nvar)::Int64
│   %4  = (%1 != %3)::Bool
└──       goto #3 if not %4
2 ─ %6  = Base.getproperty(nlp, :meta)::NLPModelMeta{Float64, Vector{Float64}}
│   %7  = Base.getproperty(%6, :nvar)::Int64
│   %8  = NLPModels.length(x)::Int64
│   %9  = NLPModels.DimensionError("x", %7, %8)::Core.PartialStruct(DimensionError, Any[String, Int64, Int64])
└──       NLPModels.throw(%9)
3 ┄       ADNLPModels.increment!(nlp, :neval_obj)
│   %12 = Base.getproperty(nlp, :f)::Any
│   %13 = (%12)(x)::Any
└──       return %13
=#

  #println("@code_warntype gradient")
  #@code_warntype grad(nlp, x0)
end

if (jet == true)
  println("@report_opt obj")
  @report_opt obj(nlp, x0)
  #═════ 53 possible errors found ═════
end


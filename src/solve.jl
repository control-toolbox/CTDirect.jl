"""
$(TYPEDSIGNATURES)

Solve the optimal control problem

Input : 
ocp : functional description of the optimal control problem (cf. ocp.jl)
N   : number of time steps for the discretization
      Int

Output
sol : solution of the discretized problem
      ...
```@examples
julia> using CTDirect
julia> using CTProblems
julia> ocp =  Problem((:integrator, :dim2, :energy)
julia> solve(ocp)
```
"""

function solve(ocp::OptimalControlModel, 
  description...;
  grid_size::Integer=__grid_size_direct(),
  print_level::Integer=__print_level_ipopt(),
  mu_strategy::String=__mu_strategy_ipopt(),
  display::Bool=__display(),
  init=nothing,
  kwargs...)

  # from OCP to NLP: call ADNLPModel constructor
  ctd = CTDirect_data(ocp, grid_size, init)
  xu0 = initial_guess(ctd)
  l_var, u_var = variables_bounds(ctd)
  lb, ub = constraints_bounds(ctd)
  nlp = ADNLPModel(xu -> ipopt_objective(xu, ctd), xu0, l_var, u_var, xu -> ipopt_constraint(xu, ctd), lb, ub) 

  # solve by IPOPT: more info at 
  # https://github.com/JuliaSmoothOptimizers/NLPModelsIpopt.jl/blob/main/src/NLPModelsIpopt.jl#L119
  # options of ipopt: https://coin-or.github.io/Ipopt/OPTIONS.html
  # callback: https://github.com/jump-dev/Ipopt.jl#solver-specific-callback
  # sb="yes": remove ipopt header
  print_level = display ?  print_level : 0
  ipopt_solution = ipopt(nlp, print_level=print_level, mu_strategy=mu_strategy, sb="yes"; kwargs...)

  # from NLP to OCP: call OptimaControlSolution constructor
  sol = _OptimalControlSolution(ocp, ipopt_solution, ctd)

  return sol

end

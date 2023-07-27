# availble methods by order of preference: from top to bottom
algorithmes = ()
algorithmes = add(algorithmes, (:adnlp, :ipopt))

"""
$(TYPEDSIGNATURES)

Return the list of available methods to solve the optimal control problem.
"""
function Methods()::Tuple{Tuple{Vararg{Symbol}}}
    return algorithmes
end

"""
$(TYPEDSIGNATURES)

Solve the the optimal control problem `ocp` by the method given by the (optional) description.
Return an 
[`OptimalControlSolution`](https://control-toolbox.org/CTBase.jl/stable/api-types.html#CTBase.OptimalControlSolution)
from [`CTBase`](https://github.com/control-toolbox/CTBase.jl) package, that is an approximation of 
the optimal solution if the method has converged correctly.

# The (optional) description

You can pass a partial description.
If you give a partial description, then, if several complete descriptions contains the partial one, 
then, the method with the highest priority is chosen. The higher in the list, the higher is the priority.

Keyword arguments:

- `display`: print or not information during the resolution
- `init`: an initial condition for the solver
- `grid_size`: number of time steps for the discretization
- `print_level`: print level for the `Ipopt` solver
- `mu_strategy`: mu strategy for the `Ipopt` solver

!!! warning

    There is only one available method for the moment: the direct method transforms
    the optimal control problem into a nonlinear programming problem (NLP) solved
    by [`Ipopt`](https://coin-or.github.io/Ipopt/), thanks to the package 
    [`ADNLPModels`](https://github.com/JuliaSmoothOptimizers/ADNLPModels.jl).

!!! tip

    - To see the list of available methods, simply call `Methods()`.
    - You can pass any other option by a pair `keyword=value` according to the chosen method. See for instance, [`Ipopt` options](https://coin-or.github.io/Ipopt/OPTIONS.html).
    - The default values for the keyword arguments are given [here](https://control-toolbox.org/CTDocs.jl/ctbase/stable/api-default.html).

```@examples
julia> solve(ocp)
julia> solve(ocp, :adnlp)
julia> solve(ocp, :adnlp, :ipopt)
julia> solve(ocp, display=false, init=OptimalControlInit(), grid_size=100, print_level=0, mu_strategy="adaptive")
```
"""
function solve(ocp::OptimalControlModel, 
    description...;
    display::Bool=__display(),
    init::OptimalControlInit=OptimalControlInit(),
    grid_size::Integer=__grid_size_direct(),
    print_level::Integer=__print_level_ipopt(),
    mu_strategy::String=__mu_strategy_ipopt(),
    kwargs...)

    # get full description from partial
    # throw error if description is not valid
    # thus, no else is needed below
    method = getFullDescription(description, Methods())

    # Model: from ocp to nlp
    if :adnlp in method
        ctd = CTDirect_data(ocp, grid_size, init)
        xu0 = initial_guess(ctd)
        l_var, u_var = variables_bounds(ctd)
        lb, ub = constraints_bounds(ctd)
        nlp = ADNLPModel!(xu -> ipopt_objective(xu, ctd), 
                          xu0, 
                          l_var, u_var, 
                          (c, xu) -> ipopt_constraint!(c, xu, ctd), 
                          lb, ub, 
                          backend = :optimized)
    end

    # solve
    if :ipopt in method
        # https://github.com/JuliaSmoothOptimizers/NLPModelsIpopt.jl/blob/main/src/NLPModelsIpopt.jl#L119
        # callback: https://github.com/jump-dev/Ipopt.jl#solver-specific-callback
        # sb="yes": remove ipopt header +++ make that default
        print_level = display ?  print_level : 0
        ipopt_solution = ipopt(nlp, print_level=print_level, mu_strategy=mu_strategy, sb="yes"; kwargs...)
    end

    # from NLP to OCP: call OptimaControlSolution constructor
    sol = _OptimalControlSolution(ocp, ipopt_solution, ctd)

return sol

end

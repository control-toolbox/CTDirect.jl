# Build functional OCP solution from discrete DOCP solution

"""
$(TYPEDSIGNATURES)

Check whether a collection `t` is empty or not defined.

# Arguments

- `t`: Any object that may be `nothing` or support `length`.

# Returns

- `::Bool`: `true` if `t` is `nothing` or has length zero, otherwise `false`.

# Example

```julia-repl
julia> is_empty([])
true

julia> is_empty([1, 2, 3])
false

julia> is_empty(nothing)
true
```
"""
is_empty(t) = (isnothing(t) || length(t) == 0)


"""
$(TYPEDSIGNATURES)

Build an OCP functional solution from a DOCP discrete solution given as
a `SolverCore.AbstractExecutionStats` object.

# Arguments

- `docp`: The discretized optimal control problem (`DOCP`).
- `nlp_solution`: A solver execution statistics object.

# Returns

- `solution::CTModels.Solution`: A functional OCP solution containing
  trajectories, multipliers, and solver information.

# Example

```julia-repl
julia> build_OCP_solution(docp, nlp_solution)
CTModels.Solution(...)
```
"""
# +++ todo: replace both parsing functions with series of getter calls
# to unify adnlp / exa cases. the getter takes nlp_solution and a label
function build_OCP_solution(docp::DOCP, nlp_solution::SolverCore.AbstractExecutionStats,
    objective, iterations, constraints_violation, message, status, successful;
    exa_getter = nothing)

    # retrieve time grid +++ unify this ?
    if isnothing(exa_getter)
        T = get_time_grid(nlp_solution.solution, docp)
    else
        T = get_time_grid_exa(nlp_solution, docp, exa_getter)
    end

    # primal variables X, U, v and box multipliers
    # +++ insert code directly here ?
    X, U, v, box_multipliers = parse_DOCP_solution_primal(
        docp, nlp_solution; exa_getter
    )

    # costate and constraints multipliers
    # +++ insert code directly here ?
    P, path_constraints_dual, boundary_constraints_dual = parse_DOCP_solution_dual(
        docp, nlp_solution; exa_getter
    )

    return CTModels.build_solution(
        ocp_model(docp),
        T,
        X,
        U,
        v,
        P;
        objective=objective,
        iterations=iterations,
        constraints_violation=constraints_violation,
        message=message,
        status=status,
        successful=successful,
        path_constraints_dual=path_constraints_dual,
        boundary_constraints_dual=boundary_constraints_dual,
        state_constraints_lb_dual=box_multipliers[1],
        state_constraints_ub_dual=box_multipliers[2],
        control_constraints_lb_dual=box_multipliers[3],
        control_constraints_ub_dual=box_multipliers[4],
        variable_constraints_lb_dual=box_multipliers[5],
        variable_constraints_ub_dual=box_multipliers[6],
    )
end


"""
$(TYPEDSIGNATURES)

Recover OCP state, control, and optimization variables from DOCP primal
variables. Bound multipliers are also parsed if available.

# Arguments

- `docp`: The discretized optimal control problem (`DOCP`).
- `solution`: Array of primal decision variables.
- `multipliers_L`: Lower bound multipliers.
- `multipliers_U`: Upper bound multipliers.
- `nlp_model_backend`: The NLP model backend.
- `nlp_solution`: A solver execution statistics object.

# Returns

- `(X, U, v, box_multipliers)`:  
  - `X`: State trajectory.  
  - `U`: Control trajectory.  
  - `v`: Optimization variables.  
  - `box_multipliers`: Tuple of bound multipliers for states, controls, and variables.

# Example

```julia-repl
julia> X, U, v, box_mults = parse_DOCP_solution_primal(docp, primal;
       multipliers_L=mL, multipliers_U=mU, nlp_model_backend=nlp_model_backend, nlp_solution=nlp_solution)
([...] , [...], [...], (...))
```
"""
function parse_DOCP_solution_primal(
    docp, nlp_solution; exa_getter=nothing
)
    # +++ change layout in CTModels solution builder to remove the transpose here

    # arrays (explicit conversion for GPU case)
    solution = Array(nlp_solution.solution)
    multipliers_L = Array(nlp_solution.multipliers_L)
    multipliers_U = Array(nlp_solution.multipliers_U)

    # state and control variables +++ could be done in getter instead ? check with exa
    N = docp.time.steps
    X = zeros(N + 1, docp.dims.NLP_x)
    U = zeros(N + 1, docp.dims.NLP_u)
    v = zeros(docp.dims.NLP_v)

    # multipliers for box constraints +++ could be done in getter instead ? check with exa
    mult_state_box_lower = zeros(size(X))
    mult_state_box_upper = zeros(size(X))
    mult_control_box_lower = zeros(size(U))
    mult_control_box_upper = zeros(size(U))
    mult_variable_box_lower = zeros(size(v))
    mult_variable_box_upper = zeros(size(v))

    if isnothing(exa_getter)
        # ADNLP 
        getter(nlp_solution; val) = CTDirect.getter(nlp_solution, docp; val)
    else
        # EXA
        getter = exa_getter
    end

    # retrieve state, control and optimization variables
    X[:] = getter(nlp_solution; val=:state)'
    U[:] = getter(nlp_solution; val=:control)'
    v[:] = getter(nlp_solution; val=:variable)

    # lower bounds multiplier
    if !is_empty(multipliers_L)
        mult_state_box_lower[:] = getter(nlp_solution; val=:state_l)'
        mult_control_box_lower[:] = getter(nlp_solution; val=:control_l)'
        mult_variable_box_lower[:] = getter(nlp_solution; val=:variable_l)
    end

    # upper bounds multipliers
    if !is_empty(multipliers_U)
        mult_state_box_upper[:] = getter(nlp_solution; val=:state_u)'
        mult_control_box_upper[:] = getter(nlp_solution; val=:control_u)'
        mult_variable_box_upper[:] = getter(nlp_solution; val=:variable_u)
    end

    box_multipliers = (
        mult_state_box_lower,
        mult_state_box_upper,
        mult_control_box_lower,
        mult_control_box_upper,
        mult_variable_box_lower,
        mult_variable_box_upper,
    )

    return X, U, v, box_multipliers
end

"""
$(TYPEDSIGNATURES)

Recover OCP costates and constraint multipliers from DOCP dual
variables.

# Arguments

- `docp`: The discretized optimal control problem (`DOCP`).
- `multipliers`: Array of dual variables (may be `nothing`).
- `nlp_model_backend`: The NLP model backend (default: `ADNLPBackend()`).
- `nlp_solution`: A solver execution statistics object.

# Returns

- `(P, path_constraints_dual, boundary_constraints_dual)`:  
  - `P`: Costate trajectory.  
  - `path_constraints_dual`: Path constraint multipliers.  
  - `boundary_constraints_dual`: Boundary constraint multipliers.

# Example

```julia-repl
julia> P, path_dual, bound_dual = parse_DOCP_solution_dual(docp, duals; nlp_model_backend=nlp_model_backend, nlp_solution=nlp_solution)
([...] , [...], [...])
```
"""
function parse_DOCP_solution_dual(
    docp, nlp_solution; exa_getter=nothing
)

    # arrays (explicit conversion for GPU case)
    multipliers = Array(nlp_solution.multipliers)
    #isnothing(multipliers) && (multipliers = zeros(docp.dim_NLP_constraints))

    # allocate arrays +++ could be done in getter instead ? check with exa
    N = docp.time.steps
    dpc = docp.dims.path_cons
    dbc = docp.dims.boundary_cons
    P = zeros(N, docp.dims.NLP_x)
    mult_path_constraints = zeros(N + 1, dpc)
    mult_boundary_constraints = zeros(dbc)

    if isnothing(exa_getter)
        #ADNLP
        getter(nlp_solution; val) = CTDirect.getter(nlp_solution, docp; val)
    else
        getter = exa_getter
    end

    # costate
    P[:] = getter(nlp_solution; val=:costate)'

    # path constraints multipliers (+++ not yet implemented for exa)
    isnothing(exa_getter) && (mult_path_constraints = getter(nlp_solution; val=:mult_path_constraints))
    isnothing(exa_getter) && (mult_boundary_constraints = getter(nlp_solution; val=:mult_boundary_constraints))

    return P, mult_path_constraints, mult_boundary_constraints
end

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
function build_OCP_solution(docp::DOCP, nlp_solution::SolverCore.AbstractExecutionStats,
    T, objective, iterations, constraints_violation, message, status, successful;
    exa_getter = nothing)

    # NB. we might get rid of the transpose depending on the CTModels part...

    # arrays (explicit conversion for GPU case)
    solution = Array(nlp_solution.solution)
    multipliers_L = Array(nlp_solution.multipliers_L)
    multipliers_U = Array(nlp_solution.multipliers_U)
    multipliers = Array(nlp_solution.multipliers)

    # primal variables X, U, v and box multipliers
    # state and control variables (allocs seem needed here)
    N = docp.time.steps
    X = zeros(N + 1, docp.dims.NLP_x)
    U = zeros(N + 1, docp.dims.NLP_u)
    v = zeros(docp.dims.NLP_v)

    # multipliers for box constraints
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

    # costate and constraints multipliers
    dpc = docp.dims.path_cons
    dbc = docp.dims.boundary_cons
    P = zeros(N, docp.dims.NLP_x)
    mult_path_constraints = zeros(N + 1, dpc)
    mult_boundary_constraints = zeros(dbc)

    # costate
    P[:] = getter(nlp_solution; val=:costate)'

    # path constraints multipliers (+++ not yet implemented for exa)
    isnothing(exa_getter) && (mult_path_constraints = getter(nlp_solution; val=:mult_path_constraints))
    isnothing(exa_getter) && (mult_boundary_constraints = getter(nlp_solution; val=:mult_boundary_constraints))

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
        path_constraints_dual=mult_path_constraints,
        boundary_constraints_dual=mult_boundary_constraints,
        state_constraints_lb_dual=mult_state_box_lower,
        state_constraints_ub_dual=mult_state_box_upper,
        control_constraints_lb_dual=mult_control_box_lower,
        control_constraints_ub_dual=mult_control_box_upper,
        variable_constraints_lb_dual=mult_variable_box_lower,
        variable_constraints_ub_dual=mult_variable_box_upper,
    )
end


"""
$(TYPEDSIGNATURES)

Retrieve the time grid from the given DOCP solution.

# Arguments

- `nlp_solution`: The DOCP solution.
- `docp`: The DOCP.

# Returns

- `::Vector{Float64}`: The time grid.
"""
function get_time_grid_exa(
    nlp_solution::SolverCore.AbstractExecutionStats, docp::CTDirect.DOCP, exa_getter
)
    grid = zeros(docp.time.steps+1)
    ocp = docp.ocp

    if docp.flags.freet0 || docp.flags.freetf
        v = exa_getter(nlp_solution; val=:variable)
    end

    if docp.flags.freet0
        t0 = CTModels.initial_time(ocp, v)
    else
        t0 = CTModels.initial_time(ocp)
    end
    if docp.flags.freetf
        tf = CTModels.final_time(ocp, v)
    else
        tf = CTModels.final_time(ocp)
    end

    @. grid = t0 + docp.time.normalized_grid * (tf - t0)
    return grid
end
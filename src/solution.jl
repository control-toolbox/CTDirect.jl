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
function build_OCP_solution(docp::DOCP, nlp_solution::SolverCore.AbstractExecutionStats)

    # retrieve NLP model and OCP model
    nlp = nlp_model(docp)
    ocp = ocp_model(docp)

    # retrieve NLP model backend
    nlp_model_backend = docp.nlp_model_backend

    # retrieve data from NLP solver
    objective, iterations, constraints_violation, message, status, successful = SolverInfos(nlp_solution, nlp)

    # arrays (explicit conversion for GPU case)
    solution = Array(nlp_solution.solution)
    multipliers = Array(nlp_solution.multipliers)
    multipliers_L = Array(nlp_solution.multipliers_L)
    multipliers_U = Array(nlp_solution.multipliers_U)

    # time grid
    if nlp_model_backend isa ADNLPBackend
        T = get_time_grid(solution, docp)
    else
        T = get_time_grid_exa(nlp_solution, docp)
    end

    # +++ todo: replace both parsing functions with series of getter calls
    # unify adnlp / exa cases

    # primal variables X, U, v and box multipliers
    X, U, v, box_multipliers = parse_DOCP_solution_primal(
        docp,
        solution;
        multipliers_L=multipliers_L,
        multipliers_U=multipliers_U,
        nlp_model_backend=nlp_model_backend,
        nlp_solution=nlp_solution,
    )

    # costate and constraints multipliers
    P, path_constraints_dual, boundary_constraints_dual = parse_DOCP_solution_dual(
        docp, multipliers; nlp_model_backend=nlp_model_backend, nlp_solution=nlp_solution
    )

    return CTModels.build_solution(
        ocp,
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

Return default convergence information for an NLP solution.

# Returns

- `(objective, iterations, constraints_violation, message, status, successful)`:  
  Default values representing an undefined solver state.

# Example

```julia-repl
julia> SolverInfos()
(0.0, 0, 0.0, "undefined", :undefined, true)
```
"""
function SolverInfos()
    return 0.0, 0, 0.0, "undefined", :undefined, true
end

"""
$(TYPEDSIGNATURES)

Retrieve convergence information from an NLP solution.

# Arguments

- `nlp_solution`: A solver execution statistics object.

# Returns

- `(objective, iterations, constraints_violation, message, status, successful)`:  
  A tuple containing the final objective value, iteration count,
  primal feasibility, solver message, solver status, and success flag.

# Example

```julia-repl
julia> SolverInfos(nlp_solution)
(1.23, 15, 1.0e-6, "Ipopt/generic", :first_order, true)
```
"""
function SolverInfos(nlp_solution::SolverCore.AbstractExecutionStats, ::NLPModels.AbstractNLPModel)
    objective = nlp_solution.objective
    iterations = nlp_solution.iter
    constraints_violation = nlp_solution.primal_feas
    status = nlp_solution.status
    successful = (status == :first_order) || (status == :acceptable)
    return objective, iterations, constraints_violation, "Ipopt/generic", status, successful
end

"""
$(TYPEDSIGNATURES)

Build an OCP functional solution from a DOCP discrete solution, given
explicit primal variables, and optionally dual variables and bound
multipliers.

# Arguments

- `docp`: The discretized optimal control problem (`DOCP`).
- `primal`: Array of primal decision variables.
- `dual`: Array of dual variables (default: `nothing`).
- `multipliers_L`: Lower bound multipliers (default: `nothing`).
- `multipliers_U`: Upper bound multipliers (default: `nothing`).
- `nlp_model_backend`: The NLP model backend (default: `ADNLPBackend()`).
- `nlp_solution`: A solver execution statistics object.

# Returns

- `solution::CTModels.Solution`: A functional OCP solution with
  trajectories, multipliers, and solver information.

# Example

```julia-repl
julia> build_OCP_solution(docp; primal=primal_vars, nlp_solution=nlp_solution)
CTModels.Solution(...)
```
"""
function build_OCP_solution(
    docp;
    primal,
    dual=nothing,
    multipliers_L=nothing,
    multipliers_U=nothing,
    nlp_model_backend=ADNLPBackend(),
    nlp_solution,
)
    ocp = ocp_model(docp)
    solution = primal

    # dummy info
    objective, iterations, constraints_violation, message, status, successful = SolverInfos()

    # recompute objective
    objective = DOCP_objective(solution, docp)

    # time grid
    if nlp_model_backend isa ADNLPBackend
        T = get_time_grid(solution, docp)
    else
        T = get_time_grid_exa(nlp_solution, docp)
    end

    # primal variables X, U, v and box multipliers
    X, U, v, box_multipliers = parse_DOCP_solution_primal(
        docp,
        solution;
        multipliers_L=multipliers_L,
        multipliers_U=multipliers_U,
        nlp_model_backend=nlp_model_backend,
        nlp_solution=nlp_solution,
    )

    # costate and constraints multipliers
    P, path_constraints_dual, boundary_constraints_dual = parse_DOCP_solution_dual(
        docp, dual; nlp_model_backend=nlp_model_backend, nlp_solution=nlp_solution
    )

    return CTModels.build_solution(
        ocp,
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
    docp, solution; multipliers_L, multipliers_U, nlp_model_backend, nlp_solution
)

    # state and control variables
    N = docp.time.steps
    X = zeros(N + 1, docp.dims.OCP_x)
    U = zeros(N + 1, docp.dims.NLP_u)
    v = zeros(docp.dims.NLP_v)

    # multipliers for box constraints
    mult_state_box_lower = zeros(size(X))
    mult_state_box_upper = zeros(size(X))
    mult_control_box_lower = zeros(size(U))
    mult_control_box_upper = zeros(size(U))
    mult_variable_box_lower = zeros(size(v))
    mult_variable_box_upper = zeros(size(v))

    if nlp_model_backend isa ExaBackend # Exa
        getter = docp.exa_getter
        X[:] = getter(nlp_solution; val=:state)' # transpose to match choice below for ADNLP
        U[:] = getter(nlp_solution; val=:control)'
        v[:] = getter(nlp_solution; val=:variable)
        if !is_empty(multipliers_L)
            mult_state_box_lower[:] = getter(nlp_solution; val=:state_l)'
            mult_control_box_lower[:] = getter(nlp_solution; val=:control_l)'
            mult_variable_box_lower[:] = getter(nlp_solution; val=:variable_l)
        end
        if !is_empty(multipliers_U)
            mult_state_box_upper[:] = getter(nlp_solution; val=:state_u)'
            mult_control_box_upper[:] = getter(nlp_solution; val=:control_u)'
            mult_variable_box_upper[:] = getter(nlp_solution; val=:variable_u)
        end

    else # ADNLP

        # replace ipopt 0-length arrays with full 0 arrays
        is_empty(multipliers_L) && (multipliers_L = zeros(docp.dim_NLP_variables))
        is_empty(multipliers_U) && (multipliers_U = zeros(docp.dim_NLP_variables))

        # retrieve optimization variables
        if docp.dims.NLP_v > 0
            v .= get_OCP_variable(solution, docp)
            mult_variable_box_lower .= get_OCP_variable(multipliers_L, docp)
            mult_variable_box_upper .= get_OCP_variable(multipliers_U, docp)
        end

        # state variables and box multipliers
        for i in 1:(N + 1)
            X[i, :] .= get_OCP_state_at_time_step(solution, docp, i)
            mult_state_box_lower[i, :] .= get_OCP_state_at_time_step(multipliers_L, docp, i)
            mult_state_box_upper[i, :] .= get_OCP_state_at_time_step(multipliers_U, docp, i)
        end
        # control variables and box multipliers
        for i in 1:(N + 1)
            U[i, :] .= get_OCP_control_at_time_step(solution, docp, i)
            mult_control_box_lower[i, :] .= get_OCP_control_at_time_step(
                multipliers_L, docp, i
            )
            mult_control_box_upper[i, :] .= get_OCP_control_at_time_step(
                multipliers_U, docp, i
            )
        end
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
function parse_DOCP_solution_dual(docp, multipliers; nlp_model_backend=ADNLPBackend(), nlp_solution)

    # costate
    N = docp.time.steps
    P = zeros(N, docp.dims.NLP_x)

    if nlp_model_backend isa ExaBackend # Exa
        getter = docp.exa_getter
        P[:] = getter(nlp_solution; val=:costate)' # transpose to match choice below for ADNLP
        dpc = docp.dims.path_cons
        dbc = docp.dims.boundary_cons
        mul_path_constraints = zeros(N + 1, dpc) # todo: add getters for path constraints for :exa in CTParser
        mul_boundary_constraints = zeros(dbc) # todo: add getters for boundary constraints for :exa in CTParser

    else # ADNLP
        disc = disc_model(docp)

        # if called with multipliers = nothing, fill with zeros
        isnothing(multipliers) && (multipliers = zeros(docp.dim_NLP_constraints))

        # convert GPU arrays if needed
        multipliers = Array(multipliers)

        # dimensions
        dpc = docp.dims.path_cons
        dbc = docp.dims.boundary_cons

        # constraints multipliers
        mul_path_constraints = zeros(N + 1, dpc)
        mul_boundary_constraints = zeros(dbc)

        # loop over time steps
        i_m = 1
        for i in 1:(N + 1)

            # state equation multiplier for costate
            if i <= N
                P[i, :] = multipliers[i_m:(i_m + docp.dims.NLP_x - 1)]
                # skip state / stage constraints
                i_m += disc._state_stage_eqs_block
            end

            # path constraints and multipliers
            if dpc > 0
                mul_path_constraints[i, :] = multipliers[i_m:(i_m + dpc - 1)]
                i_m += dpc
            end
        end

        # pointwise constraints: boundary then variables
        if dbc > 0
            mul_boundary_constraints[:] = multipliers[i_m:(i_m + dbc - 1)]
            i_m += dbc
        end
    end

    return P, mul_path_constraints, mul_boundary_constraints
end

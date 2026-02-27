
"""
$(TYPEDSIGNATURES)

Compute the objective value of a discretized OCP.

# Arguments

- `xu`: Vector of NLP decision variables.
- `docp::DOCP`: The discretized OCP.

# Returns

- `obj::Float64`: Objective function value.

# Example

```julia-repl
julia> DOCP_objective(xu, docp)
12.34
```
"""
function __objective(xu, docp::DOCP)

    # initialization
    if docp.flags.freet0 || docp.flags.freetf
        time_grid = get_time_grid(xu, docp)
    else
        time_grid = docp.time.fixed_grid
    end
    v = get_OCP_variable(xu, docp)
    ocp = ocp_model(docp)

    # mayer cost
    if docp.flags.mayer
        x0 = get_OCP_state_at_time_step(xu, docp, 1)
        xf = get_OCP_state_at_time_step(xu, docp, docp.time.steps+1)
        obj_mayer = CTModels.mayer(ocp)(x0, xf, v)
    else
        obj_mayer = 0.0
    end

    # lagrange cost
    if docp.flags.lagrange
        obj_lagrange = runningCost(docp, xu, v, time_grid)
    else
        obj_lagrange = 0.0
    end

    # total cost
    obj = obj_mayer + obj_lagrange

    return obj
end

"""
$(TYPEDSIGNATURES)

Compute the nonlinear constraints of a DOCP.

The constraints are modeled as `lb <= C(x) <= ub`.

# Arguments

- `c`: Preallocated constraint vector.
- `xu`: Vector of NLP decision variables.
- `docp::DOCP`: The discretized OCP.

# Returns

- `c`: The filled constraint vector.

# Example

```julia-repl
julia> DOCP_constraints!(zeros(docp.dim_NLP_constraints), xu, docp)
[0.0, 0.1, …]
```
"""
function __constraints!(c, xu, docp::DOCP)

    # initialization
    if docp.flags.freet0 || docp.flags.freetf
        time_grid = get_time_grid(xu, docp)
    else
        time_grid = docp.time.fixed_grid
    end
    v = get_OCP_variable(xu, docp)
    work = setWorkArray(docp, xu, time_grid, v)

    # main loop on time steps
    for i in 1:docp.time.steps
        # state equation (includes stage equation depending on scheme)
        stepStateConstraints!(docp, c, xu, v, time_grid, i, work)
        
        #path constraints
        (docp.dims.path_cons > 0) && stepPathConstraints!(docp, c, xu, v, time_grid, i)
    end
    # path constraints at final time
    (docp.dims.path_cons > 0) && stepPathConstraints!(docp, c, xu, v, time_grid, docp.time.steps+1)

    # boundary constraints
    if docp.dims.boundary_cons > 0
        offset = docp.dim_NLP_constraints - docp.dims.boundary_cons
        ocp = ocp_model(docp)
        x0 = get_OCP_state_at_time_step(xu, docp, 1)
        xf = get_OCP_state_at_time_step(xu, docp, docp.time.steps+1)
        CTModels.boundary_constraints_nl(ocp)[2](
            (@view c[(offset + 1):(offset + docp.dims.boundary_cons)]), x0, xf, v
        )
    end

    # NB. the function *needs* to return c for ADNLPModels
    return c
end


"""
$(TYPEDSIGNATURES)
Set path constraints at given time step
"""
function stepPathConstraints!(docp, c, xu, v, time_grid, i)

    ocp = ocp_model(docp)
    disc = disc_model(docp)

    # skip previous steps
    offset = (i-1)*(disc._state_stage_eqs_block + disc._step_pathcons_block) 
    # skip state equation except at final time
    (i <= docp.time.steps) && (offset += disc._state_stage_eqs_block)
    
    # set constraint
    ti = time_grid[i]
    xi = get_OCP_state_at_time_step(xu, docp, i)
    ui = get_OCP_control_at_time_step(xu, docp, i)
    CTModels.path_constraints_nl(ocp)[2](
    (@view c[(offset + 1):(offset + docp.dims.path_cons)]), ti, xi, ui, v
    )
    return
end


"""
$(TYPEDSIGNATURES)

Build lower and upper bounds vectors for the nonlinear constraints of a DOCP.

# Arguments

- `docp::DOCP`: The discretized OCP.

# Returns

- `(lb, ub)::Tuple{Vector{Float64},Vector{Float64}}`: Lower and upper bounds.

# Example

```julia-repl
julia> constraints_bounds!(docp)
([-1.0, …], [1.0, …])
```
"""
function __constraints_bounds!(docp::DOCP)
    lb = docp.bounds.con_l
    ub = docp.bounds.con_u
    disc = disc_model(docp)
    ocp = ocp_model(docp)

    offset = 0
    for i in 1:(docp.time.steps + 1)
        if i <= docp.time.steps
            # skip (ie leave 0) for state / stage equations 
            offset = offset + disc._state_stage_eqs_block
        end
        # path constraints
        if docp.dims.path_cons > 0
            lb[(offset + 1):(offset + docp.dims.path_cons)] = CTModels.path_constraints_nl(ocp)[1]
            ub[(offset + 1):(offset + docp.dims.path_cons)] = CTModels.path_constraints_nl(ocp)[3]
            offset = offset + docp.dims.path_cons
        end
    end

    # boundary constraints
    if docp.dims.boundary_cons > 0
        lb[(offset + 1):(offset + docp.dims.boundary_cons)] = CTModels.boundary_constraints_nl(ocp)[1]
        ub[(offset + 1):(offset + docp.dims.boundary_cons)] = CTModels.boundary_constraints_nl(ocp)[3]
        offset = offset + docp.dims.boundary_cons
    end

    return lb, ub
end

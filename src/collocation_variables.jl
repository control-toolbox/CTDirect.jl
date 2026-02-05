

"""
$(TYPEDSIGNATURES)

Build an initial guess vector for the discretized OCP.

# Arguments

- `docp::DOCP`: The discretized OCP.
- `init::CTModels.Init`: Initialization settings (default: `CTModels.Init()`).

# Returns

- `NLP_X::Vector{Float64}`: Initial guess vector.

# Example

```julia-repl
julia> DOCP_initial_guess(docp)
[0.1, 0.1, …]
```
"""
function DOCP_initial_guess(docp::DOCP, init::CTModels.OptimalControlInitialGuess=CTModels.initial_guess(docp.ocp))

    # default initialization (internal variables such as lagrange cost, k_i for RK schemes) will keep these default values 
    NLP_X = 0.1 * ones(docp.dim_NLP_variables)

    # set variables if provided (needed first in case of free times !)
    if !isnothing(init.variable)
        set_optim_variable!(NLP_X, init.variable, docp)
    end

    # set state / control variables if provided (final control case handled by setter)
    time_grid = get_time_grid(NLP_X, docp)
    for i in 1:(docp.time.steps + 1)
        ti = time_grid[i]
        set_state_at_time_step!(NLP_X, init.state(ti), docp, i)
        set_control_at_time_step!(NLP_X, init.control(ti), docp, i)
    end

    return NLP_X
end


"""
$(TYPEDSIGNATURES)

Build lower and upper bounds vectors for the variable box constraints of a DOCP.

# Arguments

- `docp::DOCP`: The discretized OCP.

# Returns

- `(var_l, var_u)::Tuple{Vector{Float64},Vector{Float64}}`: Lower and upper bounds for variables.

# Example

```julia-repl
julia> variables_bounds!(docp)
([-Inf, …], [Inf, …])
```
"""
function variables_bounds!(docp::DOCP)
    N = docp.time.steps
    var_l = docp.bounds.var_l
    var_u = docp.bounds.var_u
    ocp = ocp_model(docp)

    # build full ordered sets of bounds
    x_lb, x_ub = build_bounds_block(
        docp.dims.NLP_x,
        CTModels.dim_state_constraints_box(ocp),
        CTModels.state_constraints_box(ocp),
    )
    u_lb, u_ub = build_bounds_block(
        docp.dims.NLP_u,
        CTModels.dim_control_constraints_box(ocp),
        CTModels.control_constraints_box(ocp),
    )

    # set state / control box along time steps
    for i in 1:(N + 1)
        set_state_at_time_step!(var_l, x_lb, docp, i)
        set_state_at_time_step!(var_u, x_ub, docp, i)
        set_control_at_time_step!(var_l, u_lb, docp, i)
        set_control_at_time_step!(var_u, u_ub, docp, i)
    end

    # variable box
    if docp.dims.NLP_v > 0
        v_lb, v_ub = build_bounds_block(
            docp.dims.NLP_v,
            CTModels.dim_variable_constraints_box(ocp),
            CTModels.variable_constraints_box(ocp),
        )
        set_optim_variable!(var_l, v_lb, docp)
        set_optim_variable!(var_u, v_ub, docp)
    end

    return var_l, var_u
end


"""
$(TYPEDSIGNATURES)

Build lower and upper bound vectors for state, control, or optimization variables.

# Arguments

- `dim_var::Int`: Variable dimension.
- `dim_box::Int`: Number of box constraints.
- `box_triplet`: Triplet defining box constraints.

# Returns

- `(x_lb, x_ub)::Tuple{Vector{Float64},Vector{Float64}}`: Lower and upper bounds.

# Example

```julia-repl
julia> build_bounds_block(3, 1, ([0.0], [2], [1.0]))
([-Inf, 0.0, -Inf], [Inf, 1.0, Inf])
```
"""
function build_bounds_block(dim_var, dim_box, box_triplet)
    x_lb = -Inf * ones(dim_var)
    x_ub = Inf * ones(dim_var)
    for j in 1:(dim_box)
        indice = box_triplet[2][j]
        x_lb[indice] = box_triplet[1][j]
        x_ub[indice] = box_triplet[3][j]
    end

    return x_lb, x_ub
end


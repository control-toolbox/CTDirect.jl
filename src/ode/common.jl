#= Common parts for the discretization =#

# Getters

# Generic getter for post optimization parsing
# written for compatibility with examodels getter
function getter(nlp_solution, docp::DOCP; val::Symbol)

    N = docp.time.steps

    # A. dual variables: costate, path / boundary constraints multipliers
    if (val == :costate || val == :mult_path_constraints || val == :mult_boundary_constraints)
        data = nlp_solution.multipliers
        disc = disc_model(docp)
        dpc = docp.dims.path_cons
        dbc = docp.dims.boundary_cons
        P = zeros(docp.dims.NLP_x, N)
        mult_path_constraints = zeros(dpc, N + 1)
        mult_boundary_constraints = zeros(dbc)
        # loop over time grid and get multipliers for state dynamics equation and path constraints
        i_m = 1
        for i in 1:N
            # use state dynamics multipliers for costate
            P[:, i] = data[i_m:(i_m + docp.dims.NLP_x - 1)]
            # skip state / stage constraints 
            i_m += disc._state_stage_eqs_block
            # get path constraints multipliers
            if dpc > 0
                mult_path_constraints[:, i] = data[i_m:(i_m + dpc - 1)]
                i_m += dpc
            end
        end
        # add path constraints at final time
        if dpc > 0
            mult_path_constraints[:, N+1] = data[i_m:(i_m + dpc - 1)]
            i_m += dpc
        end

        # pointwise constraints: boundary then variables
        if dbc > 0
            mult_boundary_constraints[:] = data[i_m:(i_m + dbc - 1)]
            i_m += dbc
        end

        # return required values
        if val == :costate
            return P 
        elseif val == :mult_path_constraints
            return mult_path_constraints
        elseif val == :mult_boundary_constraints
            return mult_boundary_constraints
        else
            error("you should not be here: getter with val ", val)
        end
    end

    
    # B1. primal variables: state, control, optimization variables 
    # B2. or associated box multipliers which use the same layout
    # select data array according to required values
    if occursin("_l",String(val))
        data = nlp_solution.multipliers_L
    elseif occursin("_u",String(val))
        data = nlp_solution.multipliers_U
    else
        data = nlp_solution.solution
    end
    
    # optimization variables
    if val == :variable || val == :variable_l || val == :variable_u
        return get_OCP_variable(data, docp)
    
    # state
    elseif val == :state || val == :state_l || val == :state_u
        V = zeros(docp.dims.NLP_x, N + 1)
        for i in 1:(N + 1)
            V[:, i] .= get_OCP_state_at_time_step(data, docp, i)
        end
        return V
    
    # control
    elseif val == :control || val == :control_l || val == :control_u
        V = zeros(docp.dims.NLP_u, N + 1)
        for i in 1:(N + 1)
            V[:, i] .= get_OCP_control_at_time_step(data, docp, i)
        end
        return V

    else
        error("Unknown val for getter: ", val)
    end
end

"""
$(TYPEDSIGNATURES)

Retrieve optimization variables from the NLP variables.
Convention: stored at the end, hence not dependent on the discretization method
Vector output
"""
function get_OCP_variable(xu, docp::DOCP)
    return @view xu[(docp.dim_NLP_variables - docp.dims.NLP_v + 1):docp.dim_NLP_variables]
end

"""
$(TYPEDSIGNATURES)

Retrieve state variables at given time step from the NLP variables.
Convention: 1 <= i <= dim_NLP_steps+1
Vector output
"""
function get_OCP_state_at_time_step(xu, docp::DOCP, i)
    disc = disc_model(docp)
    offset = (i-1) * disc._step_variables_block
    return @view xu[(offset + 1):(offset + docp.dims.NLP_x)]
end


"""
$(TYPEDSIGNATURES)

Retrieve control variables at given time step from the NLP variables.
Convention: 1 <= i <= dim_NLP_steps(+1), with convention u(tf) = U_N
Vector output
"""
function get_OCP_control_at_time_step(xu, docp::DOCP, i; j=1)
    disc = disc_model(docp)

    # final time case, pick U_N unless U_N+1 is present  
    if !disc._final_control && i == docp.time.steps + 1
        i = docp.time.steps
    end
    
    # skip previous time steps
    offset = (i-1) * disc._step_variables_block + docp.dims.NLP_x
    
    # skip previous controls for this time step
    offset += (j-1) * docp.dims.NLP_u
    
    return @view xu[(offset + 1):(offset + docp.dims.NLP_u)]
end

"""
$(TYPEDSIGNATURES)

Retrieve stage variables at given time step/stage from the NLP variables.
Convention: 1 <= i <= dim_NLP_steps(+1),	1 <= j <= s
Vector output
Note that passing correct indices is up to the caller, no checks are made here.
"""
function get_stagevars_at_time_step(xu, docp::DOCP, i, j)
    disc = disc_model(docp)
    offset = (i-1) * disc._step_variables_block + docp.dims.NLP_x + docp.dims.NLP_u + (j-1)*docp.dims.NLP_x
    return @view xu[(offset + 1):(offset + docp.dims.NLP_x)]
end


# Setters

# NB. full setters for state/control may be useful for initial guess
# but would need to accept both array and functional init data
# while only removing the time steps loop in the calling code...

"""
$(TYPEDSIGNATURES)

Set optimization variables in the NLP variables (for initial guess)
"""
function set_optim_variable!(xu, v_init, docp)
    if !isnothing(v_init)
        xu[(end - docp.dims.NLP_v + 1):end] .= v_init
    end
end


"""
$(TYPEDSIGNATURES)

Set initial guess for state variables at given time step
Convention: 1 <= i <= dim_NLP_steps+1
"""
function set_state_at_time_step!(xu, x_init, docp::DOCP, i)
    if !isnothing(x_init)
        disc = disc_model(docp)
        offset = (i-1) * disc._step_variables_block
        xu[(offset + 1):(offset + docp.dims.NLP_x)] .= x_init
    end
    return
end

"""
$(TYPEDSIGNATURES)

Set initial guess for control variables at given time step
Convention: 1 <= i <= dim_NLP_steps(+1)
"""
function set_control_at_time_step!(xu, u_init, docp::DOCP, i)
    if !isnothing(u_init)
        disc = disc_model(docp)
        if i <= docp.time.steps || (disc._final_control && i <= docp.time.steps + 1)
            offset = (i-1) * disc._step_variables_block + docp.dims.NLP_x
            xu[(offset + 1):(offset + docp.dims.NLP_u)] .= u_init
        end
    end
    return
end

"""
$(TYPEDSIGNATURES)

Set work array for all dynamics evaluations
"""
function setWorkArray(docp::DOCP{<: Scheme}, xu, time_grid, v)
    return nothing
end


"""
$(TYPEDSIGNATURES)

Compute the running cost (must be implemented for each discretization scheme)
"""
function runningCost(docp::DOCP{D}, xu, v, time_grid) where {(D<:Scheme)}
    return integral(docp, xu, v, time_grid, CTModels.lagrange(ocp_model(docp)))
end

"""
$(TYPEDSIGNATURES)

Build sparsity pattern for Jacobian of constraints
(to be implemented for each discretization scheme)
"""
function DOCP_Jacobian_pattern(docp::DOCP{D}) where {(D<:Scheme)}
    error(
        "DOCP_Jacobian_pattern not implemented for discretization ",
        D,
        " Use option solve(...; adnlp_backend=:optimized)",
    )
end

"""
$(TYPEDSIGNATURES)

Build sparsity pattern for Hessian of Lagrangian
(to be implemented for each discretization scheme)
"""
function DOCP_Hessian_pattern(docp::DOCP{D}) where {(D<:Scheme)}
    error(
        "DOCP_Hessian_pattern not implemented for discretization ",
        D,
        " Use option solve(...; adnlp_backend=:optimized)",
    )
end

# utility functions for manual sparsity patterns

"""
$(TYPEDSIGNATURES)

Add block of nonzeros elements to a sparsity pattern 
Format: boolean matrix (M) or index vectors (Is, Js) 
Includes a more compact method for single element case
Option to add the symmetric block also (eg for Hessian)
Note: independent from discretization scheme
"""
function add_nonzero_block!(M, i_start, i_end, j_start, j_end; sym=false)
    M[i_start:i_end, j_start:j_end] .= true
    sym && (M[j_start:j_end, i_start:i_end] .= true)
    return nothing
end
function add_nonzero_block!(M, i, j; sym=false)
    M[i, j] = true
    sym && (M[j, i] = true)
    return nothing
end
function add_nonzero_block!(Is, Js, i_start, i_end, j_start, j_end; sym=false)
    for i in i_start:i_end
        for j in j_start:j_end
            push!(Is, i)
            push!(Js, j)
            sym && push!(Is, j)
            sym && push!(Js, i)
        end
    end
    return nothing
end
function add_nonzero_block!(Is, Js, i, j; sym=false)
    push!(Is, i)
    push!(Js, j)
    sym && push!(Is, j)
    sym && push!(Js, i)
    return nothing
end

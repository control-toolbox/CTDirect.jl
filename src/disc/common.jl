


#= Common parts for the discretization =#



"""
$(TYPEDSIGNATURES)

Retrieve optimization variables from the NLP variables.
Convention: stored at the end, hence not dependent on the discretization method
Vector output
"""
function get_OCP_variable(xu, docp::DOCP)
    return @view xu[(docp.dim_NLP_variables-docp.dims.NLP_v+1):docp.dim_NLP_variables]
end


"""
$(TYPEDSIGNATURES)

Retrieve state variables at given time step from the NLP variables.
Convention: 1 <= i <= dim_NLP_steps+1
Vector output
"""
function get_OCP_state_at_time_step(xu, docp::DOCP, i)
    offset = (i-1) * docp.discretization._step_variables_block
    return @view xu[(offset+1):(offset+docp.dims.OCP_x)]
end
"""
$(TYPEDSIGNATURES)

Retrieve state variable for lagrange cost at given time step from the NLP variables.
Convention: 1 <= i <= dim_NLP_steps+1   (no check for actual lagrange cost presence !)
"""
function get_lagrange_state_at_time_step(xu, docp::DOCP, i)
    offset = (i-1) * docp.discretization._step_variables_block
    return xu[offset+docp.dims.NLP_x]
end

"""
$(TYPEDSIGNATURES)

Retrieve control variables at given time step from the NLP variables.
Convention: 1 <= i <= dim_NLP_steps(+1), with convention u(tf) = U_N
Vector output
"""
function get_OCP_control_at_time_step(xu, docp::DOCP, i)
    # final time case  
    if !docp.discretization._final_control && i == docp.time.steps + 1
        i = docp.time.steps
    end
    offset = (i-1) * docp.discretization._step_variables_block + docp.dims.NLP_x
    return @view xu[(offset+1):(offset+docp.dims.NLP_u)]
end


"""
$(TYPEDSIGNATURES)

Retrieve stage variables at given time step/stage from the NLP variables.
Convention: 1 <= i <= dim_NLP_steps(+1),	1 <= j <= s
Vector output
Note that passing correct indices is up to the caller, no checks are made here.
"""
function get_stagevars_at_time_step(xu, docp::DOCP, i, j)
    offset = (i-1) * docp.discretization._step_variables_block + docp.dims.NLP_x + docp.dims.NLP_u + (j-1)*docp.dims.NLP_x
    return @view xu[(offset+1):(offset+docp.dims.NLP_x)]
end


"""
$(TYPEDSIGNATURES)

Set optimization variables in the NLP variables (for initial guess)
"""
function set_optim_variable!(xu, v_init, docp)
    xu[(end-docp.dims.NLP_v+1):end] .= v_init
end


"""
$(TYPEDSIGNATURES)

Set initial guess for state variables at given time step
Convention: 1 <= i <= dim_NLP_steps+1
"""
function set_state_at_time_step!(xu, x_init, docp::DOCP, i)
    # initialize only actual state variables from OCP (not lagrange state)
    if !isnothing(x_init)
        offset = (i-1) * docp.discretization._step_variables_block
        xu[(offset+1):(offset+docp.dims.OCP_x)] .= x_init
    end
end

# +++ add here a more toplevel set_variables_at_time_step ? (potentially include stage variables later ?) and call it from docp
"""
$(TYPEDSIGNATURES)

Set initial guess for control variables at given time step
Convention: 1 <= i <= dim_NLP_steps(+1)
"""
function set_control_at_time_step!(xu, u_init, docp::DOCP, i)
    if !isnothing(u_init)
        if i <= docp.time.steps || (docp.discretization._final_control && i <= docp.time.steps + 1)
            offset = (i-1) * docp.discretization._step_variables_block + docp.dims.NLP_x
            xu[(offset+1):(offset+docp.dims.NLP_u)] .= u_init
        end
    end
end


"""
$(TYPEDSIGNATURES)

Set work array for all dynamics and lagrange cost evaluations
"""
function setWorkArray(docp::DOCP{<: Discretization}, xu, time_grid, v)
    return nothing
end


"""
$(TYPEDSIGNATURES)

Build sparsity pattern for Jacobian of constraints
"""
function DOCP_Jacobian_pattern(docp::DOCP{D}) where (D<:Discretization)
    error("DOCP_Jacobian_pattern not implemented for discretization ", D,
        " Use option solve(...; adnlp_backend=:optimized)")
end


"""
$(TYPEDSIGNATURES)

Build sparsity pattern for Hessian of Lagrangian
"""
function DOCP_Hessian_pattern(docp::DOCP{D}) where (D<:Discretization)
    error("DOCP_Hessian_pattern not implemented for discretization ", D,
        " Use option solve(...; adnlp_backend=:optimized)")
end


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
    return
end
function add_nonzero_block!(M, i, j; sym=false)
    M[i, j] = true
    sym && (M[j, i] = true)
    return
end
function add_nonzero_block!(Is, Js, i_start, i_end, j_start, j_end; sym=false)
    for i=i_start:i_end
        for j=j_start:j_end
            push!(Is, i)
            push!(Js, j)
            sym && push!(Is, j)
            sym && push!(Js, i)
        end
    end
    return
end
function add_nonzero_block!(Is, Js, i, j; sym=false)
    push!(Is, i)
    push!(Js, j)
    sym && push!(Is, j)
    sym && push!(Js, i)
    return
end

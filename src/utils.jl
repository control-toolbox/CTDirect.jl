"""
$(TYPEDSIGNATURES)

Retrieve optimization variables from the NLP variables.
Internal layout: [X0,U0, X1,U1, .., XN,UN,V]
"""
function get_variable(xu, docp)
    if docp.has_variable
        nx = docp.dim_NLP_x
        m = docp.dim_NLP_u
        N = docp.dim_NLP_steps
        offset = (nx+m) * (N+1)

        if docp.dim_NLP_v == 1
            return xu[offset + 1]
            #return xu[end]
        else
            #return xu[end-docp.dim_NLP_v+1:end]
            return xu[offset + 1: offset + docp.dim_NLP_v]
        end
    else
        return Float64[]
    end
end

"""
$(TYPEDSIGNATURES)

Retrieve a single optimization variable (no dim check).
Internal layout: [X0,U0, X1,U1, .., XN,UN,V]
"""
function get_single_variable(xu, docp, i::Int)
    if docp.has_variable
        #return xu[end-docp.dim_NLP_v+i]
        nx = docp.dim_NLP_x
        m = docp.dim_NLP_u
        N = docp.dim_NLP_steps
        offset = (nx+m) * (N+1)
        return xu[offset + i]
    else
        error("Tring to access variable in variable independent problem")
    end
end

"""
$(TYPEDSIGNATURES)

Retrieve state and control variables at given time step from the NLP variables. 
Internal layout: [X0,U0, X1,U1, .., XN,UN,V]
"""
# +++ this one would depend on the discretization scheme
function get_variables_at_time_step(xu, docp, i)

    nx = docp.dim_NLP_x
    n = docp.dim_OCP_x
    m = docp.dim_NLP_u
    offset = (nx + m) * i
    
    # retrieve scalar/vector OCP state (w/o lagrange state) 
    if n == 1
        xi = xu[offset+1]
    else
        xi = xu[offset+1:offset+n]
    end
    # NB. meaningful ONLY if has_lagrange is true !
    # add check without performance cost ?
    if docp.has_lagrange
        xli = xu[offset+nx]
    else
        xli = 0.
    end

    # retrieve scalar/vector control
    if m == 1
        ui = xu[offset+nx+1]
    else
        ui = xu[offset+nx+1:offset+nx+m]
    end

    return xi, ui, xli
end

# internal NLP version for solution parsing
function get_NLP_variables_at_time_step(xu, docp, i)

    nx = docp.dim_NLP_x
    m = docp.dim_NLP_u
    offset = (nx + m) * i
    
    xi = xu[offset+1:offset+nx]
    ui = xu[offset+nx+1:offset+nx+m]

    return xi, ui
end


#=
"""
$(TYPEDSIGNATURES)

Retrieve scalar/vector state variables at given time step from the NLP variables 
"""
function get_state_at_time_step(xu, docp, i)
    nx = docp.dim_NLP_x
    n = docp.dim_OCP_x
    if n == 1
        return xu[i*nx + 1]
    else
        return xu[i*nx + 1 : i*nx + n]
    end
end

"""
$(TYPEDSIGNATURES)

Retrieve the additional state variable corresponding to the lagrange (running) cost at given time step from the NLP variables
"""
function get_lagrange_cost_at_time_step(xu, docp, i)
    nx = docp.dim_NLP_x
    return xu[(i+1)*nx]
end

# internal NLP version for solution parsing
function get_NLP_state_at_time_step(xu, docp, i)
    nx = docp.dim_NLP_x
    return xu[i*nx + 1 : (i+1)*nx]
end


"""
$(TYPEDSIGNATURES)

Retrieve scalar/vector control variables at given time step from the NLP variables
"""
function get_control_at_time_step(xu, docp, i)
    nx = docp.dim_NLP_x
    m = docp.dim_NLP_u
    N = docp.dim_NLP_steps
    offset = (N+1)*nx

    if docp.dim_NLP_u == 1
        return xu[offset + i*m + 1]
    else
        return xu[offset + i*m + 1 : offset + (i+1)*m]
    end
end
=#

"""
$(TYPEDSIGNATURES)

Retrieve initial time for OCP (may be fixed or variable)
"""
function get_initial_time(xu, docp)
    if docp.has_free_t0
        return get_single_variable(xu, docp, Base.to_index(docp.ocp.initial_time))
    else
        return docp.ocp.initial_time
    end
end


"""
$(TYPEDSIGNATURES)

Retrieve final time for OCP (may be fixed or variable)
"""
function get_final_time(xu, docp)
    if docp.has_free_tf
        return get_single_variable(xu, docp, Base.to_index(docp.ocp.final_time))
    else
        return docp.ocp.final_time
    end
end


"""
$(TYPEDSIGNATURES)

Get actual (un-normalized) time value
"""
function get_unnormalized_time(xu, docp, t_normalized)
    t0 = get_initial_time(xu, docp)
    tf = get_final_time(xu, docp)    
    return t0 + t_normalized * (tf - t0)
end


"""
$(TYPEDSIGNATURES)

Get actual (un-normalized) time at give time step
"""
function get_time_at_time_step(xu, docp, i)
    return get_unnormalized_time(xu, docp, docp.NLP_normalized_time_grid[i+1])
end


#=
"""
$(TYPEDSIGNATURES)

Set state variables at given time step in the NLP variables (for initial guess)
"""
function set_state_at_time_step!(xu, x_init, docp, i)
    nx = docp.dim_NLP_x
    n = docp.dim_OCP_x
    # NB. only set first the actual state variables from the OCP (not the possible additional state for lagrange cost)
    xu[i*nx + 1 : i*nx + n] .= x_init
end


"""
$(TYPEDSIGNATURES)

Set control variables at given time step in the NLP variables (for initial guess)
"""
function set_control_at_time_step!(xu, u_init, docp, i)
    nx = docp.dim_NLP_x
    m = docp.dim_NLP_u
    N = docp.dim_NLP_steps
    offset = (N+1)*nx
    xu[offset + i*m + 1 : offset + i*m + m] .= u_init
end
=#

function set_variables_at_time_step!(xu, x_init, u_init, docp, i)
    nx = docp.dim_NLP_x
    n = docp.dim_OCP_x
    m = docp.dim_NLP_u
    N = docp.dim_NLP_steps
    offset = (nx + m) * i

    # NB. only set first the actual state variables from the OCP (not the possible additional state for lagrange cost)
    if !isnothing(x_init)
        xu[offset + 1 : offset + n] .= x_init
    end
    if !isnothing(u_init)
        xu[offset + nx + 1 : offset + nx + m] .= u_init
    end
end

"""
$(TYPEDSIGNATURES)

Set optimization variables in the NLP variables (for initial guess)
"""
function set_variable!(xu, v_init, docp)
    xu[end-docp.dim_NLP_v+1 : end] .= v_init
end


"""
$(TYPEDSIGNATURES)

Build initial guess for discretized problem
"""
function DOCP_initial_guess(docp,
    init::OptimalControlInit=OptimalControlInit())

    # default initialization
    # note: internal variables (lagrange cost, k_i for RK schemes) will keep these default values 
    NLP_X = 0.1 * ones(docp.dim_NLP_variables)

    # set variables if provided
    # (needed first in case of free times !)
    if !isnothing(init.variable_init)
        set_variable!(NLP_X, init.variable_init, docp)
    end

    # set state / control variables if provided
    for i in 0:docp.dim_NLP_steps
        ti = get_time_at_time_step(NLP_X, docp, i)
        set_variables_at_time_step!(NLP_X, init.state_init(ti), init.control_init(ti), docp, i)
    end

    return NLP_X
end

# placeholders (see CTDirectExt)
function export_ocp_solution end
function import_ocp_solution end


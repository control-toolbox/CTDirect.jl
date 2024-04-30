"""
$(TYPEDSIGNATURES)

Retrieve optimization variables from the NLP variables
"""
function get_variable(xu, docp)
    if docp.has_variable
        if docp.variable_dimension == 1
            return xu[end]
        else
            return xu[end-docp.variable_dimension+1:end]
        end
    else
        return Float64[]
    end
end


"""
$(TYPEDSIGNATURES)

Retrieve state variables at given time step from the NLP variables
"""
function get_state_at_time_step(xu, docp, i::Int64)
    """
        return
        x(t_i)
    """
    nx = docp.dim_NLP_state
    n = docp.ocp.state_dimension
    N = docp.dim_NLP_steps
    @assert i <= N "trying to get x(t_i) for i > N"
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
    nx = docp.dim_NLP_state
    N = docp.dim_NLP_steps
    @assert i <= N "trying to get lagrange cost at t_i for i > N"
    return xu[(i+1)*nx]
end

# internal vector version
function vget_state_at_time_step(xu, docp, i)
    nx = docp.dim_NLP_state
    N = docp.dim_NLP_steps
    @assert i <= N "trying to get x(t_i) for i > N"
    return xu[i*nx + 1 : (i+1)*nx]
end


"""
$(TYPEDSIGNATURES)

Retrieve control variables at given time step from the NLP variables
"""
function get_control_at_time_step(xu, docp, i)
    """
        return
        u(t_i)
    """
    nx = docp.dim_NLP_state
    m = docp.ocp.control_dimension
    N = docp.dim_NLP_steps
    @assert i <= N "trying to get u(t_i) for i > N"
    if m == 1
        return xu[(N+1)*nx + i*m + 1]
    else
        return xu[(N+1)*nx + i*m + 1 : (N+1)*nx + (i+1)*m]
    end
end

# internal vector version
function vget_control_at_time_step(xu, docp, i)
    nx = docp.dim_NLP_state
    m = docp.ocp.control_dimension
    N = docp.dim_NLP_steps
    @assert i <= N "trying to get u(t_i) for i > N"
    return xu[(N+1)*nx + i*m + 1 : (N+1)*nx + (i+1)*m]
end


"""
$(TYPEDSIGNATURES)

Retrieve initial time for OCP (may be fixed or variable)
"""
function get_initial_time(xu, docp)
    if docp.has_free_initial_time
        v = get_variable(xu, docp)
        return v[docp.ocp.initial_time]
    else
        return docp.ocp.initial_time
    end
end


"""
$(TYPEDSIGNATURES)

Retrieve final time for OCP (may be fixed or variable)
"""
function get_final_time(xu, docp)
    if docp.has_free_final_time
        v = get_variable(xu, docp)
        return v[docp.ocp.final_time]
    else
        return docp.ocp.final_time
    end
end


"""
$(TYPEDSIGNATURES)

Set state variables at given time step in the NLP variables (for initial guess)
"""
function set_state_at_time_step!(xu, x_init, docp, i)
    nx = docp.dim_NLP_state
    n = docp.ocp.state_dimension
    N = docp.dim_NLP_steps
    @assert i <= N "trying to set init for x(t_i) with i > N"
    # NB. only set first n components of state variable (ie the possible additional state for lagrange cost keeps the default init since it is not available from the OCP solution)
    if n == 1
        xu[i*n + 1] = x_init[]
    else
        xu[i*nx + 1 : i*nx + n] = x_init
    end
end


"""
$(TYPEDSIGNATURES)

Set control variables at given time step in the NLP variables (for initial guess)
"""
function set_control_at_time_step!(xu, u_init, docp, i)
    nx = docp.dim_NLP_state
    m = docp.ocp.control_dimension
    N = docp.dim_NLP_steps
    @assert i <= N "trying to set init for u(t_i) with i > N"
    offset = (N+1)*nx
    if m == 1
        xu[offset + i*m + 1] = u_init[]
    else        
        xu[offset + i*m + 1 : offset + i*m + m] = u_init
    end
end


"""
$(TYPEDSIGNATURES)

Set optimization variables in the NLP variables (for initial guess)
"""
function set_variable!(xu, v_init, docp)
    if docp.variable_dimension == 1
        xu[end] = v_init[]
    else
        xu[end-docp.variable_dimension+1 : end] = v_init
    end
end


"""
$(TYPEDSIGNATURES)

Build initial guess for discretized problem
"""
function initial_guess(docp, init::OCPInit=OCPInit())

    # default initialization
    # note: internal variables (lagrange cost, k_i for RK schemes) will keep these default values 
    xuv = 0.1 * ones(docp.dim_NLP_variables)

    # set variables if provided
    # (needed first in case of free times !)
    if !isnothing(init.variable_init)
        set_variable!(xuv, init.variable_init, docp)
    end

    # get time grid for state / control variables
    N = docp.dim_NLP_steps
    t0 = get_initial_time(xuv, docp)
    tf = get_final_time(xuv, docp)
    h = (tf - t0) / N 

    # set state / control variables if provided
    for i in 0:N
        ti = t0 + i * h
        if !isnothing(init.state_init(ti))
            set_state_at_time_step!(xuv, init.state_init(ti), docp, i)
        end
        if !isnothing(init.control_init(ti))
            set_control_at_time_step!(xuv, init.control_init(ti), docp, i)
        end
    end

    return xuv
end
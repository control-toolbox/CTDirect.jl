"""
$(TYPEDSIGNATURES)

Retrieve optimization variables from the NLP variables.
Internal layout: 
- Trapeze: [X_0,U_0, X_1,U_1, .., X_N,U_N,V]
- Midpoint: [X_0,U_0,K_0, X_1,U_1,K_1 .., X_N-1,U_N-1,K_N-1, XN, V]
with the convention u([t_i,t_i+1[) = U_i and u(tf) = U_N-1
"""
function get_optim_variable(xu, docp)

    if docp.has_variable
        #=nx = docp.dim_NLP_x
        m = docp.dim_NLP_u
        N = docp.dim_NLP_steps
        offset = (nx + m) * (N + 1)=#

        if docp.dim_NLP_v == 1
            #return xu[offset + 1]
            return xu[end]
        else
            #return xu[(offset + 1):(offset + docp.dim_NLP_v)]
            return xu[(end - docp.dim_NLP_v + 1):end]
        end
    else
        return similar(xu, 0) #Float64[]
    end
end


"""
$(TYPEDSIGNATURES)

Retrieve state and control variables at given time step from the NLP variables.
Internal layout: 
- Trapeze: [X_0,U_0, X_1,U_1, .., X_N,U_N,V]
- Midpoint: [X_0,U_0,K_0, X_1,U_1,K_1 .., X_N-1,U_N-1,K_N-1, XN, V]
with the convention u([t_i,t_i+1[) = U_i and u(tf) = U_N-1
"""
function get_variables_at_time_step(xu, docp, i)

    nx = docp.dim_NLP_x
    n = docp.dim_OCP_x
    m = docp.dim_NLP_u
    N = docp.dim_NLP_steps
    offset = (nx*2 + m) * i

    # retrieve scalar/vector OCP state (w/o lagrange state) 
    if n == 1
        xi = xu[offset + 1]
    else
        xi = xu[(offset + 1):(offset + n)]
    end
    # NB. meaningful ONLY if has_lagrange is true !
    # add check without performance cost ?
    if docp.has_lagrange
        xli = xu[offset + nx]
    else
        xli = 0.0
    end

    # retrieve scalar/vector control (convention u(tf) = U_N-1)
    if i < N
        if m == 1
            ui = xu[offset + nx + 1]
        else
            ui = xu[(offset + nx + 1):(offset + nx + m)]
        end
    else
        # final time: pick previous control
        offset2 = (nx*2 + m) * (i-1)
        if m == 1
            ui = xu[offset2 + nx + 1]
        else
            ui = xu[(offset2 + nx + 1):(offset2 + nx + m)]
        end
    end

    # retrieve vector stage variable (except at final time)
    if i < N
        ki = xu[(offset + nx + m + 1):(offset + nx + m + nx) ]
    else
        ki = nothing
    end

    return xi, ui, xli, ki
end


# internal NLP version for solution parsing
function get_NLP_variables_at_time_step(xu, docp, i)

    nx = docp.dim_NLP_x
    m = docp.dim_NLP_u
    N = docp.dim_NLP_steps
    offset = (nx*2 + m) * i

    # state
    xi = xu[(offset + 1):(offset + nx)]
    # control
    if i < N
        ui = xu[(offset + nx + 1):(offset + nx + m)]
    else
        offset2 = (nx*2 + m) * (i-1)
        ui = xu[(offset2 + nx + 1):(offset2 + nx + m)]
    end 
    # stage
    if i < N
        ki = xu[(offset + nx + m + 1):(offset + nx + m + nx) ]
    else
        ki = nothing
    end

    return xi, ui, ki
end


"""
$(TYPEDSIGNATURES)

Retrieve initial time for OCP (may be fixed or variable)
"""
function get_initial_time(xu, docp)
    if docp.has_free_t0
        #return get_single_variable(xu, docp, Base.to_index(docp.ocp.initial_time))
        #return get_optim_variable(xu, docp)[Base.to_index(docp.ocp.initial_time)]
        return get_optim_variable(xu, docp)[docp.ocp.initial_time]
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
        #return get_single_variable(xu, docp, Base.to_index(docp.ocp.final_time))
        return get_optim_variable(xu, docp)[docp.ocp.final_time]
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
    return get_unnormalized_time(xu, docp, docp.NLP_normalized_time_grid[i + 1])
end


function set_variables_at_time_step!(xu, x_init, u_init, docp, i)

    nx = docp.dim_NLP_x
    n = docp.dim_OCP_x
    m = docp.dim_NLP_u
    N = docp.dim_NLP_steps
    offset = (nx*2 + m) * i

    # NB. only set the actual state variables from the OCP 
    # - skip the possible additional state for lagrange cost
    # - skip internal discretization variables (K_i)
    if !isnothing(x_init)
        xu[(offset + 1):(offset + n)] .= x_init
    end
    if (i < N) && !isnothing(u_init)
        xu[(offset + nx + 1):(offset + nx + m)] .= u_init
    end
end

"""
$(TYPEDSIGNATURES)

Set optimization variables in the NLP variables (for initial guess)
"""
function set_optim_variable!(xu, v_init, docp)
    xu[(end - docp.dim_NLP_v + 1):end] .= v_init
end

"""
$(TYPEDSIGNATURES)

Build initial guess for discretized problem
"""
function DOCP_initial_guess(docp, init::OptimalControlInit = OptimalControlInit())

    # default initialization (internal variables such as lagrange cost, k_i for RK schemes) will keep these default values 
    NLP_X = 0.1 * ones(docp.dim_NLP_variables)

    # set variables if provided (needed first in case of free times !)
    if !isnothing(init.variable_init)
        set_variable!(NLP_X, init.variable_init, docp)
    end

    # set state / control variables if provided
    for i = 0:(docp.dim_NLP_steps)
        ti = get_time_at_time_step(NLP_X, docp, i)
        set_variables_at_time_step!(NLP_X, init.state_init(ti), init.control_init(ti), docp, i)
    end

    return NLP_X
end

# placeholders (see CTDirectExt)
function export_ocp_solution end
function import_ocp_solution end

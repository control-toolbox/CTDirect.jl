#= Functions for implicit midpoint discretization scheme
Internal layout for NLP variables: 
[X_0,U_0,K_0, X_1,U_1,K_1 .., X_N-1,U_N-1,K_N-1, XN, V]
with the convention u([t_i,t_i+1[) = U_i and u(tf) = U_N-1
=#


struct MidpointTag <: DiscretizationTag 
    stage::Int
    additional_controls::Int
    MidpointTag() = new(1, 0)
end


"""
$(TYPEDSIGNATURES)

Retrieve state and control variables at given time step from the NLP variables.
"""
function get_variables_at_time_step(xu, docp, i, tag::MidpointTag)

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
    if docp.has_lagrange
        xli = xu[offset + nx]
    else
        xli = nothing # dummy. use xu type ?
    end

    # retrieve scalar/vector control (convention u(tf) = U_N-1)
    if i < N
        offset_u = offset
    else
        offset_u = (nx*2 + m) * (i-1)
    end
    if m == 1
        ui = xu[offset_u + nx + 1]
    else
        ui = xu[(offset_u + nx + 1):(offset_u + nx + m)]
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
# could be fused with one above if 
# - using extended dynamics that include lagrange cost
# - scalar case is handled at OCP level
function get_NLP_variables_at_time_step(xu, docp, i, tag::MidpointTag)

    nx = docp.dim_NLP_x
    m = docp.dim_NLP_u
    N = docp.dim_NLP_steps
    offset = (nx*2 + m) * i

    # state
    xi = xu[(offset + 1):(offset + nx)]
    # control
    if i < N
        offset_u = offset
    else
        offset_u = (nx*2 + m) * (i-1)
    end 
    ui = xu[(offset_u + nx + 1):(offset_u + nx + m)]
    # stage
    if i < N
        ki = xu[(offset + nx + m + 1):(offset + nx + m + nx) ]
    else
        ki = nothing
    end

    return xi, ui, ki
end


function set_variables_at_time_step!(xu, x_init, u_init, docp, i, tag::MidpointTag)

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


# trivial version for now...
function initArgs(xu, docp, tag::MidpointTag)
    return xu, get_optim_variable(xu, docp)
end
function updateArgs(args, xu, v, docp, i, tag::MidpointTag)
    return args
end


"""
$(TYPEDSIGNATURES)

Set the constraints corresponding to the state equation
"""
function setStateEquation!(docp::DOCP, c, index::Int, xu, v, i, tag::MidpointTag)

    ocp = docp.ocp

    # variables
    ti = get_time_at_time_step(xu, docp, i)
    tip1 = get_time_at_time_step(xu, docp, i+1)
    hi = tip1 - ti

    xi, ui, xli, ki = get_variables_at_time_step(xu, docp, i, tag)
    xip1, uip1, xlip1 = get_variables_at_time_step(xu, docp, i+1, tag)
    
    v = get_optim_variable(xu, docp)

    # midpoint rule
    @. c[index:(index + docp.dim_OCP_x - 1)] =
        xip1 - (xi + hi * ki[1:docp.dim_OCP_x])
    # +++ just define extended dynamics !
    if docp.has_lagrange
        c[index + docp.dim_OCP_x] = xlip1 - (xli + hi * ki[end])
    end
    index = index + docp.dim_NLP_x
    
    # stage equation at mid-step
    t_s = 0.5 * (ti + tip1)
    x_s = 0.5 * (xi + xip1)
    c[index:(index + docp.dim_OCP_x - 1)] .=
        ki[1:docp.dim_OCP_x] .- ocp.dynamics(t_s, x_s, ui, v)
    # +++ just define extended dynamics !
    if docp.has_lagrange
        c[index + docp.dim_OCP_x] = ki[end] - ocp.lagrange(t_s, x_s, ui, v) 
    end
    index = index + docp.dim_NLP_x

    return index
end


"""
$(TYPEDSIGNATURES)

Set the path constraints at given time step
"""
function setPathConstraints!(docp::DOCP, c, index::Int, xu, v, i::Int, tag::MidpointTag)

    ocp = docp.ocp
    ti = get_time_at_time_step(xu, docp, i)
    xi, ui = get_variables_at_time_step(xu, docp, i, tag)

    # NB. using .= below *doubles* the allocations oO
    # pure control constraints
    if docp.dim_u_cons > 0
        c[index:(index + docp.dim_u_cons - 1)] = docp.control_constraints[2](ti, ui, v)
        index = index + docp.dim_u_cons
    end

    # pure state constraints
    if docp.dim_x_cons > 0
        c[index:(index + docp.dim_x_cons - 1)] = docp.state_constraints[2](ti, xi, v)
        index = index + docp.dim_x_cons
    end

    # mixed state / control constraints
    if docp.dim_mixed_cons > 0
        c[index:(index + docp.dim_mixed_cons - 1)] = docp.mixed_constraints[2](ti, xi, ui, v)
        index = index + docp.dim_mixed_cons
    end

    return index
end
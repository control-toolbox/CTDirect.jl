#= Functions for trapeze discretization scheme
Internal layout for NLP variables: 
[X_0,U_0, X_1,U_1, .., X_N,U_N, V]
=#

# NB. could be defined as a generic IRK
struct Trapeze <: Discretization

    stage::Int
    additional_controls::Int  # add control at tf
    info::String

    Trapeze() = new(0, 1, "Implicit Trapeze aka Crank-Nicolson, 2nd order, A-stable")
end


"""
$(TYPEDSIGNATURES)

Retrieve state and control variables at given time step from the NLP variables.
Convention: 1 <= i <= dim_NLP_steps+1
"""
#= split getters: same allocations overall...
function get_variables_at_time_step(xu, docp::DOCP{Trapeze}, i)

    nx = docp.dim_NLP_x
    n = docp.dim_OCP_x
    m = docp.dim_NLP_u
    offset = (nx + m) * (i-1)

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

    # retrieve scalar/vector control
    if m == 1
        ui = xu[offset + nx + 1]
    else
        ui = xu[(offset + nx + 1):(offset + nx + m)]
    end

    return xi, ui, xli
end
=#

function get_state_at_time_step(xu, docp::DOCP{Trapeze}, i)

    nx = docp.dim_NLP_x
    n = docp.dim_OCP_x
    m = docp.dim_NLP_u
    offset = (nx + m) * (i-1)

    # retrieve scalar/vector OCP state (w/o lagrange state) 
    if n == 1
        return xu[offset + 1]
    else
        return xu[(offset + 1):(offset + n)]
    end
end
#= using Val on dimension is much worse...
function get_state_at_time_step(xu, docp::DOCP{Trapeze}, ::Val{n}, i) where {n}

    nx = docp.dim_NLP_x
    #n = docp.dim_OCP_x
    m = docp.dim_NLP_u
    offset = (nx + m) * (i-1)

    # retrieve scalar/vector OCP state (w/o lagrange state) 
    if n == 1
        return xu[offset + 1]
    else
        return xu[(offset + 1):(offset + n)]
    end
end
=#

function get_lagrange_state_at_time_step(xu, docp::DOCP{Trapeze}, i)    

    nx = docp.dim_NLP_x
    n = docp.dim_OCP_x
    m = docp.dim_NLP_u
    offset = (nx + m) * (i-1)

    if docp.has_lagrange
        return xu[offset + nx]
    else
        error("problem has no lagrange cost")
    end
end

function get_control_at_time_step(xu, docp::DOCP{Trapeze}, i)

    nx = docp.dim_NLP_x
    n = docp.dim_OCP_x
    m = docp.dim_NLP_u
    offset = (nx + m) * (i-1)

    # retrieve scalar/vector control
    if m == 1
        return xu[offset + nx + 1]
    else
        return xu[(offset + nx + 1):(offset + nx + m)]
    end

end

# internal NLP version for solution parsing
# could be fused with one above if 
# - using extended dynamics that include lagrange cost
# - scalar case is handled at OCP level
function get_NLP_variables_at_t_i(xu, docp::DOCP{Trapeze}, i::Int)

    nx = docp.dim_NLP_x
    m = docp.dim_NLP_u
    offset = (nx + m) * i

    # state
    xi = xu[(offset + 1):(offset + nx)]
    # control
    ui = xu[(offset + nx + 1):(offset + nx + m)]

    return xi, ui
end


function set_variables_at_t_i!(xu, x_init, u_init, docp::DOCP{Trapeze}, i::Int)

    nx = docp.dim_NLP_x
    n = docp.dim_OCP_x
    m = docp.dim_NLP_u
    N = docp.dim_NLP_steps
    offset = (nx + m) * i

    # NB. only set the actual state variables from the OCP 
    # - skip the possible additional state for lagrange cost
    if !isnothing(x_init)
        xu[(offset + 1):(offset + n)] .= x_init
    end
    if !isnothing(u_init)
        xu[(offset + nx + 1):(offset + nx + m)] .= u_init
    end
end


# can later contain vectors for inplace getters ?
# for trapeze the dynamics (init at t0, compute and store at i+1)
function setWorkArray(docp::DOCP{Trapeze}, xu, time_grid, v)
   
    work = similar(xu, docp.dim_NLP_x)

    ocp = docp.ocp
    t0 = time_grid[1]
    x0 = get_state_at_time_step(xu, docp, 1)
    u0 = get_control_at_time_step(xu, docp, 1)

    if docp.has_inplace
        docp.dynamics((@view work[1:docp.dim_OCP_x]), t0, x0, u0, v)
    else
        work[1:docp.dim_OCP_x] .= docp.dynamics(t0, x0, u0, v)
    end
    
    if docp.has_lagrange
        if docp.has_inplace
            l = similar(xu, 1)
            docp.lagrange(l, t0, x0, u0, v) # +++cannot pass work[end] directly ?
            work[docp.dim_NLP_x] = l[1]
        else
            work[docp.dim_NLP_x] = docp.lagrange(t0, x0, u0, v)
        end
    end

    return work
end


"""
$(TYPEDSIGNATURES)

Set the constraints corresponding to the state equation
Convention: 1 <= i <= dim_NLP_steps
"""
function setConstraintBlock!(docp::DOCP{Trapeze}, c, xu, v, time_grid, i, work)

    # offset for previous steps
    offset = (i-1)*(docp.dim_NLP_x + docp.dim_path_cons)

    # variables
    ocp = docp.ocp
    ti = time_grid[i]
    xi = get_state_at_time_step(xu, docp, i)
    ui = get_control_at_time_step(xu, docp, i)
    fi = work[1:docp.dim_OCP_x] # copy !

    tip1 = time_grid[i+1]
    xip1 = get_state_at_time_step(xu, docp, i+1)
    uip1 = get_control_at_time_step(xu, docp, i+1)

    if docp.has_inplace
        docp.dynamics((@view work[1:docp.dim_OCP_x]), tip1, xip1, uip1, v)
    else
        work[1:docp.dim_OCP_x] .= docp.dynamics(tip1, xip1, uip1, v)
    end
    hi = tip1 - ti

    # trapeze rule with 'smart' update for dynamics
    @. c[offset+1:offset+docp.dim_OCP_x] = xip1 - (xi + 0.5 * hi * (fi + work[1:docp.dim_OCP_x]))
    
    if docp.has_lagrange
        xli = get_lagrange_state_at_time_step(xu, docp, i)
        xlip1 = get_lagrange_state_at_time_step(xu, docp, i+1)
        li = work[docp.dim_OCP_x+1]
        if docp.has_inplace
            lip1 = similar(xu, 1)
            docp.lagrange(lip1, tip1, xip1, uip1, v) # +++cannot pass work[end] directly ?
            work[docp.dim_OCP_x+1] = lip1[1]
        else
            work[docp.dim_OCP_x+1] = docp.lagrange(tip1, xip1, uip1, v)
        end
        c[offset+docp.dim_OCP_x+1] = xlip1 - (xli + 0.5 * hi * (li + work[docp.dim_OCP_x+1]))
    end
    offset += docp.dim_NLP_x    

    # path constraints
    setPathConstraints!(docp, c, ti, xi, ui, v, offset)

end

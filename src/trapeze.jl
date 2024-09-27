#= Functions for trapeze discretization scheme
Internal layout for NLP variables: 
[X_0,U_0, X_1,U_1, .., X_N,U_N, V]
=#

# NB. could be defined as a generic IRK
struct Trapeze <: Discretization

    stage::Int
    additional_controls::Int  # add control at tf
    info::String

    # constructor
    function Trapeze(dim_NLP_x, dim_NLP_u)
        return new(0, 1, "Implicit Trapeze aka Crank-Nicolson, 2nd order, A-stable")
    end
end

# not only for Trapeze, but may be redefined if needed
function get_OCP_variable_param(xu, docp::DOCP{<: Discretization, <: ScalVect, <: ScalVect, ScalVariable})
    return xu[end]
end
function get_OCP_variable_param(xu, docp::DOCP{<: Discretization, <: ScalVect, <: ScalVect, VectVariable})
    return @view xu[(end - docp.dim_NLP_v + 1):end]
end


"""
$(TYPEDSIGNATURES)

Retrieve state and control variables at given time step from the NLP variables.
Convention: 1 <= i <= dim_NLP_steps+1
"""
function get_state_at_time_step(xu, docp::DOCP{Trapeze}, i)
    offset = (i-1) * (docp.dim_NLP_x + docp.dim_NLP_u)
    return @view xu[(offset + 1):(offset + docp.dim_NLP_x)]
end
   
function get_OCP_state_at_time_step(xu, docp::DOCP{Trapeze}, i)
    offset = (i-1) * (docp.dim_NLP_x + docp.dim_NLP_u)
    if docp.dim_OCP_x == 1
        return xu[offset+1]
    else
        return @view xu[(offset + 1):(offset + docp.dim_OCP_x)]
    end
end

function get_OCP_state_at_time_step_param(xu, docp::DOCP{Trapeze, ScalVariable, <: ScalVect, <: ScalVect}, i)
    offset = (i-1) * (docp.dim_NLP_x + docp.dim_NLP_u)
    return xu[offset+1]
end
function get_OCP_state_at_time_step_param(xu, docp::DOCP{Trapeze, VectVariable, <: ScalVect, <: ScalVect}, i)
    offset = (i-1) * (docp.dim_NLP_x + docp.dim_NLP_u)
    return @view xu[(offset + 1):(offset + docp.dim_OCP_x)]
end

function get_lagrange_state_at_time_step(xu, docp, i)
    offset = (i-1) * (docp.dim_NLP_x + docp.dim_NLP_u)
    return xu[offset + docp.dim_NLP_x]
end

function get_control_at_time_step(xu, docp::DOCP{Trapeze}, i)
    offset = (i-1) * (docp.dim_NLP_x + docp.dim_NLP_u) + docp.dim_NLP_x
    return @view xu[(offset + 1):(offset + docp.dim_NLP_u)]
end

function get_OCP_control_at_time_step(xu, docp::DOCP{Trapeze}, i)
    offset = (i-1) * (docp.dim_NLP_x + docp.dim_NLP_u) + docp.dim_NLP_x
    if docp.dim_NLP_u == 1
        return xu[offset+1]
    else
        return @view xu[(offset + 1):(offset + docp.dim_NLP_u)]
    end
end

function get_OCP_control_at_time_step_param(xu, docp::DOCP{Trapeze, <: ScalVect, ScalVariable, <: ScalVect}, i)
    offset = (i-1) * (docp.dim_NLP_x + docp.dim_NLP_u) + docp.dim_NLP_x
    return xu[offset+1]
end
function get_OCP_control_at_time_step_param(xu, docp::DOCP{Trapeze, <: ScalVect, VectVariable, <: ScalVect}, i)
    offset = (i-1) * (docp.dim_NLP_x + docp.dim_NLP_u) + docp.dim_NLP_x
    return @view xu[(offset + 1):(offset + docp.dim_NLP_u)]
end


function set_state_at_time_step!(xu, x_init, docp::DOCP{Trapeze}, i)
    nx = docp.dim_NLP_x
    n = docp.dim_OCP_x
    m = docp.dim_NLP_u
    offset = (nx + m) * (i-1)
    # initialize only actual state variables from OCP (not lagrange state)
    if !isnothing(x_init)
        xu[(offset + 1):(offset + n)] .= x_init
    end
end

function set_control_at_time_step!(xu, u_init, docp::DOCP{Trapeze}, i)
    nx = docp.dim_NLP_x
    m = docp.dim_NLP_u
    offset = (nx + m) * (i-1)
    if !isnothing(u_init)
        xu[(offset + nx + 1):(offset + nx + m)] .= u_init
    end
end


function setWorkArray_param(docp::DOCP{Trapeze}, xu, time_grid, v)
    if docp.is_lagrange || !docp.is_inplace
        error("setworkarray to complete")
    end
    work = similar(xu, docp.dim_NLP_x)
    t0 = time_grid[1]
    x0 = get_OCP_state_at_time_step_param(xu, docp, 1)
    u0 = get_OCP_control_at_time_step_param(xu, docp, 1)
    docp.ocp.dynamics(work, t0, x0, u0, v)
    return work
end
function setWorkArray(docp::DOCP{Trapeze}, xu, time_grid, v)
    work = similar(xu, docp.dim_NLP_x)
    t0 = time_grid[1]
    x0 = get_state_at_time_step(xu, docp, 1)
    u0 = get_control_at_time_step(xu, docp, 1)
    docp.dynamics_ext!(work, t0, x0, u0, v)
    return work
end

"""
$(TYPEDSIGNATURES)

Set the constraints corresponding to the state equation
Convention: 1 <= i <= dim_NLP_steps (+1)
"""
function setConstraintBlock_param!(docp::DOCP{Trapeze}, c, xu, v, time_grid, i, work)
    if docp.is_lagrange || !docp.is_inplace
        error("setconstraintsblock to complete")
    end

    disc = docp.discretization

    # offset for previous steps
    offset = (i-1)*(docp.dim_NLP_x + docp.dim_path_cons)

    # 0. variables
    ti = time_grid[i]
    xi = get_OCP_state_at_time_step_param(xu, docp, i)
    ui = get_OCP_control_at_time_step_param(xu, docp, i)

    #1. state equation
    if i <= docp.dim_NLP_steps
        # more variables
        fi = copy(work) # create new copy, not just a reference
        tip1 = time_grid[i+1]
        xip1 = get_OCP_state_at_time_step_param(xu, docp, i+1)
        uip1 = get_OCP_control_at_time_step_param(xu, docp, i+1)
        docp.ocp.dynamics(work, tip1, xip1, uip1, v)

        # trapeze rule with 'smart' update for dynamics (need @. for scalar ?)
        c[offset+1:offset+docp.dim_NLP_x] = xip1 - (xi + 0.5 * (tip1 - ti) * (fi + work)) 
        offset += docp.dim_NLP_x
    end

    # 2. path constraints
    # Notes on allocations:.= seems similar
    if docp.dim_u_cons > 0
        if docp.is_inplace
            docp.control_constraints[2]((@view c[offset+1:offset+docp.dim_u_cons]),ti, ui, v)
        else
            c[offset+1:offset+docp.dim_u_cons] = docp.control_constraints[2](ti, ui, v)
        end
    end
    if docp.dim_x_cons > 0 
        if docp.is_inplace
            docp.state_constraints[2]((@view c[offset+docp.dim_u_cons+1:offset+docp.dim_u_cons+docp.dim_x_cons]),ti, xi, v)
        else
            c[offset+docp.dim_u_cons+1:offset+docp.dim_u_cons+docp.dim_x_cons] = docp.state_constraints[2](ti, xi, v)
        end
    end
    if docp.dim_mixed_cons > 0 
        if docp.is_inplace
            docp.mixed_constraints[2]((@view c[offset+docp.dim_u_cons+docp.dim_x_cons+1:offset+docp.dim_u_cons+docp.dim_x_cons+docp.dim_mixed_cons]), ti, xi, ui, v)
        else
            c[offset+docp.dim_u_cons+docp.dim_x_cons+1:offset+docp.dim_u_cons+docp.dim_x_cons+docp.dim_mixed_cons] = docp.mixed_constraints[2](ti, xi, ui, v)
        end
    end

end


function setConstraintBlock!(docp::DOCP{Trapeze}, c, xu, v, time_grid, i, work)

    disc = docp.discretization

    # offset for previous steps
    offset = (i-1)*(docp.dim_NLP_x + docp.dim_path_cons)

    # 0. variables
    ti = time_grid[i]
    xi = get_state_at_time_step(xu, docp, i)
    ui = get_control_at_time_step(xu, docp, i)

    #1. state equation
    if i <= docp.dim_NLP_steps
        # more variables
        fi = copy(work) # create new copy, not just a reference
        tip1 = time_grid[i+1]
        xip1 = get_state_at_time_step(xu, docp, i+1)
        uip1 = get_control_at_time_step(xu, docp, i+1)
        docp.dynamics_ext!(work, tip1, xip1, uip1, v)

        # trapeze rule with 'smart' update for dynamics (similar with @.)
        c[offset+1:offset+docp.dim_NLP_x] = xip1 - (xi + 0.5 * (tip1 - ti) * (fi + work)) #+++ jet runtime dispatch here even with explicit index ranges
        offset += docp.dim_NLP_x
    end

    # 2. path constraints
    # Notes on allocations:.= seems similar
    if docp.dim_u_cons > 0
        if docp.is_inplace
            docp.control_constraints[2]((@view c[offset+1:offset+docp.dim_u_cons]),ti, docp._u(ui), v)
        else
            c[offset+1:offset+docp.dim_u_cons] = docp.control_constraints[2](ti, docp._u(ui), v)
        end
    end
    if docp.dim_x_cons > 0 
        if docp.is_inplace
            docp.state_constraints[2]((@view c[offset+docp.dim_u_cons+1:offset+docp.dim_u_cons+docp.dim_x_cons]),ti, docp._x(xi), v)
        else
            c[offset+docp.dim_u_cons+1:offset+docp.dim_u_cons+docp.dim_x_cons] = docp.state_constraints[2](ti, docp._x(xi), v)
        end
    end
    if docp.dim_mixed_cons > 0 
        if docp.is_inplace
            docp.mixed_constraints[2]((@view c[offset+docp.dim_u_cons+docp.dim_x_cons+1:offset+docp.dim_u_cons+docp.dim_x_cons+docp.dim_mixed_cons]), ti, docp._x(xi), docp._u(ui), v)
        else
            c[offset+docp.dim_u_cons+docp.dim_x_cons+1:offset+docp.dim_u_cons+docp.dim_x_cons+docp.dim_mixed_cons] = docp.mixed_constraints[2](ti, docp._x(xi), docp._u(ui), v)
        end
    end

end

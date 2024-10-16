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
function get_OCP_variable(xu, docp::DOCP{<: Discretization, <: ScalVect, <: ScalVect, ScalVariable})
    return xu[docp.dim_NLP_variables]
end
function get_OCP_variable(xu, docp::DOCP{<: Discretization, <: ScalVect, <: ScalVect, VectVariable})
    return @view xu[(docp.dim_NLP_variables - docp.dim_NLP_v + 1):docp.dim_NLP_variables]
end


"""
$(TYPEDSIGNATURES)

Retrieve state and control variables at given time step from the NLP variables.
Convention: 1 <= i <= dim_NLP_steps+1
"""
function get_OCP_state_at_time_step(xu, docp::DOCP{Trapeze, ScalVariable, <: ScalVect, <: ScalVect}, i)
    offset = (i-1) * (docp.dim_NLP_x + docp.dim_NLP_u)
    return xu[offset+1]
end
function get_OCP_state_at_time_step(xu, docp::DOCP{Trapeze, VectVariable, <: ScalVect, <: ScalVect}, i)
    offset = (i-1) * (docp.dim_NLP_x + docp.dim_NLP_u)
    return @view xu[(offset + 1):(offset + docp.dim_OCP_x)]
end
function get_lagrange_state_at_time_step(xu, docp::DOCP{Trapeze}, i)
    offset = (i-1) * (docp.dim_NLP_x + docp.dim_NLP_u)
    return xu[offset + docp.dim_NLP_x]
end


function get_OCP_control_at_time_step(xu, docp::DOCP{Trapeze, <: ScalVect, ScalVariable, <: ScalVect}, i)
    offset = (i-1) * (docp.dim_NLP_x + docp.dim_NLP_u) + docp.dim_NLP_x
    return xu[offset+1]
end
function get_OCP_control_at_time_step(xu, docp::DOCP{Trapeze, <: ScalVect, VectVariable, <: ScalVect}, i)
    offset = (i-1) * (docp.dim_NLP_x + docp.dim_NLP_u) + docp.dim_NLP_x
    return @view xu[(offset + 1):(offset + docp.dim_NLP_u)]
end


function set_state_at_time_step!(xu, x_init, docp::DOCP{Trapeze}, i)
    offset = (i-1) * (docp.dim_NLP_x + docp.dim_NLP_u)
    # initialize only actual state variables from OCP (not lagrange state)
    if !isnothing(x_init)
        xu[(offset + 1):(offset + docp.dim_OCP_x)] .= x_init
    end
end

function set_control_at_time_step!(xu, u_init, docp::DOCP{Trapeze}, i)
    offset = (i-1) * (docp.dim_NLP_x + docp.dim_NLP_u) + docp.dim_NLP_x
    if !isnothing(u_init)
        xu[(offset + 1):(offset + docp.dim_NLP_u)] .= u_init
    end
end


function setWorkArray(docp::DOCP{Trapeze}, xu, time_grid, v)
    work = similar(xu, docp.dim_NLP_x)
    t0 = time_grid[1]
    x0 = get_OCP_state_at_time_step(xu, docp, 1)
    u0 = get_OCP_control_at_time_step(xu, docp, 1)
    if docp.is_inplace
        docp.ocp.dynamics((@view work[1:docp.dim_OCP_x]), t0, x0, u0, v)
    else
        work[1:docp.dim_OCP_x] .= docp.ocp.dynamics(t0, x0, u0, v)
    end
    if docp.is_lagrange
        if docp.is_inplace
            docp.ocp.lagrange((@view work[docp.dim_NLP_x:docp.dim_NLP_x]), t0, x0, u0, v)
        else
            work[docp.dim_NLP_x] = docp.ocp.lagrange(t0, x0, u0, v)
        end
    end
    return work
end


"""
$(TYPEDSIGNATURES)

Set the constraints corresponding to the state equation
Convention: 1 <= i <= dim_NLP_steps (+1)
"""
function setConstraintBlock!(docp::DOCP{Trapeze}, c, xu, v, time_grid, i, work)

    # offset for previous steps
    offset = (i-1)*(docp.dim_NLP_x + docp.dim_path_cons)

    # 0. variables
    ti = time_grid[i]
    xi = get_OCP_state_at_time_step(xu, docp, i)
    ui = get_OCP_control_at_time_step(xu, docp, i)

    # 1. state equation
    if i <= docp.dim_NLP_steps
        # more variables
        fi = copy(work) # create new copy, not just a reference
        tip1 = time_grid[i+1]
        xip1 = get_OCP_state_at_time_step(xu, docp, i+1)
        uip1 = get_OCP_control_at_time_step(xu, docp, i+1)
        if docp.is_inplace
            docp.ocp.dynamics((@view work[1:docp.dim_OCP_x]), tip1, xip1, uip1, v)
        else
            work[1:docp.dim_OCP_x] .= docp.ocp.dynamics(tip1, xip1, uip1, v)
        end
        if docp.is_lagrange
            if docp.is_inplace
                docp.ocp.lagrange((@view work[docp.dim_NLP_x:docp.dim_NLP_x]), tip1, xip1, uip1, v)
            else
                work[docp.dim_NLP_x] = docp.ocp.lagrange(tip1, xip1, uip1, v)
            end
        end

        # trapeze rule with 'smart' update for dynamics (@. allocs more and is slower -_-)
        # or split dyn and lag ?
        if docp.dim_OCP_x == 1
            c[offset+1] = xip1 - (xi + 0.5 * (tip1 - ti) * (fi[1] + work[1]))
        else
            c[offset+1:offset+docp.dim_OCP_x] = xip1 - (xi + 0.5 * (tip1 - ti) * (fi[1:docp.dim_OCP_x] + work[1:docp.dim_OCP_x]))
        end
        if docp.is_lagrange
            c[offset+docp.dim_NLP_x] = get_lagrange_state_at_time_step(xu, docp, i+1) - (get_lagrange_state_at_time_step(xu, docp, i) + 0.5 * (tip1 - ti) * (fi[docp.dim_NLP_x] + work[docp.dim_NLP_x]))
        end
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

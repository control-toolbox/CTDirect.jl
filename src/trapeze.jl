#= Functions for trapeze discretization scheme
Internal layout for NLP variables: 
[X_0,U_0, X_1,U_1, .., X_N,U_N, V]
=#
#using StaticArrays

# NB. could be defined as a generic IRK
struct Trapeze <: Discretization

    stage::Int
    additional_controls::Int  # add control at tf
    info::String

    #= work arrays (worse allocs even with inplace getters...)
    work::AbstractVector{<:Real}
    xi::AbstractVector{<:Real}
    xip1::AbstractVector{<:Real}=#

    # we could pass some ocp field to constructor (eg dims, lagranage)
    Trapeze() = new(0, 1, "Implicit Trapeze aka Crank-Nicolson, 2nd order, A-stable")
end

#=function initWork(discretization::Trapeze, dim_NLP_x::Int)
    resize!(discretization.work, dim_NLP_x)
    resize!(discretization.xi, dim_NLP_x)
    resize!(discretization.xip1, dim_NLP_x)
end=#

"""
$(TYPEDSIGNATURES)

Retrieve state and control variables at given time step from the NLP variables.
Convention: 1 <= i <= dim_NLP_steps+1
"""
function get_state_at_time_step(xu, docp::DOCP{Trapeze}, i)
    nx = docp.dim_NLP_x
    m = docp.dim_NLP_u
    offset = (nx + m) * (i-1)
    return @view xu[(offset + 1):(offset + nx)]
end

#=function get_state_at_time_step!(x, xu, docp::DOCP{Trapeze}, i)
    nx = docp.dim_NLP_x
    m = docp.dim_NLP_u
    offset = (nx + m) * (i-1)
    x[1:nx] = @view xu[(offset + 1):(offset + nx)] #view ?
    return
end=#

function get_control_at_time_step(xu, docp::DOCP{Trapeze}, i)
    nx = docp.dim_NLP_x
    m = docp.dim_NLP_u
    offset = (nx + m) * (i-1)
    return @view xu[(offset + nx + 1):(offset + nx + m)]
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


function setWorkArray(docp::DOCP{Trapeze}, xu, time_grid, v)
   
    work = similar(xu, docp.dim_NLP_x)
    t0 = time_grid[1]
    x0 = get_state_at_time_step(xu, docp, 1)
    u0 = get_control_at_time_step(xu, docp, 1)

    if docp.has_inplace
        docp.dynamics_ext(work, t0, x0, u0, v)
    else
        # NB. work = will create a new variable ;-) (work .= is fine)
        work[:] = docp.dynamics_ext(t0, x0, u0, v)
    end
    return work
end

#=function setWorkArray(docp::DOCP{Trapeze}, xu, time_grid, v)

    t0 = time_grid[1]
    x0 = get_state_at_time_step(xu, docp, 1)
    u0 = get_control_at_time_step(xu, docp, 1)

    if docp.has_inplace
        docp.dynamics_ext(docp.discretization.work, t0, x0, u0, v)
    else
        # NB. work = will create a new variable ;-) (work .= is fine)
        docp.discretization.work[:] = docp.dynamics_ext(t0, x0, u0, v)
    end
    return
end=#

"""
$(TYPEDSIGNATURES)

Set the constraints corresponding to the state equation
Convention: 1 <= i <= dim_NLP_steps (+1)
"""
function setConstraintBlock!(docp::DOCP{Trapeze}, c, xu, v, time_grid, i, work)
#function setConstraintBlock!(docp::DOCP{Trapeze}, c, xu, v, time_grid, i)

    # offset for previous steps
    offset = (i-1)*(docp.dim_NLP_x + docp.dim_path_cons)

    #=trying to use work array in discretization struct seems much worse...
    work = docp.discretization.work
    xi = docp.discretization.xi
    xip1 = docp.discretization.xip1=#

    # 0. variables
    ti = time_grid[i]
    xi = get_state_at_time_step(xu, docp, i)
    #get_state_at_time_step!(xi, xu, docp, i)
    ui = get_control_at_time_step(xu, docp, i)

    #1. state equation
    if i <= docp.dim_NLP_steps
        # more variables
        fi = copy(work) # create new copy, not just a reference
        tip1 = time_grid[i+1]
        xip1 = get_state_at_time_step(xu, docp, i+1)
        #get_state_at_time_step!(xip1, xu, docp, i+1)
        uip1 = get_control_at_time_step(xu, docp, i+1)
        if docp.has_inplace
            docp.dynamics_ext(work, tip1, xip1, uip1, v)
        else
            # copy, do not create a new variable !
            work[:] = docp.dynamics_ext(tip1, xip1, uip1, v)
        end

        # trapeze rule with 'smart' update for dynamics (similar with @.)
        c[offset+1:offset+docp.dim_NLP_x] = xip1 - (xi + 0.5 * (tip1 - ti) * (fi + work))
        offset += docp.dim_NLP_x
    end

    # 2. path constraints
    #setPathConstraints!(docp, c, ti, xi, ui, v, offset)
    #Apparently function calls always seem to add some overhead allocation...

    # Notes on allocations:.= seems similar
    if docp.dim_u_cons > 0
        if docp.has_inplace
            docp.control_constraints[2]((@view c[offset+1:offset+docp.dim_u_cons]),ti, docp._u(ui), v)
        else
            c[offset+1:offset+docp.dim_u_cons] = docp.control_constraints[2](ti, docp._u(ui), v)
        end
    end
    if docp.dim_x_cons > 0 
        if docp.has_inplace
            docp.state_constraints[2]((@view c[offset+docp.dim_u_cons+1:offset+docp.dim_u_cons+docp.dim_x_cons]),ti, docp._x(xi), v)
        else
            c[offset+docp.dim_u_cons+1:offset+docp.dim_u_cons+docp.dim_x_cons] = docp.state_constraints[2](ti, docp._x(xi), v)
        end
    end
    if docp.dim_mixed_cons > 0 
        if docp.has_inplace
            docp.mixed_constraints[2]((@view c[offset+docp.dim_u_cons+docp.dim_x_cons+1:offset+docp.dim_u_cons+docp.dim_x_cons+docp.dim_mixed_cons]), ti, docp._x(xi), docp._u(ui), v)
        else
            c[offset+docp.dim_u_cons+docp.dim_x_cons+1:offset+docp.dim_u_cons+docp.dim_x_cons+docp.dim_mixed_cons] = docp.mixed_constraints[2](ti, docp._x(xi), docp._u(ui), v)
        end
    end

end

#= Functions for trapeze discretization scheme
Internal layout for NLP variables: 
[X_0,U_0, X_1,U_1, .., X_N,U_N, V]
=#

# NB. could be defined as a generic IRK
struct Trapeze <: Discretization

    stage::Int
    additional_controls::Int  # add control at tf
    info::String

    # work arrays (x, u etc also ?)
    work::AbstractVector

    # we could pass some ocp field to constructor (eg dims, lagranage)
    Trapeze() = new(0, 1, "Implicit Trapeze aka Crank-Nicolson, 2nd order, A-stable", [])
end

function initWork(discretization::Trapeze, dim_NLP_x::Int)
    resize!(discretization.work, dim_NLP_x)
end

"""
$(TYPEDSIGNATURES)

Retrieve state and control variables at given time step from the NLP variables.
Convention: 1 <= i <= dim_NLP_steps+1
"""
function vget_state_at_time_step(xu, docp::DOCP{Trapeze}, i)

    nx = docp.dim_NLP_x
    m = docp.dim_NLP_u
    offset = (nx + m) * (i-1)

    return @view xu[(offset + 1):(offset + nx)]
end

function vget_control_at_time_step(xu, docp::DOCP{Trapeze}, i)

    nx = docp.dim_NLP_x
    m = docp.dim_NLP_u
    offset = (nx + m) * (i-1)

    return @view xu[(offset + nx + 1):(offset + nx + m)]

end

# internal NLP version for solution parsing
# could be fused with one above if 
# - using extended dynamics that include lagrange cost
# - scalar case is handled at OCP level
# +++ remove
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

# +++ split
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
    t0 = time_grid[1]
    x0 = vget_state_at_time_step(xu, docp, 1)
    u0 = vget_control_at_time_step(xu, docp, 1)

    if docp.has_inplace
        # passing just work fails some tests oO ?
        #docp.dynamics_ext((@view work[1:docp.dim_NLP_x]), t0, x0, u0, v)
        docp.dynamics_ext(work, t0, x0, u0, v)
    else
        # NB. work = will create a new variable ;-) (work .= is fine)
        work[:] = docp.dynamics_ext(t0, x0, u0, v)
    end
    
    return work
end


"""
$(TYPEDSIGNATURES)

Set the constraints corresponding to the state equation
Convention: 1 <= i <= dim_NLP_steps
"""
function setConstraintBlock!(docp::DOCP{Trapeze}, c, xu, v, time_grid, i, work)

    # NB. the handling of the work array is a bit tricky
    # some formulations of the copy/affectation give incorrect results...

    # offset for previous steps
    offset = (i-1)*(docp.dim_NLP_x + docp.dim_path_cons)

    # variables
    ti = time_grid[i]
    xi = vget_state_at_time_step(xu, docp, i)
    ui = vget_control_at_time_step(xu, docp, i)
    fi = copy(work) # create new copy, not just a reference

    tip1 = time_grid[i+1]
    xip1 = vget_state_at_time_step(xu, docp, i+1)
    uip1 = vget_control_at_time_step(xu, docp, i+1)

    if docp.has_inplace
        # passing just work fails some tests oO...
        #docp.dynamics_ext((@view work[1:docp.dim_NLP_x]), tip1, xip1, uip1, v)
        docp.dynamics_ext(work, tip1, xip1, uip1, v)
    else
        # copy, do not create a new variable !
        work[:] = docp.dynamics_ext(tip1, xip1, uip1, v)
    end

    # trapeze rule with 'smart' update for dynamics (similar with @.)
    c[offset+1:offset+docp.dim_NLP_x] = xip1 - (xi + 0.5 * (tip1 - ti) * (fi + work))
    offset += docp.dim_NLP_x

    # path constraints
    setPathConstraints!(docp, c, ti, xi, ui, v, offset)

end

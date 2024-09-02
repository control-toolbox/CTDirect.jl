#= Functions for trapeze discretization scheme
Internal layout for NLP variables: 
[X_0,U_0, X_1,U_1, .., X_N,U_N, V]
=#


struct TrapezeTag <: DiscretizationTag 
    stage::Int
    additional_controls::Int
    # add control at tf
    TrapezeTag() = new(0, 1)
end


"""
$(TYPEDSIGNATURES)

Retrieve state and control variables at given time step from the NLP variables. 
"""
function get_variables_at_time_step(xu, docp::DOCP{TrapezeTag}, i)

    nx = docp.dim_NLP_x
    n = docp.dim_OCP_x
    m = docp.dim_NLP_u
    offset = (nx + m) * i

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


# internal NLP version for solution parsing
# could be fused with one above if 
# - using extended dynamics that include lagrange cost
# - scalar case is handled at OCP level
function get_NLP_variables_at_time_step(xu, docp, i, tag::TrapezeTag)

    nx = docp.dim_NLP_x
    m = docp.dim_NLP_u
    offset = (nx + m) * i

    # state
    xi = xu[(offset + 1):(offset + nx)]
    # control
    ui = xu[(offset + nx + 1):(offset + nx + m)]

    return xi, ui
end


function set_variables_at_time_step!(xu, x_init, u_init, docp, i, tag::TrapezeTag)

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


# ? use abstract type for args ?
"""
$(TYPEDSIGNATURES)

Useful values at a time step: time, state, control, dynamics...
"""
struct ArgsAtTimeStep_Trapeze
    time::Any
    state::Any
    control::Any
    dynamics::Any
    lagrange_state::Any
    lagrange_cost::Any

    function ArgsAtTimeStep_Trapeze(xu, docp::DOCP{TrapezeTag}, v, time_grid, i::Int)

        # variables
        ti = time_grid[i+1]
        xi, ui, xli = get_variables_at_time_step(xu, docp, i)

        # dynamics and lagrange cost
        fi = docp.ocp.dynamics(ti, xi, ui, v)

        if docp.has_lagrange
            li = docp.ocp.lagrange(ti, xi, ui, v)
            args = new(ti, xi, ui, fi, xli, li)
        else
            args = new(ti, xi, ui, fi)
        end

        return args
    end
end
# +++multiple dispatch here seems to cause more allocations !
function initArgs(xu, docp::DOCP{TrapezeTag}, time_grid)
    # optimization variables
    v = Float64[]
    docp.has_variable && (v = get_optim_variable(xu, docp))
    args_i = ArgsAtTimeStep_Trapeze(xu, docp, v, time_grid, 0)
    args_ip1 = ArgsAtTimeStep_Trapeze(xu, docp, v, time_grid, 1)
    return (args_i, args_ip1), v
end
function updateArgs(args, xu, docp::DOCP{TrapezeTag}, v, time_grid, i)
    args_i, args_ip1 = args
    if i < docp.dim_NLP_steps - 1
        # are we allocating more than one args here ?
        return (args_ip1, ArgsAtTimeStep_Trapeze(xu, docp, v, time_grid, i+2))
    else
        return (args_ip1, args_ip1)
    end
end


"""
$(TYPEDSIGNATURES)

Set the constraints corresponding to the state equation
"""
function setStateEquation!(docp::DOCP{TrapezeTag}, c, index::Int, args, v, i)

    # NB. arguments v,i are unused here but present for unified call
    ocp = docp.ocp
    args_i, args_ip1 = args
    hi = args_ip1.time - args_i.time

    # trapeze rule (NB. @. allocates more ...)
    c[index:(index + docp.dim_OCP_x - 1)] .=
        args_ip1.state .- (args_i.state .+ 0.5 * hi * (args_i.dynamics .+ args_ip1.dynamics))
    # +++ just define extended dynamics !
    if docp.has_lagrange
        c[index + docp.dim_OCP_x] =
            args_ip1.lagrange_state -
            (args_i.lagrange_state + 0.5 * hi * (args_i.lagrange_cost + args_ip1.lagrange_cost))
    end
    index += docp.dim_NLP_x

    return index
end


"""
$(TYPEDSIGNATURES)

Set the path constraints at given time step
"""
function setPathConstraints!(docp::DOCP{TrapezeTag}, c, index::Int, args, v, i::Int)

    # note: i is unused but passed for call compatibility
    ocp = docp.ocp
    args_i, args_ip1 = args
    ti = args_i.time
    xi = args_i.state
    ui = args_i.control

    # NB. using .= below *doubles* the allocations oO
    # pure control constraints
    if docp.dim_u_cons > 0
        c[index:(index + docp.dim_u_cons - 1)] = docp.control_constraints[2](ti, ui, v)
        index += docp.dim_u_cons
    end

    # pure state constraints
    if docp.dim_x_cons > 0
        c[index:(index + docp.dim_x_cons - 1)] = docp.state_constraints[2](ti, xi, v)
        index += docp.dim_x_cons
    end

    # mixed state / control constraints
    if docp.dim_mixed_cons > 0
        c[index:(index + docp.dim_mixed_cons - 1)] = docp.mixed_constraints[2](ti, xi, ui, v)
        index += docp.dim_mixed_cons
    end

    return index
end

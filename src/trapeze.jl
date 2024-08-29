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
function get_variables_at_time_step(xu, docp, i, tag::TrapezeTag)

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


function initArgs(xu, docp, tag::TrapezeTag)
    v = get_optim_variable(xu, docp)
    args_i = ArgsAtTimeStep(xu, docp, 0, v)
    args_ip1 = ArgsAtTimeStep(xu, docp, 1, v)
    return (args_i, args_ip1), v
end


function updateArgs((args_i, args_ip1), xu, v, docp, i, tag::TrapezeTag)
    if i < docp.dim_NLP_steps - 1
        # are we allocating more than one args here ?
        return (args_ip1, ArgsAtTimeStep(xu, docp, i + 2, v))
    else
        return (args_ip1, args_ip1)
    end
end


"""
$(TYPEDSIGNATURES)

Useful values at a time step: time, state, control, dynamics...
"""
struct ArgsAtTimeStep
    time::Any
    state::Any
    control::Any
    dynamics::Any
    lagrange_state::Any
    lagrange_cost::Any

    function ArgsAtTimeStep(xu, docp::DOCP, i::Int, v)

        # variables
        ti = get_time_at_time_step(xu, docp, i)
        xi, ui, xli = get_variables_at_time_step(xu, docp, i, docp.discretization)

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


"""
$(TYPEDSIGNATURES)

Set the constraints corresponding to the state equation
"""
function setStateEquation!(docp::DOCP, c, index::Int, (args_i, args_ip1), v, i, tag::TrapezeTag)

    # NB. arguments v,i are unused here but present for unified call
    ocp = docp.ocp
    hi = args_ip1.time - args_i.time

    # trapeze rule
    c[index:(index + docp.dim_OCP_x - 1)] .=
        args_ip1.state .- (args_i.state .+ 0.5 * hi * (args_i.dynamics .+ args_ip1.dynamics))
    # +++ just define extended dynamics !
    if docp.has_lagrange
        c[index + docp.dim_OCP_x] =
            args_ip1.lagrange_state -
            (args_i.lagrange_state + 0.5 * hi * (args_i.lagrange_cost + args_ip1.lagrange_cost))
    end
    index = index + docp.dim_NLP_x

    return index
end


"""
$(TYPEDSIGNATURES)

Set the path constraints at given time step
"""
function setPathConstraints!(docp::DOCP, c, index::Int, (args_i, args_ip1), v, tag::TrapezeTag)

    ocp = docp.ocp
    ti = args_i.time
    xi = args_i.state
    ui = args_i.control

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

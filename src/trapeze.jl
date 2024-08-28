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


"""
$(TYPEDSIGNATURES)

Set the constraints corresponding to the state equation
"""
function setStateEquation!(docp::DOCP, c, index::Int, args_trapeze, tag::TrapezeTag)
    ocp = docp.ocp
    args_i = args_trapeze[1]
    args_ip1 = args_trapeze[2]
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

#= Functions for implicit midpoint discretization scheme
Internal layout for NLP variables: 
[X_0,U_0,K_0, X_1,U_1,K_1 .., X_N-1,U_N-1,K_N-1, XN, V]
with the convention u([t_i,t_i+1[) = U_i and u(tf) = U_N-1
=#

# +++ TODO: use args
# NB. could be defined as a generic IRK
struct Midpoint <: Discretization

    stage::Int
    additional_controls::Int
    butcher_a::Matrix{Float64}
    butcher_b::Vector{Float64}
    butcher_c::Vector{Float64}
    info::String

    Midpoint() = new(1, 0, hcat(0.5), [1], [0.5], "Implicit Midpoint aka Gauss-Legendre collocation for s=1, 2nd order, symplectic")
end


"""
$(TYPEDSIGNATURES)

Retrieve state and control variables at given time step from the NLP variables.
"""
function get_variables_at_t_i(xu, docp::DOCP{Midpoint}, i)

    nx = docp.dim_NLP_x
    n = docp.dim_OCP_x
    m = docp.dim_NLP_u
    N = docp.dim_NLP_steps
    offset = (nx*(1+docp.discretization.stage) + m) * i

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
        offset_u = (nx * (1 + docp.discretization.stage) + m) * (i - 1)
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
function get_NLP_variables_at_t_i(xu, docp::DOCP{Midpoint}, i)

    nx = docp.dim_NLP_x
    m = docp.dim_NLP_u
    N = docp.dim_NLP_steps
    offset = (nx * (1 + docp.discretization.stage) + m) * i

    # state
    xi = xu[(offset + 1):(offset + nx)]
    # control
    if i < N
        offset_u = offset
    else
        offset_u = (nx * (1 + docp.discretization.stage) + m) * (i - 1)
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


function set_variables_at_t_i!(xu, x_init, u_init, docp::DOCP{Midpoint}, i)

    nx = docp.dim_NLP_x
    n = docp.dim_OCP_x
    m = docp.dim_NLP_u
    N = docp.dim_NLP_steps
    offset = (nx * (1 + docp.discretization.stage) + m) * i

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

Useful values at a time step: time, state, control, dynamics...
"""
struct Midpoint_Args <: ArgsAtStep
    variable
    time
    state
    control
    stage_k
    next_time
    next_state
    lagrange_state
    next_lagrange_state
end

function initArgs(docp::DOCP{Midpoint}, xu)

    args = Vector{Midpoint_Args}(undef, docp.dim_NLP_steps + 1)
    dummy = similar(xu,0)

    # get time grid
    time_grid = get_time_grid(xu, docp)

    # get optim variable
    if docp.has_variable
        v = get_optim_variable(xu, docp)
    else
        v = dummy
    end

    # loop over time steps
    for i = 1:docp.dim_NLP_steps
        t_i = time_grid[i]
        t_ip1 = time_grid[i+1]
        x_i, u_i, xl_i, k_i = get_variables_at_t_i(xu, docp, i-1)
        x_ip1, u_ip1, xl_ip1 = get_variables_at_t_i(xu, docp, i)
        args[i] = Midpoint_Args(v, t_i, x_i, u_i, k_i, t_ip1, x_ip1, xl_i, xl_ip1)
    end

    # final time: for path constraints only
    # useful fields are: v, ti, xi, ui
    t_f = time_grid[docp.dim_NLP_steps+1]
    x_f, u_f = get_variables_at_t_i(xu, docp, docp.dim_NLP_steps)
    args[docp.dim_NLP_steps+1] = Midpoint_Args(v, t_f, x_f, u_f, dummy, dummy, dummy, dummy, dummy)

    return args

end


"""
$(TYPEDSIGNATURES)

Set the constraints corresponding to the state equation
"""
function setStateEquation!(docp::DOCP{Midpoint}, c_block, args::Midpoint_Args)

    ocp = docp.ocp
    disc = docp.discretization

    # variables
    v = args.variable
    ti = args.time
    xi = args.state
    ui = args.control
    xli = args.lagrange_state
    ki = args.stage_k
    tip1 = args.next_time
    xip1 = args.next_state
    xlip1 = args.next_lagrange_state
    hi = tip1 - ti

    # midpoint rule
    h_sum_bk = hi * disc.butcher_b[1] * ki[1:docp.dim_NLP_x]
    c_block[1:docp.dim_OCP_x] .= xip1 .- (xi .+ h_sum_bk[1:docp.dim_OCP_x])
    # +++ just define extended dynamics !
    docp.has_lagrange && (c_block[docp.dim_NLP_x] = xlip1 - (xli + h_sum_bk[end]))

    # stage equation at mid-step
    t_s = ti + hi * disc.butcher_c[1]
    if docp.dim_OCP_x == 1
        x_s = xi + hi * disc.butcher_a[1][1] * ki[1] #FFS
    else
        x_s = xi .+ hi .* (disc.butcher_a[1][1] .* ki[1:docp.dim_OCP_x])
    end
    c_block[docp.dim_NLP_x+1:docp.dim_NLP_x+docp.dim_OCP_x] .= ki[1:docp.dim_OCP_x] .- ocp.dynamics(t_s, x_s, ui, v)
    # +++ just define extended dynamics !
    docp.has_lagrange && (c_block[docp.dim_NLP_x*2] = ki[end] - ocp.lagrange(t_s, x_s, ui, v))

end

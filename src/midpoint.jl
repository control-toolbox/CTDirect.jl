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
Convention: 1 <= i <= dim_NLP_steps+1
"""
function get_variables_at_time_step(xu, docp::DOCP{Midpoint}, i)

    nx = docp.dim_NLP_x
    n = docp.dim_OCP_x
    m = docp.dim_NLP_u
    N = docp.dim_NLP_steps
    offset = (nx*(1+docp.discretization.stage) + m) * (i-1)

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
    if i < N+1
        offset_u = offset
    else
        offset_u = (nx*(1+docp.discretization.stage) + m) * (i-2)
    end
    if m == 1
        ui = xu[offset_u + nx + 1]
    else
        ui = xu[(offset_u + nx + 1):(offset_u + nx + m)]
    end

    # retrieve vector stage variable (except at final time)
    if i < N+1
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
    offset = (nx*(1+docp.discretization.stage) + m) * i

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
        ki = xu[(offset + nx + m + 1):(offset + nx + m + nx)]
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
    offset = (nx*(1+docp.discretization.stage) + m) * i

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


# can later contain vectors for inplace getters ?
# for trapeze the dynamics (init at t0, compute and store at i+1)
function setWorkArray(docp::DOCP{Midpoint}, xu, time_grid, v)
    return nothing
end


"""
$(TYPEDSIGNATURES)

Set the constraints corresponding to the state equation
Convention: 1 <= i <= dim_NLP_steps
"""
function setConstraintBlock!(docp::DOCP{Midpoint}, c, xu, v, time_grid, i, work)

    # offset for previous steps
    offset = (i-1)*(docp.dim_NLP_x * (1+docp.discretization.stage) + docp.dim_path_cons)

    # variables
    ocp = docp.ocp
    disc = docp.discretization
    ti = time_grid[i]
    xi, ui, xli, ki = get_variables_at_time_step(xu, docp, i)

    tip1 = time_grid[i+1]
    xip1, _, xlip1 = get_variables_at_time_step(xu, docp, i+1)
    
    hi = tip1 - ti

    # midpoint rule
    h_sum_bk = hi * disc.butcher_b[1] * ki[1:docp.dim_NLP_x]
    c[offset+1:offset+docp.dim_OCP_x] .= xip1 .- (xi .+ h_sum_bk[1:docp.dim_OCP_x])

    if docp.has_lagrange
        c[offset+docp.dim_NLP_x] = xlip1 - (xli + h_sum_bk[end])
    end
    offset += docp.dim_NLP_x

    # stage equation at mid-step
    t_s = ti + hi * disc.butcher_c[1]
    if docp.dim_OCP_x == 1
        x_s = xi + hi * disc.butcher_a[1][1] * ki[1] #FFS
    else
        x_s = xi .+ hi .* (disc.butcher_a[1][1] .* ki[1:docp.dim_OCP_x])
    end
    c[offset+1:offset+docp.dim_OCP_x] .= ki[1:docp.dim_OCP_x] .- ocp.dynamics(t_s, x_s, ui, v)
    # +++ just define extended dynamics !
    if docp.has_lagrange
        c[offset+docp.dim_NLP_x] = ki[end] - ocp.lagrange(t_s, x_s, ui, v)
    end
    offset += docp.dim_NLP_x

    # path constraints
    setPathConstraints!(docp, c, ti, xi, ui, v, offset)

end

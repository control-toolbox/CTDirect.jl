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
<<<<<<< HEAD

    # constructor    
    function Midpoint(dim_NLP_x, dim_NLP_u, dim_NLP_steps) 
        return new(1, 0, hcat(0.5), [1], [0.5], "Implicit Midpoint aka Gauss-Legendre collocation for s=1, 2nd order, symplectic")
    end
end
=======
>>>>>>> main

    Midpoint() = new(1, 0, hcat(0.5), [1], [0.5], "Implicit Midpoint aka Gauss-Legendre collocation for s=1, 2nd order, symplectic")
end

"""
$(TYPEDSIGNATURES)

Retrieve state and control variables at given time step from the NLP variables.
Convention: 1 <= i <= dim_NLP_steps+1
"""
<<<<<<< HEAD
function get_state_at_time_step(xu, docp::DOCP{Midpoint}, i)
    nx = docp.dim_NLP_x
    m = docp.dim_NLP_u
    offset = (nx*(1+docp.discretization.stage) + m) * (i-1)
    return @view xu[(offset + 1):(offset + nx)]
end

function get_control_at_time_step(xu, docp::DOCP{Midpoint}, i)
=======
function get_variables_at_time_step(xu, docp::DOCP{Midpoint}, i)
>>>>>>> main
    nx = docp.dim_NLP_x
    m = docp.dim_NLP_u
    N = docp.dim_NLP_steps
<<<<<<< HEAD
    if i < N+1
        offset = (nx*(1+docp.discretization.stage) + m) * (i-1)
    else
        offset = (nx*(1+docp.discretization.stage) + m) * (i-2)
    end
    return @view xu[(offset + nx + 1):(offset + nx + m)]
end
=======
    offset = (nx * (1 + docp.discretization.stage) + m) * i

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
>>>>>>> main

function get_stagevars_at_time_step(xu, docp::DOCP{Midpoint}, i)
    nx = docp.dim_NLP_x
    m = docp.dim_NLP_u
    N = docp.dim_NLP_steps
    offset = (nx*(1+docp.discretization.stage) + m) * (i-1)
    # retrieve vector stage variable (except at final time)
<<<<<<< HEAD
    if i < N+1
        return @view xu[(offset + nx + m + 1):(offset + nx + m + nx) ]
=======
    if i < N
        ki = xu[(offset + nx + m + 1):(offset + nx + m + nx)]
>>>>>>> main
    else
        return nothing
    end
end

<<<<<<< HEAD
function set_state_at_time_step!(xu, x_init, docp::DOCP{Midpoint}, i)
    nx = docp.dim_NLP_x
    n = docp.dim_OCP_x
    m = docp.dim_NLP_u
    offset = (nx*(1+docp.discretization.stage) + m) * (i-1)
    # initialize only actual state variables from OCP (not lagrange state)
=======
# internal NLP version for solution parsing
# could be fused with one above if 
# - using extended dynamics that include lagrange cost
# - scalar case is handled at OCP level
function get_NLP_variables_at_time_step(xu, docp, i, disc::Midpoint)
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
        ki = xu[(offset + nx + m + 1):(offset + nx + m + nx)]
    else
        ki = nothing
    end

    return xi, ui, ki
end

function set_variables_at_time_step!(xu, x_init, u_init, docp, i, disc::Midpoint)
    nx = docp.dim_NLP_x
    n = docp.dim_OCP_x
    m = docp.dim_NLP_u
    N = docp.dim_NLP_steps
    offset = (nx * (1 + docp.discretization.stage) + m) * i

    # NB. only set the actual state variables from the OCP 
    # - skip the possible additional state for lagrange cost
    # - skip internal discretization variables (K_i)
>>>>>>> main
    if !isnothing(x_init)
        xu[(offset + 1):(offset + n)] .= x_init
    end
end

function set_control_at_time_step!(xu, u_init, docp::DOCP{Midpoint}, i)
    nx = docp.dim_NLP_x
    m = docp.dim_NLP_u
    N = docp.dim_NLP_steps
    if i < N+1
        offset = (nx*(1+docp.discretization.stage) + m) * (i-1)
    else
        offset = (nx*(1+docp.discretization.stage) + m) * (i-2)
    end
    if !isnothing(u_init)
        xu[(offset + nx + 1):(offset + nx + m)] .= u_init
    end
end

<<<<<<< HEAD

function setWorkArray(docp::DOCP{Midpoint}, xu, time_grid, v)
    return nothing
=======
# trivial version for now...
# +++multiple dispatch here seems to cause more allocations !
# +++? use abstract type for all Args ?
"""
$(TYPEDSIGNATURES)

Useful values at a time step: time, state, control, dynamics...
"""
struct ArgsAtTimeStep_Midpoint
    time::Any
    state::Any
    control::Any
    lagrange_state::Any
    stage_k::Any
    next_time::Any
    next_state::Any
    next_lagrange_state::Any

    function ArgsAtTimeStep_Midpoint(xu, docp::DOCP{Midpoint}, v, time_grid, i::Int)
        disc = docp.discretization

        # variables
        ti = time_grid[i + 1]
        xi, ui, xli, ki = get_variables_at_time_step(xu, docp, i)

        if i == docp.dim_NLP_steps
            return new(ti, xi, ui, xli, ki, disc)
        else
            tip1 = time_grid[i + 2]
            xip1, uip1, xlip1 = get_variables_at_time_step(xu, docp, i + 1)
            return new(ti, xi, ui, xli, ki, tip1, xip1, xlip1)
        end
    end
end
function initArgs(xu, docp::DOCP{Midpoint}, time_grid)
    v = Float64[]
    docp.has_variable && (v = get_optim_variable(xu, docp))
    args = ArgsAtTimeStep_Midpoint(xu, docp, v, time_grid, 0)
    return args, v
end
function updateArgs(args, xu, docp::DOCP{Midpoint}, v, time_grid, i::Int)
    return ArgsAtTimeStep_Midpoint(xu, docp, v, time_grid, i + 1)
>>>>>>> main
end

"""
$(TYPEDSIGNATURES)

Set the constraints corresponding to the state equation
Convention: 1 <= i <= dim_NLP_steps (+1)
"""
<<<<<<< HEAD
function setConstraintBlock!(docp::DOCP{Midpoint}, c, xu, v, time_grid, i, work)

    disc = docp.discretization

    # offset for previous steps
    offset = (i-1)*(docp.dim_NLP_x * (1+docp.discretization.stage) + docp.dim_path_cons)
=======
function setStateEquation!(docp::DOCP{Midpoint}, c, index::Int, args, v, i)
    
    ocp = docp.ocp
    disc = docp.discretization
>>>>>>> main

    # variables
    ti = time_grid[i]
    xi = get_state_at_time_step(xu, docp, i)
    ui = get_control_at_time_step(xu, docp, i) 

<<<<<<< HEAD
    if i <= docp.dim_NLP_steps
        # more variables
        ki = get_stagevars_at_time_step(xu, docp, i)
        tip1 = time_grid[i+1]
        xip1 = get_state_at_time_step(xu, docp, i+1)
        hi = tip1 - ti

        # midpoint rule
        disc = docp.discretization
        h_sum_bk = hi * disc.butcher_b[1] * ki
        c[offset+1:offset+docp.dim_NLP_x] = xip1 - (xi + h_sum_bk)
        offset += docp.dim_NLP_x

        # stage equation at mid-step
        ts = ti + hi * disc.butcher_c[1]
        #xs = xi + hi * (disc.butcher_a[1][1] * ki)
        xs = 0.5 * (xi + xip1) #compare bench
        if docp.has_inplace
            docp.dynamics_ext((@view c[offset+1:offset+docp.dim_NLP_x]), ts, xs, ui, v)
            @views c[offset+1:offset+docp.dim_NLP_x] = -c[offset+1:offset+docp.dim_NLP_x] + ki
        else
            c[offset+1:offset+docp.dim_NLP_x] = ki - docp.dynamics_ext(ts, xs, ui, v)
        end
        offset += docp.dim_NLP_x
    end

    # path constraints
    # Notes on allocations:.= seems similar
=======
    # midpoint rule
    h_sum_bk = hi * disc.butcher_b[1] * ki[1:docp.dim_NLP_x]
    c[index:(index + docp.dim_OCP_x - 1)] .= xip1 .- (xi .+ h_sum_bk[1:docp.dim_OCP_x])
    # +++ just define extended dynamics !
    if docp.has_lagrange
        c[index + docp.dim_OCP_x] = xlip1 - (xli + h_sum_bk[end])
    end
    index += docp.dim_NLP_x

    # stage equation at mid-step
    t_s = ti + hi * disc.butcher_c[1]
    if docp.dim_OCP_x == 1
        x_s = xi + hi * disc.butcher_a[1][1] * ki[1] #FFS
    else
        x_s = xi .+ hi .* (disc.butcher_a[1][1] .* ki[1:docp.dim_OCP_x])
    end
    c[index:(index + docp.dim_OCP_x - 1)] .= ki[1:(docp.dim_OCP_x)] .- ocp.dynamics(t_s, x_s, ui, v)
    # +++ just define extended dynamics !
    if docp.has_lagrange
        c[index + docp.dim_OCP_x] = ki[end] - ocp.lagrange(t_s, x_s, ui, v)
    end
    index += docp.dim_NLP_x

    return index
end

"""
$(TYPEDSIGNATURES)

Set the path constraints at given time step
"""
function setPathConstraints!(docp::DOCP{Midpoint}, c, index::Int, args, v, i::Int)
    ocp = docp.ocp
    ti = args.time
    xi = args.state
    ui = args.control

    # NB. using .= below *doubles* the allocations oO ??
    # pure control constraints
>>>>>>> main
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

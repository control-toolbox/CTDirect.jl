#= Functions for implicit midpoint discretization scheme
Internal layout for NLP variables: 
[X_0,U_0,K_0, X_1,U_1,K_1 .., X_N-1,U_N-1,K_N-1, XN, V]
with the convention u([t_i,t_i+1[) = U_i and u(tf) = U_N-1
=#

# NB. could be defined as a generic IRK in future IRK.jl for comparison purposes
#butcher_a::Matrix{Float64}
#butcher_b::Vector{Float64}
#butcher_c::Vector{Float64}
#return new(1, 0, hcat(0.5), [1], [0.5], "Implicit Midpoint aka Gauss-Legendre collocation for s=1, 2nd order, symplectic")
struct Midpoint <: Discretization

    stage::Int
    additional_controls::Int  # add control at tf
    info::String

    # constructor
    function Midpoint()
        return new(1, 0, "Implicit Midpoint aka Gauss-Legendre collocation for s=1, 2nd order, symplectic")
    end
end


"""
$(TYPEDSIGNATURES)

Retrieve state and control variables at given time step from the NLP variables.
Convention: 1 <= i <= dim_NLP_steps+1
"""
function get_OCP_state_at_time_step(xu, docp::DOCP{Midpoint,ScalVariable,<:ScalVect,<:ScalVect}, i)
    offset = (i - 1) * (docp.dim_NLP_x * (1 + docp.discretization.stage) + docp.dim_NLP_u)
    return xu[offset+1]
end
function get_OCP_state_at_time_step(xu, docp::DOCP{Midpoint,VectVariable,<:ScalVect,<:ScalVect}, i)
    offset = (i - 1) * (docp.dim_NLP_x * (1 + docp.discretization.stage) + docp.dim_NLP_u)
    return @view xu[(offset+1):(offset+docp.dim_OCP_x)]
end
function get_lagrange_state_at_time_step(xu, docp::DOCP{Midpoint}, i)
    offset = (i - 1) * (docp.dim_NLP_x * (1 + docp.discretization.stage) + docp.dim_NLP_u)
    return xu[offset+docp.dim_NLP_x]
end


function get_OCP_control_at_time_step(xu, docp::DOCP{Midpoint,<:ScalVect,ScalVariable,<:ScalVect}, i)
    if i < docp.dim_NLP_steps + 1
        offset = (i - 1) * (docp.dim_NLP_x * (1 + docp.discretization.stage) + docp.dim_NLP_u) + docp.dim_NLP_x
    else
        offset = (i - 2) * (docp.dim_NLP_x * (1 + docp.discretization.stage) + docp.dim_NLP_u) + docp.dim_NLP_x
    end
    return xu[offset+1]
end
function get_OCP_control_at_time_step(xu, docp::DOCP{Midpoint,<:ScalVect,VectVariable,<:ScalVect}, i)
    if i < docp.dim_NLP_steps + 1
        offset = (i - 1) * (docp.dim_NLP_x * (1 + docp.discretization.stage) + docp.dim_NLP_u) + docp.dim_NLP_x
    else
        offset = (i - 2) * (docp.dim_NLP_x * (1 + docp.discretization.stage) + docp.dim_NLP_u) + docp.dim_NLP_x
    end
    return @view xu[(offset+1):(offset+docp.dim_NLP_u)]
end


function get_stagevars_at_time_step(xu, docp::DOCP{Midpoint}, i)
    offset = (i - 1) * (docp.dim_NLP_x * (1 + docp.discretization.stage) + docp.dim_NLP_u) + docp.dim_NLP_x + docp.dim_NLP_u
    # retrieve vector stage variable (except at final time)
    if i < docp.dim_NLP_steps + 1
        return @view xu[(offset+1):(offset+docp.dim_NLP_x)]
    else
        return @view xu[1:docp.dim_NLP_x] # unused but keep same type !
    end
end


function set_state_at_time_step!(xu, x_init, docp::DOCP{Midpoint}, i)
    offset = (i - 1) * (docp.dim_NLP_x * (1 + docp.discretization.stage) + docp.dim_NLP_u)
    # initialize only actual state variables from OCP (not lagrange state)
    if !isnothing(x_init)
        xu[(offset+1):(offset+docp.dim_OCP_x)] .= x_init
    end
end


function set_control_at_time_step!(xu, u_init, docp::DOCP{Midpoint}, i)
    if i < docp.dim_NLP_steps + 1
        offset = (i - 1) * (docp.dim_NLP_x * (1 + docp.discretization.stage) + docp.dim_NLP_u) + docp.dim_NLP_x
    else
        offset = (i - 2) * (docp.dim_NLP_x * (1 + docp.discretization.stage) + docp.dim_NLP_u) + docp.dim_NLP_x
    end
    if !isnothing(u_init)
        xu[(offset+1):(offset+docp.dim_NLP_u)] .= u_init
    end
end


function setWorkArray(docp::DOCP{Midpoint}, xu, time_grid, v)
    # use work array to store all dynamics + lagrange costs
    work = similar(xu, docp.dim_NLP_x * (docp.dim_NLP_steps))
    if docp.dim_OCP_x > 1
        xs = similar(xu, docp.dim_OCP_x)
    end

    # loop over time steps
    for i = 1:docp.dim_NLP_steps
        offset = (i - 1) * docp.dim_NLP_x
        ti = time_grid[i]
        xi = get_OCP_state_at_time_step(xu, docp, i)
        ui = get_OCP_control_at_time_step(xu, docp, i)
        tip1 = time_grid[i+1]
        xip1 = get_OCP_state_at_time_step(xu, docp, i + 1)
        ts = 0.5 * (ti + tip1)
        # add a new getter ?
        if docp.dim_OCP_x == 1
            xs = 0.5 * (xi + xip1)
        else
            @. xs = 0.5 * (xi + xip1) # not ok for dim 1...
        end
        # OCP dynamics
        docp.ocp.dynamics((@view work[offset+1:offset+docp.dim_OCP_x]), ts, xs, ui, v)
        # lagrange cost
        if docp.is_lagrange
            docp.ocp.lagrange((@view work[offset+docp.dim_NLP_x:offset+docp.dim_NLP_x]), ts, xs, ui, v)
        end
    end
    return work
end


"""
$(TYPEDSIGNATURES)

Set the constraints corresponding to the state equation
Convention: 1 <= i <= dim_NLP_steps (+1)
"""
function setConstraintBlock!(docp::DOCP{Midpoint}, c, xu, v, time_grid, i, work)

    # offset for previous steps
    offset = (i - 1) * (docp.dim_NLP_x * (1 + docp.discretization.stage) + docp.dim_path_cons)

    # 0. variables
    ti = time_grid[i]
    xi = get_OCP_state_at_time_step(xu, docp, i)
    ui = get_OCP_control_at_time_step(xu, docp, i)

    # 1. state equation
    if i <= docp.dim_NLP_steps
        # more variables
        ki = get_stagevars_at_time_step(xu, docp, i)
        tip1 = time_grid[i+1]
        xip1 = get_OCP_state_at_time_step(xu, docp, i + 1)
        hi = tip1 - ti
        offset_dyn_i = (i - 1) * docp.dim_NLP_x

        # midpoint rule
        @views @. c[offset+1:offset+docp.dim_OCP_x] = xip1 - (xi + hi * ki[1:docp.dim_OCP_x])

        if docp.is_lagrange
            c[offset+docp.dim_NLP_x] = get_lagrange_state_at_time_step(xu, docp, i + 1) - (get_lagrange_state_at_time_step(xu, docp, i) + hi * ki[docp.dim_NLP_x])
        end
        offset += docp.dim_NLP_x

        # stage equation at mid-step
        @views @. c[offset+1:offset+docp.dim_OCP_x] = ki[1:docp.dim_OCP_x] - work[offset_dyn_i+1:offset_dyn_i+docp.dim_OCP_x]
        if docp.is_lagrange
            c[offset+docp.dim_NLP_x] = ki[docp.dim_NLP_x] - work[offset_dyn_i+docp.dim_NLP_x]
        end
        offset += docp.dim_NLP_x
    end

    # 2. path constraints
    if docp.dim_u_cons > 0
        docp.control_constraints[2]((@view c[offset+1:offset+docp.dim_u_cons]), ti, ui, v)
    end
    if docp.dim_x_cons > 0
        docp.state_constraints[2]((@view c[offset+docp.dim_u_cons+1:offset+docp.dim_u_cons+docp.dim_x_cons]), ti, xi, v)
    end
    if docp.dim_mixed_cons > 0
        docp.mixed_constraints[2]((@view c[offset+docp.dim_u_cons+docp.dim_x_cons+1:offset+docp.dim_u_cons+docp.dim_x_cons+docp.dim_mixed_cons]), ti, xi, ui, v)
    end

end

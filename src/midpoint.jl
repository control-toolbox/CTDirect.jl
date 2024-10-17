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

    # constructor    
    function Midpoint(dim_NLP_x, dim_NLP_u, dim_NLP_steps) 
        return new(1, 0, hcat(0.5), [1], [0.5], "Implicit Midpoint aka Gauss-Legendre collocation for s=1, 2nd order, symplectic")
    end
end


"""
$(TYPEDSIGNATURES)

Retrieve state and control variables at given time step from the NLP variables.
Convention: 1 <= i <= dim_NLP_steps+1
"""
function get_OCP_state_at_time_step(xu, docp::DOCP{Midpoint, ScalVariable, <: ScalVect, <: ScalVect}, i)
    offset = (i-1) * (docp.dim_NLP_x*(1+docp.discretization.stage) + docp.dim_NLP_u)
    return xu[offset+1]
end
function get_OCP_state_at_time_step(xu, docp::DOCP{Midpoint, VectVariable, <: ScalVect, <: ScalVect}, i)
    offset = (i-1) * (docp.dim_NLP_x*(1+docp.discretization.stage) + docp.dim_NLP_u)
    return @view xu[(offset + 1):(offset + docp.dim_OCP_x)]
end
function get_lagrange_state_at_time_step(xu, docp::DOCP{Midpoint}, i)
    offset = (i-1) * (docp.dim_NLP_x*(1+docp.discretization.stage) + docp.dim_NLP_u)
    return xu[offset + docp.dim_NLP_x]
end


function get_OCP_control_at_time_step(xu, docp::DOCP{Midpoint, <: ScalVect, ScalVariable, <: ScalVect}, i)
    if i < docp.dim_NLP_steps+1
        offset = (i-1) * (docp.dim_NLP_x*(1+docp.discretization.stage) + docp.dim_NLP_u) + docp.dim_NLP_x
    else
        offset = (i-2) * (docp.dim_NLP_x*(1+docp.discretization.stage) + docp.dim_NLP_u) + docp.dim_NLP_x
    end
    return xu[offset+1]
end
function get_OCP_control_at_time_step(xu, docp::DOCP{Midpoint, <: ScalVect, VectVariable, <: ScalVect}, i)
    if i < docp.dim_NLP_steps+1
        offset = (i-1) * (docp.dim_NLP_x*(1+docp.discretization.stage) + docp.dim_NLP_u) + docp.dim_NLP_x
    else
        offset = (i-2) * (docp.dim_NLP_x*(1+docp.discretization.stage) + docp.dim_NLP_u) + docp.dim_NLP_x
    end
    return @view xu[(offset + 1):(offset + docp.dim_NLP_u)]
end


function get_stagevars_at_time_step(xu, docp::DOCP{Midpoint}, i)
    nx = docp.dim_NLP_x
    m = docp.dim_NLP_u
    N = docp.dim_NLP_steps
    offset = (nx*(1+docp.discretization.stage) + m) * (i-1)
    # retrieve vector stage variable (except at final time)
    if i < N+1
        return @view xu[(offset + nx + m + 1):(offset + nx + m + nx) ]
    else
        return nothing
    end
end

function set_state_at_time_step!(xu, x_init, docp::DOCP{Midpoint}, i)
    nx = docp.dim_NLP_x
    n = docp.dim_OCP_x
    m = docp.dim_NLP_u
    offset = (nx*(1+docp.discretization.stage) + m) * (i-1)
    # initialize only actual state variables from OCP (not lagrange state)
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


function setWorkArray(docp::DOCP{Midpoint}, xu, time_grid, v)
    return nothing
end

"""
$(TYPEDSIGNATURES)

Set the constraints corresponding to the state equation
Convention: 1 <= i <= dim_NLP_steps (+1)
"""
function setConstraintBlock!(docp::DOCP{Midpoint}, c, xu, v, time_grid, i, work)

    disc = docp.discretization

    # offset for previous steps
    offset = (i-1)*(docp.dim_NLP_x * (1+docp.discretization.stage) + docp.dim_path_cons)

    # variables
    ti = time_grid[i]
    xi = get_OCP_state_at_time_step(xu, docp, i)
    ui = get_OCP_control_at_time_step(xu, docp, i) 

    # 1. state equation
    if i <= docp.dim_NLP_steps
        # more variables
        ki = get_stagevars_at_time_step(xu, docp, i)
        tip1 = time_grid[i+1]
        xip1 = get_OCP_state_at_time_step(xu, docp, i+1)
        hi = tip1 - ti

        # midpoint rule
        h_sum_bk = hi * disc.butcher_b[1] * ki
        if docp.dim_OCP_x == 1
            c[offset+1] = xip1 - (xi + h_sum_bk[1])
        else
            c[offset+1:offset+docp.dim_OCP_x] = xip1 - (xi + h_sum_bk[1:docp.dim_OCP_x])
        end
        if docp.is_lagrange
            c[offset+docp.dim_NLP_x] = get_lagrange_state_at_time_step(xu, docp, i+1) - (get_lagrange_state_at_time_step(xu, docp, i) + h_sum_bk[docp.dim_NLP_x])
        end
        offset += docp.dim_NLP_x

        # stage equation at mid-step
        ts = ti + hi * disc.butcher_c[1]
        #xs = xi + hi * (disc.butcher_a[1][1] * ki)
        xs = 0.5 * (xi + xip1) #compare bench
        if docp.is_inplace
            docp.ocp.dynamics((@view c[offset+1:offset+docp.dim_OCP_x]), ts, xs, ui, v)
            @views c[offset+1:offset+docp.dim_OCP_x] = -c[offset+1:offset+docp.dim_OCP_x] + ki[1:docp.dim_OCP_x]
        else
            c[offset+1:offset+docp.dim_OCP_x] .= ki[1:docp.dim_OCP_x] .- docp.ocp.dynamics(ts, xs, ui, v)
        end
        if docp.is_lagrange
            if docp.is_inplace
                docp.ocp.lagrange((@view c[offset+docp.dim_NLP_x:offset+docp.dim_NLP_x]), ts, xs, ui, v)
                @views c[offset+docp.dim_NLP_x] = -c[offset+docp.dim_NLP_x] + ki[docp.dim_NLP_x]
            else
                c[offset+docp.dim_NLP_x] = ki[docp.dim_NLP_x] - docp.ocp.lagrange(ts, xs, ui, v)
            end
        end
        offset += docp.dim_NLP_x
    end

    # 2. path constraints
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

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

    get_state_at_time_step::Function
    get_control_at_time_step::Function
    get_stagevars_at_time_step::Function

    # constructor    
    function Midpoint(dim_NLP_x, dim_NLP_u, dim_NLP_steps) 
        
        stage = 1

        # getters for state and control variables
        get_state_at_time_step = function (xu, i)
            offset = (dim_NLP_x*(1+stage) + dim_NLP_u) * (i-1)
            return @view xu[(offset + 1):(offset + dim_NLP_x)]
        end

        get_control_at_time_step = function (xu, i)
            if i < dim_NLP_steps+1
                offset = (dim_NLP_x*(1+stage) + dim_NLP_u) * (i-1) + dim_NLP_x
            else
                offset = (dim_NLP_x*(1+stage) + dim_NLP_u) * (i-2) + dim_NLP_x
            end
            return @view xu[(offset + 1):(offset + dim_NLP_u)]
        end

        get_stagevars_at_time_step = function (xu, i)
            # retrieve vector stage variable (except at final time)
            if i < dim_NLP_steps+1
                offset = (dim_NLP_x *(1+stage) + dim_NLP_u) * (i-1) + dim_NLP_x  + dim_NLP_u
                return @view xu[(offset + 1):(offset + dim_NLP_x ) ]
            else
                return nothing
            end
        end

        return new(stage, 0, hcat(0.5), [1], [0.5], "Implicit Midpoint aka Gauss-Legendre collocation for s=1, 2nd order, symplectic", get_state_at_time_step, get_control_at_time_step, get_stagevars_at_time_step)
    end
end

#=
"""
$(TYPEDSIGNATURES)

Retrieve state and control variables at given time step from the NLP variables.
Convention: 1 <= i <= dim_NLP_steps+1
"""
function get_state_at_time_step(xu, docp::DOCP{Midpoint}, i)
    nx = docp.dim_NLP_x
    m = docp.dim_NLP_u
    offset = (nx*(1+docp.discretization.stage) + m) * (i-1)
    return @view xu[(offset + 1):(offset + nx)]
end

function get_control_at_time_step(xu, docp::DOCP{Midpoint}, i)
    nx = docp.dim_NLP_x
    m = docp.dim_NLP_u
    N = docp.dim_NLP_steps
    if i < N+1
        offset = (nx*(1+docp.discretization.stage) + m) * (i-1)
    else
        offset = (nx*(1+docp.discretization.stage) + m) * (i-2)
    end
    return @view xu[(offset + nx + 1):(offset + nx + m)]
end

function get_ki_at_time_step(xu, docp::DOCP{Midpoint}, i)
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
end=#

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
    xi = disc.get_state_at_time_step(xu, i)
    ui = disc.get_control_at_time_step(xu, i)
  
    if i <= docp.dim_NLP_steps
        # more variables
        ki = disc.get_stagevars_at_time_step(xu, i)
        tip1 = time_grid[i+1]
        xip1 = disc.get_state_at_time_step(xu, i+1)
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
    #setPathConstraints!(docp, c, ti, xi, ui, v, offset)
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

#= Functions for trapeze discretization scheme
Internal layout for NLP variables: 
[X_0,U_0, X_1,U_1, .., X_N,U_N, V]
=#

# NB. could be defined as a generic IRK
struct Trapeze <: Discretization

    stage::Int
    additional_controls::Int  # add control at tf
    info::String

    Trapeze() = new(0, 1, "Implicit Trapeze aka Crank-Nicolson, 2nd order, A-stable")
end


"""
$(TYPEDSIGNATURES)

Retrieve state and control variables at given time step from the NLP variables.
Convention: 1 <= i <= dim_NLP_steps+1
"""
function get_variables_at_time_step(xu, docp::DOCP{Trapeze}, i)

    nx = docp.dim_NLP_x
    n = docp.dim_OCP_x
    m = docp.dim_NLP_u
    offset = (nx + m) * (i-1)

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
   
    f = similar(xu, docp.dim_NLP_x)

    ocp = docp.ocp
    t0 = time_grid[1]
    x0, u0 = get_variables_at_time_step(xu, docp, 1)
    
    if docp.has_inplace
        #+++ inplace view
    else
        f[1:docp.dim_OCP_x] .= ocp.dynamics(t0, x0, u0, v)
    end
    
    if docp.has_lagrange
        if docp.has_inplace
            #+++ inplace
        else
            f[docp.dim_NLP_x] = ocp.lagrange(t0, x0, u0, v)
        end
    end

    return f
end


"""
$(TYPEDSIGNATURES)

Set the constraints corresponding to the state equation
Convention: 1 <= i <= dim_NLP_steps
"""
function setConstraintBlock!(docp::DOCP{Trapeze}, c, xu, v, time_grid, i, work)

    # offset for previous steps
    offset = (i-1)*(docp.dim_NLP_x + docp.dim_path_cons)

    # variables
    ocp = docp.ocp
    ti = time_grid[i]
    xi, ui, xli = get_variables_at_time_step(xu, docp, i)
    fi = work[1:docp.dim_OCP_x]

    tip1 = time_grid[i+1]
    xip1, uip1, xlip1 = get_variables_at_time_step(xu, docp, i+1)
    if docp.has_inplace
        #+++ inplace view ?
    else
        fip1 = ocp.dynamics(tip1, xip1, uip1, v)
    end
    hi = tip1 - ti

    # trapeze rule with 'smart' update for dynamics
    c[offset+1:offset+docp.dim_OCP_x] .= xip1 .- (xi .+ 0.5 * hi * (fi .+ fip1))
    work[1:docp.dim_OCP_x] .= fip1

    if docp.has_lagrange
        li = work[docp.dim_OCP_x+1]
        if docp.has_inplace
            #+++ inplace
        else
            lip1 = ocp.lagrange(tip1, xip1, uip1, v)
        end
        c[offset+docp.dim_OCP_x+1] = xlip1 - (xli + 0.5 * hi * (li + lip1))
        work[docp.dim_NLP_x] = lip1
    end
    offset += docp.dim_NLP_x    

    # path constraints
    setPathConstraints!(docp, c, ti, xi, ui, v, offset)

end

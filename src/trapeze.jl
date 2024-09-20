#= Functions for trapeze discretization scheme
Internal layout for NLP variables: 
[X_0,U_0, X_1,U_1, .., X_N,U_N, V]
=#

# NB. could be defined as a generic IRK
struct Trapeze <: Discretization
<<<<<<< HEAD

    stage::Int
    additional_controls::Int  # add control at tf
    info::String

    # constructor
    function Trapeze(dim_NLP_x, dim_NLP_u)
        return new(0, 1, "Implicit Trapeze aka Crank-Nicolson, 2nd order, A-stable")
    end
end
=======
>>>>>>> main

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
<<<<<<< HEAD
function get_state_at_time_step(xu, docp::DOCP{Trapeze}, i)
    nx = docp.dim_NLP_x
    m = docp.dim_NLP_u
    offset = (nx + m) * (i-1)
    return @view xu[(offset + 1):(offset + nx)]
end

function get_control_at_time_step(xu, docp::DOCP{Trapeze}, i)
    nx = docp.dim_NLP_x
    m = docp.dim_NLP_u
    offset = (nx + m) * (i-1)
    return @view xu[(offset + nx + 1):(offset + nx + m)]
end

function set_state_at_time_step!(xu, x_init, docp::DOCP{Trapeze}, i)
    nx = docp.dim_NLP_x
    n = docp.dim_OCP_x
    m = docp.dim_NLP_u
    offset = (nx + m) * (i-1)
    # initialize only actual state variables from OCP (not lagrange state)
=======
function get_variables_at_time_step(xu, docp::DOCP{Trapeze}, i)
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
function get_NLP_variables_at_time_step(xu, docp, i, disc::Trapeze)
    nx = docp.dim_NLP_x
    m = docp.dim_NLP_u
    offset = (nx + m) * i

    # state
    xi = xu[(offset + 1):(offset + nx)]
    # control
    ui = xu[(offset + nx + 1):(offset + nx + m)]

    return xi, ui
end

function set_variables_at_time_step!(xu, x_init, u_init, docp, i, disc::Trapeze)
    nx = docp.dim_NLP_x
    n = docp.dim_OCP_x
    m = docp.dim_NLP_u
    N = docp.dim_NLP_steps
    offset = (nx + m) * i

    # NB. only set the actual state variables from the OCP 
    # - skip the possible additional state for lagrange cost
>>>>>>> main
    if !isnothing(x_init)
        xu[(offset + 1):(offset + n)] .= x_init
    end
end

function set_control_at_time_step!(xu, u_init, docp::DOCP{Trapeze}, i)
    nx = docp.dim_NLP_x
    m = docp.dim_NLP_u
    offset = (nx + m) * (i-1)
    if !isnothing(u_init)
        xu[(offset + nx + 1):(offset + nx + m)] .= u_init
    end
end

<<<<<<< HEAD

function setWorkArray(docp::DOCP{Trapeze}, xu, time_grid, v)
   
    disc = docp.discretization
=======
# ? use abstract type for args ?
"""
$(TYPEDSIGNATURES)
>>>>>>> main

    work = similar(xu, docp.dim_NLP_x)
    t0 = time_grid[1]
    x0 = get_state_at_time_step(xu, docp, 1)
    u0 = get_control_at_time_step(xu, docp, 1)

<<<<<<< HEAD
    if docp.has_inplace
        docp.dynamics_ext(work, t0, x0, u0, v)
=======
    function ArgsAtTimeStep_Trapeze(xu, docp::DOCP{Trapeze}, v, time_grid, i::Int)

        # variables
        ti = time_grid[i + 1]
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
function initArgs(xu, docp::DOCP{Trapeze}, time_grid)
    # optimization variables
    v = Float64[]
    docp.has_variable && (v = get_optim_variable(xu, docp))
    args_i = ArgsAtTimeStep_Trapeze(xu, docp, v, time_grid, 0)
    args_ip1 = ArgsAtTimeStep_Trapeze(xu, docp, v, time_grid, 1)
    return (args_i, args_ip1), v
end
function updateArgs(args, xu, docp::DOCP{Trapeze}, v, time_grid, i)
    args_i, args_ip1 = args
    if i < docp.dim_NLP_steps - 1
        # are we allocating more than one args here ?
        return (args_ip1, ArgsAtTimeStep_Trapeze(xu, docp, v, time_grid, i + 2))
>>>>>>> main
    else
        # NB. work = will create a new variable ;-) (work .= is fine)
        work[:] = docp.dynamics_ext(t0, x0, u0, v)
    end
    return work
end

"""
$(TYPEDSIGNATURES)

Set the constraints corresponding to the state equation
Convention: 1 <= i <= dim_NLP_steps (+1)
"""
function setConstraintBlock!(docp::DOCP{Trapeze}, c, xu, v, time_grid, i, work)

    disc = docp.discretization

    # offset for previous steps
    offset = (i-1)*(docp.dim_NLP_x + docp.dim_path_cons)

    # 0. variables
    ti = time_grid[i]
    xi = get_state_at_time_step(xu, docp, i)
    ui = get_control_at_time_step(xu, docp, i)

    #1. state equation
    if i <= docp.dim_NLP_steps
        # more variables
        fi = copy(work) # create new copy, not just a reference
        tip1 = time_grid[i+1]
        xip1 = get_state_at_time_step(xu, docp, i+1)
        uip1 = get_control_at_time_step(xu, docp, i+1)
        if docp.has_inplace
            docp.dynamics_ext(work, tip1, xip1, uip1, v)
        else
            # copy, do not create a new variable !
            work[:] = docp.dynamics_ext(tip1, xip1, uip1, v)
        end

        # trapeze rule with 'smart' update for dynamics (similar with @.)
        c[offset+1:offset+docp.dim_NLP_x] = xip1 - (xi + 0.5 * (tip1 - ti) * (fi + work))
        offset += docp.dim_NLP_x
    end

<<<<<<< HEAD
    # 2. path constraints
    # Notes on allocations:.= seems similar
=======
    return index
end

"""
$(TYPEDSIGNATURES)

Set the path constraints at given time step
"""
function setPathConstraints!(docp::DOCP{Trapeze}, c, index::Int, args, v, i::Int)

    # note: i is unused but passed for call compatibility
    ocp = docp.ocp
    args_i, args_ip1 = args
    ti = args_i.time
    xi = args_i.state
    ui = args_i.control

    # NB. using .= below *doubles* the allocations oO
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

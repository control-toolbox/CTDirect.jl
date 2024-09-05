#= Functions for trapeze discretization scheme
Internal layout for NLP variables: 
[X_0,U_0, X_1,U_1, .., X_N,U_N, V]
=#

# +++ todo change arguments order: docp first, then xu


struct Trapeze <: Discretization
    stage::Int
    additional_controls::Int
    # add control at tf
    Trapeze() = new(0, 1)
end


"""
$(TYPEDSIGNATURES)

Retrieve state and control variables at given time step from the NLP variables. 
"""
function get_variables_at_t_i(xu, docp::DOCP{Trapeze}, i)

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


#=
# ? use abstract type for args ?
"""
$(TYPEDSIGNATURES)

Useful values at a time step: time, state, control, dynamics...
"""
struct ArgsAtTimeStep_Trapeze
    time::Any
    state::Any
    control::Any
    dynamics::Any
    lagrange_state::Any
    lagrange_cost::Any

    function ArgsAtTimeStep_Trapeze(xu, docp::DOCP{Trapeze}, v, time_grid, i::Int)

        # variables
        ti = time_grid[i+1]
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
    v = Float64[][]
    docp.has_variable && (v = get_optim_variable(xu, docp))
    args_i = ArgsAtTimeStep_Trapeze(xu, docp, v, time_grid, 0)
    args_ip1 = ArgsAtTimeStep_Trapeze(xu, docp, v, time_grid, 1)
    return (args_i, args_ip1), v
end
function updateArgs(args, xu, docp::DOCP{Trapeze}, v, time_grid, i)
    args_i, args_ip1 = args
    if i < docp.dim_NLP_steps - 1
        # are we allocating more than one args here ?
        return (args_ip1, ArgsAtTimeStep_Trapeze(xu, docp, v, time_grid, i+2))
    else
        return (args_ip1, args_ip1)
    end
end
=#

struct Trapeze_Args
    variable
    time
    state
    control
    dynamics
    next_time
    next_state
    next_dynamics
    lagrange_state
    lagrange_cost
    next_lagrange_state
    next_lagrange_cost
end

function initArgs(docp::DOCP{Trapeze}, xu)
    
    # dynamics / lagange costs are present to avoid recomputation
    # in trapeze method
    args = Vector{Trapeze_Args}(undef, docp.dim_NLP_steps + 1)
    dummy = similar(xu,0)

    # get time grid
    time_grid = get_time_grid(xu, docp)

    #+++NB. we could compute all the path constraints here, so that passing v and u in args could be removed ? but then we have to pass the path constraint block to setstateequation...
    
    # get optim variable
    if docp.has_variable
        v = get_optim_variable(xu, docp)
    else
        v = dummy
    end

    # evaluate all dynamics and lagrange costs
    dynamics_vec = Matrix{eltype(xu[1])}(undef, docp.dim_OCP_x, docp.dim_NLP_steps+1)
    docp.has_lagrange && (lagrange_vec = Vector{eltype(xu[1])}(undef, docp.dim_NLP_steps+1))
    for i = 1:docp.dim_NLP_steps+1
        # NB we call the getters twice for t, x and u ...
        t_i = time_grid[i]
        x_i, u_i, = get_variables_at_t_i(xu, docp, i-1)
        dynamics_vec[:,i] = docp.ocp.dynamics(t_i, x_i, u_i, v)
        docp.has_lagrange && (lagrange_vec[i] = docp.ocp.lagrange(t_i, x_i, u_i, v))
    end

    # loop over time steps
    for i = 1:docp.dim_NLP_steps
        t_i = time_grid[i]
        t_ip1 = time_grid[i+1]
        x_i, u_i, xl_i = get_variables_at_t_i(xu, docp, i-1)
        x_ip1, u_ip1, xl_ip1 = get_variables_at_t_i(xu, docp, i)

        # NB. view are not better...
        f_i = dynamics_vec[:,i]
        f_ip1 = dynamics_vec[:,i+1]
        if docp.has_lagrange
            l_i = lagrange_vec[i]
            l_ip1 = lagrange_vec[i+1]
        else
            l_i = dummy
            l_ip1 = dummy
        end

        # set args
        args[i] = Trapeze_Args(v, t_i, x_i, u_i, f_i, t_ip1, x_ip1, f_ip1, xl_i, l_i, xl_ip1, l_ip1)

    end

    # final time: for path constraints only
    # useful fields are: v, ti, xi, ui
    t_f = time_grid[docp.dim_NLP_steps+1]
    x_f, u_f = get_variables_at_t_i(xu, docp, docp.dim_NLP_steps)
    args[docp.dim_NLP_steps+1] = Trapeze_Args(v, t_f, x_f, u_f, dummy, dummy, dummy, dummy, dummy, dummy, dummy, dummy)

    return args
end

"""
$(TYPEDSIGNATURES)

Set the constraints corresponding to the state equation
"""
function setStateEquation!(docp::DOCP{Trapeze}, c_block, args::Trapeze_Args)

    step_size = args.next_time - args.time

    # trapeze rule (NB. @. allocates more ...)
    c_block[1:docp.dim_OCP_x] .=
        args.next_state .- (args.state .+ 0.5 * step_size * (args.dynamics .+ args.next_dynamics))
    
        # +++ just define extended dynamics !
    if docp.has_lagrange
        c_block[docp.dim_OCP_x+1] =
            args.next_lagrange_state -
            (args.lagrange_state + 0.5 * step_size * (args.lagrange_cost + args.next_lagrange_cost))
    end
    
end


"""
$(TYPEDSIGNATURES)

Set the path constraints at given time step
"""
function setPathConstraints!(docp::DOCP{Trapeze}, c_block, args::Trapeze_Args)    

    # Arguments list for one time step: 
    v = args.variable
    t_i = args.time
    x_i = args.state
    u_i = args.control

    ocp = docp.ocp

    # NB. using .= below *doubles* the allocations oO
    # pure control constraints
    # WAIT FOR INPLACE VERSION IN OCP !
    if docp.dim_u_cons > 0
        c_block[1:docp.dim_u_cons] = docp.control_constraints[2](t_i, u_i, v)
    end

    # pure state constraints
    if docp.dim_x_cons > 0
        c_block[docp.dim_u_cons+1:docp.dim_u_cons+docp.dim_x_cons] = docp.state_constraints[2](t_i, x_i, v)
    end

    # mixed state / control constraints
    if docp.dim_mixed_cons > 0
        c_block[docp.dim_u_cons+docp.dim_x_cons+1:docp.dim_u_cons+docp.dim_x_cons+docp.dim_mixed_cons] = docp.mixed_constraints[2](t_i, x_i, u_i, v)
    end

end

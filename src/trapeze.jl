#= Functions for trapeze discretization scheme
Internal layout for NLP variables: 
[X_0,U_0, X_1,U_1, .., X_N,U_N, V]
=#


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
        return (args_ip1, ArgsAtTimeStep_Trapeze(xu, docp, v, time_grid, i+2))
    else
        return (args_ip1, args_ip1)
    end
end
=#

function initArgs(docp::DOCP{Trapeze}, xu)
    
    # Arguments list for one time step: 
    # t_i, x_i, u_i, f_i, t_ip1, x_ip1, f_ip1 
    # [, xl_i, l_i, xl_ip1, l_ip1]
    # dynamics / lagange costs are present to avoid recomputation
    # in trapeze method
    args = Vector{}(undef, docp.dim_NLP_steps + 1)

    # get time grid
    time_grid = get_time_grid(xu, docp)

    #+++NB. we could compute all the path constraints here, so that passing v and u in args could be removed ? but then we have to pass the path constraint block to setstateequation...
    
    +++ computing the whole dynamics may reduce the update part a bit
    +++ later use sub function, maybe pass args too to avoid multiple calls to getters for t,x,u ?    
    dynamics_vec vect(vect) N+1 x dim_ocp_x, eltype(xu[0]) ?

    # get optim variable
    if docp.has_variable
        v = get_optim_variables(xu, docp)
    else
        v = Float64[]
    end

    # fill args vector


    # loop over time steps
    for i = 1:docp.dim_NLP_steps
        t_i = time_grid[i]
        t_ip1 = time_grid[i+1]
        if docp.has_lagrange
            x_i, u_i, xl_i = get_variables_at_t_i(xu, docp, i-1)
            l_i = docp.ocp.lagrange(t_i, x_i, u_i, v)
            x_ip1, u_ip1, xl_ip1 = get_variables_at_t_i(xu, docp, i)
            l_ip1 = docp.ocp.lagrange(t_ip1, x_ip1, u_ip1, v)
        else
            x_i, u_i = get_variables_at_t_i(xu, docp, i-1)
            x_ip1, u_ip1 = get_variables_at_t_i(xu, docp, i)
        end
        #+++ compute whole vector before loop
        #f_i = docp.ocp.dynamics(t_i, x_i, u_i, v)
        #f_ip1 = docp.ocp.dynamics(t_ip1, x_ip1, u_ip1, v)

        # set args
        h_i = t_ip1 - t_i
        f_i = dynamics_vec[i]
        f_ip1 = dynamics_vec[i+1]
        if docp.has_lagrange
            args[i] = (v, t_i, x_i, u_i, f_i, t_ip1, x_ip1, f_ip1, xl_i, l_i, xl_ip1, l_ip1)
        else
            args[i] = (v, t_i, x_i, u_i, f_i, t_ip1, x_ip1, f_ip1)
        end

        #+++ remove this 'smart' update
        #t_i, x_i, u_i, f_i = t_ip1, x_ip1, u_ip1, f_ip1
        #docp.has_lagrange && (xl_i, l_i = xl_ip1, l_ip1)
    end

    # final time: for path constraints only
    # useful fields are: v, ti, xi, ui 
    args[docp.dim_NLP_steps+1] = (v, t_i, x_i, u_i, f_i, t_i, x_i, f_i) 

    return args
end

"""
$(TYPEDSIGNATURES)

Set the constraints corresponding to the state equation
"""
#function setStateEquation!(docp::DOCP{Trapeze}, c, index::Int, args, v, i)
function setStateEquation!(docp::DOCP{Trapeze}, c, index::Int, args)

    # Arguments list for one time step:
    # v, t_i, x_i, u_i, f_i, t_ip1, x_ip1, f_ip1 
    # [, xl_i, l_i, xl_ip1, l_ip1]
    if docp.has_lagrange
        _, time, state, _, dynamics, next_time, next_state, next_dynamics, lagrange_state, lagrange_cost, next_lagrange_state, next_lagrange_cost = args
    else
        time, state, control, dynamics, next_time, next_state, next_dynamics = args
    end

    # trapeze rule (NB. @. allocates more ...)
    c[index:(index + docp.dim_OCP_x - 1)] .=
        next_state .- (state .+ 0.5 * (next_time - time) * (dynamics .+ next_dynamics))
    
        # +++ just define extended dynamics !
    if docp.has_lagrange
        c[index + docp.dim_OCP_x] =
            next_lagrange_state -
            (lagrange_state + 0.5 * (next_time - time) * (lagrange_cost + next_lagrange_cost))
    end
    
    index += docp.dim_NLP_x
    return index
end


"""
$(TYPEDSIGNATURES)

Set the path constraints at given time step
"""
#function setPathConstraints!(docp::DOCP{Trapeze}, c, index::Int, args, v, i::Int)
function setPathConstraints!(docp::DOCP{Trapeze}, c, index::Int, args)    

    # Arguments list for one time step: 
    # v, t_i, x_i, u_i, f_i, t_ip1, x_ip1, f_ip1 
    # [, xl_i, l_i, xl_ip1, l_ip1]
    v, t_i, x_i, u_i, _, _, _, _  = args

    ocp = docp.ocp

    # NB. using .= below *doubles* the allocations oO
    # pure control constraints
    # WAIT FOR INPLACE VERSION IN OCP !
    if docp.dim_u_cons > 0
        c[index:(index + docp.dim_u_cons - 1)] = docp.control_constraints[2](t_i, u_i, v)
        index += docp.dim_u_cons
    end

    # pure state constraints
    if docp.dim_x_cons > 0
        c[index:(index + docp.dim_x_cons - 1)] = docp.state_constraints[2](t_i, x_i, v)
        index += docp.dim_x_cons
    end

    # mixed state / control constraints
    if docp.dim_mixed_cons > 0
        c[index:(index + docp.dim_mixed_cons - 1)] = docp.mixed_constraints[2](t_i, x_i, u_i, v)
        index += docp.dim_mixed_cons
    end

    return index
end

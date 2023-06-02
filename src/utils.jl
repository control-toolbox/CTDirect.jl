function get_variable(xu, ctd)
    if ctd.has_variable
        if ctd.variable_dimension == 1
            return xu[end]
        else
            return xu[end-ctd.variable_dimension+1:end]
        end
    else
        return Real[]
    end
end

# return original ocp state
function get_state_at_time_step(xu, ctd, i)
    """
        return
        x(t_i)
    """
    nx = ctd.dim_NLP_state
    n = ctd.state_dimension
    N = ctd.dim_NLP_steps
    @assert i <= N "trying to get x(t_i) for i > N"
    if n == 1
        return xu[i*nx + 1]
    else
        return xu[i*nx + 1 : i*nx + n]
    end
end

function get_lagrange_cost_at_time_step(xu, ctd, i)
    nx = ctd.dim_NLP_state
    N = ctd.dim_NLP_steps
    @assert i <= N "trying to get lagrange cost at t_i for i > N"
    return xu[(i+1)*nx]
end

function vget_state_at_time_step(xu, ctd, i)
    nx = ctd.dim_NLP_state
    N = ctd.dim_NLP_steps
    @assert i <= N "trying to get x(t_i) for i > N"
    return xu[i*nx + 1 : (i+1)*nx]
end

function get_control_at_time_step(xu, ctd, i)
    """
        return
        u(t_i)
    """
    nx = ctd.dim_NLP_state
    m = ctd.control_dimension
    N = ctd.dim_NLP_steps
    @assert i <= N "trying to get u(t_i) for i > N"
    if m == 1
        return xu[(N+1)*nx + i*m + 1]
    else
        return xu[(N+1)*nx + i*m + 1 : (N+1)*nx + (i+1)*m]
    end
end

function vget_control_at_time_step(xu, ctd, i)
    nx = ctd.dim_NLP_state
    m = ctd.control_dimension
    N = ctd.dim_NLP_steps
    @assert i <= N "trying to get u(t_i) for i > N"
    return xu[(N+1)*nx + i*m + 1 : (N+1)*nx + (i+1)*m]
end

function get_initial_time(xu, ctd)
    if ctd.has_free_initial_time
        v = get_variable(xu, ctd)
        return v[ctd.initial_time]
    else
        return ctd.initial_time
    end
end

function get_final_time(xu, ctd)
    if ctd.has_free_final_time
        v = get_variable(xu, ctd)
        return v[ctd.final_time]
    else
        return ctd.final_time
    end
end

#=
function get_time_step(xu, ctd, i)
    N = ctd.dim_NLP_steps
    @assert i <= N "trying to get time step for i > N"
    t0 = get_initial_time(xu, ctd)
    tf = get_final_time(xu, ctd)
    return t0 + i * (tf - t0) / N 
end
=#

## Initialization for the NLP problem

function set_state_at_time_step!(xu, x_init, ctd, i)
    nx = ctd.dim_NLP_state
    n = ctd.state_dimension
    N = ctd.dim_NLP_steps
    @assert i <= N "trying to set init for x(t_i) with i > N"
    # NB. only set first n components of state variable (nx = n+1 for lagrange cost)
    if n == 1
        xu[i*nx + 1] = x_init[]
    else
        xu[i*nx + 1 : i*nx + n] = x_init
    end
end
    
function set_control_at_time_step!(xu, u_init, ctd, i)
    nx = ctd.dim_NLP_state
    m = ctd.control_dimension
    N = ctd.dim_NLP_steps
    @assert i <= N "trying to set init for u(t_i) with i > N"
    offset = (N+1)*nx
    if m == 1
        xu[offset + i*m + 1] = u_init[]
    else        
        xu[offset + i*m + 1 : offset + i*m + m] = u_init
    end
end

function set_variable!(xu, v_init, ctd)
    if ctd.variable_dimension == 1
        xu[end] = v_init[]
    else
        xu[end-ctd.variable_dimension+1 : end] = v_init
    end
end

function initial_guess(ctd)

    # default initialization
    # note: internal variables (lagrange cost, k_i for RK schemes) will keep these default values 
    xu0 = 0.1 * ones(ctd.dim_NLP_variables)

    init = ctd.NLP_init
    if init.info != :undefined
        N = ctd.dim_NLP_steps
        t0 = get_initial_time(xu0, ctd)
        tf = get_final_time(xu0, ctd)
        h = (tf - t0) / N 

        # set state / control variables
        for i in 0:N
            ti = t0 + i * h
            set_state_at_time_step!(xu0, init.state_init(ti), ctd, i)
            set_control_at_time_step!(xu0, init.control_init(ti), ctd, i)
        end

        # set variables
        if (ctd.variable_dimension > 0)
            set_variable!(xu0, init.variable_init, ctd)
        end
    end

    return xu0
end
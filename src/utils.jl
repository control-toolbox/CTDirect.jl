function get_variable(xu, ctd)
    # put variable at the end of NLP unknown ?
    if ctd.has_variable
        v = xu[end - ctd.variable_dimension + 1 : end]
    else
        v = Real[]
    end
    return v
end

# return augmented state including potential component for lagrange objective
function get_augmented_state_at_time_step(xu, ctd, i)
    """
        return
        x(t_i)
    """
    nx = ctd.dim_NLP_state
    N = ctd.dim_NLP_steps
    @assert i <= N "trying to get x(t_i) for i > N"
    if nx == 1
        return xu[i*nx + 1]
    else
        return xu[i*nx + 1 : (i+1)*nx]
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
        v = get_variable(xu, ctd.variable_dimension)
        return v[ctd.initial_time]
    else
        return ctd.initial_time
    end
end

function get_final_time(xu, ctd)
    if ctd.has_free_final_time
        v = get_variable(xu, ctd.variable_dimension)
        return v[ctd.final_time]
    else
        return ctd.final_time
    end
end


## Initialization for the NLP problem

function set_state_at_time_step!(xu, x_init, ctd, i)
    nx = ctd.dim_NLP_state
    N = ctd.dim_NLP_steps
    @assert i <= N "trying to set init for x(t_i) with i > N"
    xu[1+i*nx:(i+1)*nx] = x_init[1:nx]
end
    
function set_control_at_time_step!(xu, u_init, ctd, i)
    nx = ctd.dim_NLP_state
    m = ctd.control_dimension
    N = ctd.dim_NLP_steps
    @assert i <= N "trying to set init for u(t_i) with i > N"
    xu[1+(N+1)*nx+i*m:m+(N+1)*nx+i*m] = u_init[1:m]
end

function set_variable!(xu, v_init, ctd)
    xu[end-ctd.variable_dimension+1:end] = v_init[1:ctd.variable_dimension]
end

function initial_guess(ctd)

    N = ctd.dim_NLP_steps
    init = ctd.NLP_init

    if init === nothing
        # default initialization
        xu0 = 0.1*ones(ctd.dim_NLP_variables)
    else
        # split state / control init values
        if length(init) != (ctd.state_dimension + ctd.control_dimension + ctd.variable_dimension)
            error("vector for initialization should be of size dim_x + dim_u + dim_v ie:",ctd.state_dimension+ctd.control_dimension+ctd.variable_dimension)
        end
        x_init = zeros(ctd.dim_NLP_state)
        x_init[1:ctd.state_dimension] = init[1:ctd.state_dimension]
        u_init = zeros(ctd.control_dimension)
        u_init[1:ctd.control_dimension] = init[ctd.state_dimension+1:ctd.state_dimension+ctd.control_dimension]

        # mayer -> lagrange additional state
        if ctd.has_lagrange_cost
            x_init[ctd.dim_NLP_state] = 0.1
        end

        # set constant initialization for state / control variables
        xu0 = zeros(ctd.dim_NLP_variables)
        for i in 0:N
            set_state_at_time_step!(xu0, x_init, ctd, i)
            set_control_at_time_step!(xu0, u_init, ctd, i)
        end

        # set variables
        if (ctd.variable_dimension > 0)
            v_init = zeros(ctd.control_variable)
            v_init[1:ctd.variable_dimension] = init[ctd.state_dimension+ctd.control_dimension+1:ctd.state_dimension+ctd.control_dimension+ctd.variable_dimension]    
            set_variable!(xu0, v_init, ctd)
        end
    end

    return xu0
end
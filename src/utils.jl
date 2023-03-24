function get_state_at_time_step(xu, i, nx, N)
    """
        return
        x(t_i)
    """
    @assert i <= N "trying to get x(t_i) for i > N"
    return xu[i*nx + 1 : (i+1)*nx]
end

function get_control_at_time_step(xu, i, nx, N, m)
    """
        return
        u(t_i)
    """
    @assert i <= N "trying to get u(t_i) for i > N"
    return xu[(N+1)*nx + i*m + 1 : (N+1)*nx + (i+1)*m]
end

get_final_time(xu, fixed_final_time, has_free_final_time) = has_free_final_time ? xu[end] : fixed_final_time

## Initialization for the NLP problem

function set_state_at_time_step!(x, i, nx, N, xu)
    @assert i <= N "trying to set init for x(t_i) with i > N"
    xu[1+i*nx:(i+1)*nx] = x[1:nx]
end
    
function set_control_at_time_step!(u, i, nx, N, m, xu)
    @assert i <= N "trying to set init for u(t_i) with i > N"
    xu[1+(N+1)*nx+i*m:m+(N+1)*nx+i*m] = u[1:m]
end

function initial_guess(ctd)

    N = ctd.dim_NLP_steps
    init = ctd.NLP_init

    if init === nothing
        # default initialization (put back O.1 here ?)
        xu0 = 1.1*ones(ctd.dim_NLP_variables)
    else
        if length(init) != (ctd.state_dimension + ctd.control_dimension)
            error("vector for initialization should be of size n+m",ctd.state_dimension+ctd.control_dimension)
        end
        # split state / control init values
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
            set_state_at_time_step!(x_init, i, ctd.dim_NLP_state, N, xu0)
            set_control_at_time_step!(u_init, i, ctd.dim_NLP_state, N, ctd.control_dimension, xu0)
        end
    end

    # free final time case, put back 0.1 here ?
    # +++todo: add a component in init vector for tf and put this part in main if/then above
    if ctd.has_free_final_time
        xu0[end] = 1.0
    end

    return xu0
end
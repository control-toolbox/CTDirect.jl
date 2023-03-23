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

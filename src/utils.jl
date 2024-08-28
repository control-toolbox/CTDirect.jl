#= Functions that are not dependent on the discretization scheme
Internal layout for NLP variables: 
- Trapeze: [X_0,U_0, X_1,U_1, .., X_N,U_N, V]
- Midpoint: [X_0,U_0,K_0, X_1,U_1,K_1 .., X_N-1,U_N-1,K_N-1, XN, V]
with the convention u([t_i,t_i+1[) = U_i and u(tf) = U_N-1
=#


"""
$(TYPEDSIGNATURES)

Retrieve optimization variables from the NLP variables.
"""
function get_optim_variable(xu, docp)

    if docp.has_variable
        if docp.dim_NLP_v == 1
            #return xu[offset + 1]
            return xu[end]
        else
            #return xu[(offset + 1):(offset + docp.dim_NLP_v)]
            return xu[(end - docp.dim_NLP_v + 1):end]
        end
    else
        return similar(xu, 0)
    end
end


"""
$(TYPEDSIGNATURES)

Retrieve initial time for OCP (may be fixed or variable)
"""
function get_initial_time(xu, docp)
    if docp.has_free_t0
        return get_optim_variable(xu, docp)[docp.ocp.initial_time]
    else
        return docp.ocp.initial_time
    end
end


"""
$(TYPEDSIGNATURES)

Retrieve final time for OCP (may be fixed or variable)
"""
function get_final_time(xu, docp)
    if docp.has_free_tf
        return get_optim_variable(xu, docp)[docp.ocp.final_time]
    else
        return docp.ocp.final_time
    end
end


"""
$(TYPEDSIGNATURES)

Get actual (un-normalized) time value
"""
function get_unnormalized_time(xu, docp, t_normalized)
    t0 = get_initial_time(xu, docp)
    tf = get_final_time(xu, docp)
    return t0 + t_normalized * (tf - t0)
end


"""
$(TYPEDSIGNATURES)

Get actual (un-normalized) time at give time step
"""
function get_time_at_time_step(xu, docp, i)
    return get_unnormalized_time(xu, docp, docp.NLP_normalized_time_grid[i + 1])
end


"""
$(TYPEDSIGNATURES)

Set optimization variables in the NLP variables (for initial guess)
"""
function set_optim_variable!(xu, v_init, docp)
    xu[(end - docp.dim_NLP_v + 1):end] .= v_init
end


"""
$(TYPEDSIGNATURES)

Build full, ordered sets of bounds for state, control or optimization variables
"""
function build_bounds(dim_var, dim_box, box_triplet)

    x_lb = -Inf * ones(dim_var)
    x_ub = Inf * ones(dim_var)
    for j = 1:(dim_box)
        indice = box_triplet[2][j]
        x_lb[indice] = box_triplet[1][j]
        x_ub[indice] = box_triplet[3][j]
    end

    return x_lb, x_ub
end


"""
$(TYPEDSIGNATURES)

Build initial guess for discretized problem
"""
function DOCP_initial_guess(docp, init::OptimalControlInit = OptimalControlInit())

    # default initialization (internal variables such as lagrange cost, k_i for RK schemes) will keep these default values 
    NLP_X = 0.1 * ones(docp.dim_NLP_variables)

    # set variables if provided (needed first in case of free times !)
    if !isnothing(init.variable_init)
        set_optim_variable!(NLP_X, init.variable_init, docp)
    end

    # set state / control variables if provided
    for i = 0:(docp.dim_NLP_steps)
        ti = get_time_at_time_step(NLP_X, docp, i)
        set_variables_at_time_step!(NLP_X, init.state_init(ti), init.control_init(ti), docp, i)
    end

    return NLP_X
end

# placeholders (see CTDirectExt)
function export_ocp_solution end
function import_ocp_solution end

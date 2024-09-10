#= Functions that are not dependent on the discretization scheme
Internal layout for NLP variables: 
- Trapeze: [X_0,U_0, X_1,U_1, .., X_N,U_N, V]
- Midpoint: [X_0,U_0,K_0, X_1,U_1,K_1 .., X_N-1,U_N-1,K_N-1, XN, V]
with the convention u([t_i,t_i+1[) = U_i and u(tf) = U_N-1
=#

# +++ todo change arguments order: docp first, then xu

function vectorize(fun, dim_x, dim_u)
    if dim_x == 1
        if dim_u == 1
            fun2 = (t, x, u, v) -> fun(t, x[1], u[1], v)
        else
            fun2 = (t, x, u, v) -> fun(t, x[1], u, v)
        end
    else
        if dim_u == 1
            fun2 = (t, x, u, v) -> fun(t, x[1:dim_x], u[1], v)
        else
            fun2 = (t, x, u, v) -> fun(t, x[1:dim_x], u, v)
        end
    end
    if length(fun2) == 1
        return [fun2]
    else
        return fun2
    end
end


"""
$(TYPEDSIGNATURES)

Retrieve optimization variables from the NLP variables.
"""
function get_optim_variable(xu, docp)

    if docp.has_variable
        if docp.dim_NLP_v == 1
            return xu[end]
        else
            return xu[(end - docp.dim_NLP_v + 1):end]
        end
    else
        error("Problem is not variable dependent")
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

Get full (un-normalized) time grid
"""
function get_time_grid(xu, docp)
    t0 = get_initial_time(xu, docp)
    tf = get_final_time(xu, docp)
    return @. t0 + docp.NLP_normalized_time_grid * (tf - t0)
end
function get_time_grid!(time_grid, xu, docp)
    t0 = get_initial_time(xu, docp)
    tf = get_final_time(xu, docp)
    @. time_grid = t0 + docp.NLP_normalized_time_grid * (tf - t0)
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


# placeholders (see CTDirectExt) +++ can be removed if functions moved to ctbase
function export_ocp_solution end
function import_ocp_solution end

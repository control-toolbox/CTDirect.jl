#= Functions that are not dependent on the discretization scheme
Internal layout for NLP variables: 
- Trapeze: [X_0,U_0, X_1,U_1, .., X_N,U_N, V]
- Midpoint: [X_0,U_0,K_0, X_1,U_1,K_1 .., X_N-1,U_N-1,K_N-1, XN, V]
with the convention u([t_i,t_i+1[) = U_i and u(tf) = U_N-1
=#

# getter for optimization variables
function get_optim_variable(xu, docp)
    if docp.has_variable
        if docp.dim_NLP_v == 1
            return xu[end]
        else
            return @view xu[(end - docp.dim_NLP_v + 1):end]
        end
    else
        return Float64[]
    end
end

# getters for initial and final time
function get_initial_time(xu, docp)
    if 
    return get_initial_time(xu, docp, docp.ocp.initial_time)
end
function get_initial_time(xu, docp, ti::Real)
    return ti
end
function get_initial_time(xu, docp, ti_index::Index)
    return get_optim_variable(xu, docp)[ti_index]
end

function get_final_time(xu, docp)
    return get_final_time(xu, docp, docp.ocp.final_time)           
end
function get_final_time(xu, docp, tf::Real)
    return tf
end
function get_final_time(xu, docp, tf_index::Index)
    return get_optim_variable(xu, docp)[tf_index]
end

# time grid
function get_time_grid!(xu, docp)
    t0 = get_initial_time(xu, docp)
    tf = get_final_time(xu, docp)
    @. docp.NLP_time_grid = t0 + docp.NLP_normalized_time_grid * (tf - t0)
    return
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

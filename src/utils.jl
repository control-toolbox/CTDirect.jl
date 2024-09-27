#= Functions that are not dependent on the discretization scheme
Internal layout for NLP variables: 
- Trapeze: [X_0,U_0, X_1,U_1, .., X_N,U_N, V]
- Midpoint: [X_0,U_0,K_0, X_1,U_1,K_1 .., X_N-1,U_N-1,K_N-1, XN, V]
with the convention u([t_i,t_i+1[) = U_i and u(tf) = U_N-1
=#

# getter for optimization variables
function get_optim_variable(xu, docp)
    # empty for dim 0
    return @view xu[(end - docp.dim_NLP_v + 1):end]
end

function get_OCP_variable(xu, docp)
    # empty for dim 0
    if docp.dim_NLP_v == 1
        return xu[end]
    else
        return @view xu[(end - docp.dim_NLP_v + 1):end]
    end
end
#=function get_OCP_variable2(xu, docp, dim_v::Val{1})
    return xu[end]
end
function get_OCP_variable2(xu, docp, dim_v)
    return @view xu[(end - docp.dim_NLP_v + 1):end]
end=#

# +++ in problem.jl ?
# getters for initial and final time
function get_initial_time(xu, docp)
    if docp.is_free_initial_time
        return get_optim_variable(xu, docp)[docp.index_initial_time]
    else
        return docp.fixed_initial_time
    end
end

function get_final_time(xu, docp)
    if docp.is_free_final_time
        return get_optim_variable(xu, docp)[docp.index_final_time]
    else
        return docp.fixed_final_time
    end
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

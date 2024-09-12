#= Functions that are not dependent on the discretization scheme
Internal layout for NLP variables: 
- Trapeze: [X_0,U_0, X_1,U_1, .., X_N,U_N, V]
- Midpoint: [X_0,U_0,K_0, X_1,U_1,K_1 .., X_N-1,U_N-1,K_N-1, XN, V]
with the convention u([t_i,t_i+1[) = U_i and u(tf) = U_N-1
=#

# +++ todo change arguments order: docp first, then xu ?


# Vectorize OCP functions (in and out arguments)
#= NB try the simpler ?
_vec(x::Number) = [x]
_vec(x::AbstractVector) = x

function f(x, y, z)  # could specify ::Union{Number, AbstractVector}
  xv = _vec(x)
  yv = _vec(y)
  ...
end=#

# +++ add v case v[1] if dim1, v else (including empty v case)
# vectorization for state or control constraints
# NB. result is vectorized
function vectorize_1(fun, dim_arg, dim_f)
    if dim_arg == 1
        funv = (x1, x2, v) -> [fun(x1[1], x2[1], v)]
    else
        funv = fun
    end
    if dim_f == 1
        return (x1, x2, v) -> [funv(x1, x2, v)]
    else
        return funv
    end
end
# vectorization for mayer cost 
# NB. keep scalar result
function vectorize_xx(fun, dim_x)
    if dim_x == 1
        return (x1, x2, v) -> [fun(x1[1], x2[1], v)]
    else
        return fun
    end
end
# vectorization for dynamics, lagrange cost, mixed constraints
# NB. result is vectorized
function vectorize_xu(fun, dim_x, dim_u, dim_f)
    if dim_x == 1
        if dim_u == 1
            funv = (t, x, u, v) -> fun(t, x[1], u[1], v)
        else
            funv = (t, x, u, v) -> fun(t, x[1], u, v)
        end
    else
        if dim_u == 1
            funv = (t, x, u, v) -> fun(t, x[1:dim_x], u[1], v)
        else
            funv = (t, x, u, v) -> fun(t, x[1:dim_x], u, v)
        end
    end
    if dim_f == 1
        return (t, x, u, v) -> [funv(t, x, u, v)]
    else
        return funv
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

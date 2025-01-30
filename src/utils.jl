
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


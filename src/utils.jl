function get_state_at_time_step(xu, i, dim_x, N)
    """
        return
        x(t_i)
    """
    if i > N
        error("trying to get x(t_i) for i > N")
    end  
    return xu[1+i*dim_x:(i+1)*dim_x]
end

function get_control_at_time_step(xu, i, dim_x, N, m)
    """
        return
        u(t_i)
    """
    if i > N
        error("trying to get (t_i) for i > N")
    end
    return xu[1+(N+1)*dim_x+i*m:m+(N+1)*dim_x+i*m]
end

get_final_time(xu, tf_, has_free_final_time) = has_free_final_time ? xu[end] : tf_

function direct_infos(ocp::OptimalControlModel, N::Integer)

    # Parameters of the Optimal Control Problem
    # times
    t0 = ocp.initial_time
    tf = ocp.final_time
    has_free_final_time = isnothing(tf)
        # multiplier la dynamique par (tf-t0)
        # travailler avec le nouveau temps s dans (0., 1.)
        # une fonction t(s)
    # dimensions
    n_x = ocp.state_dimension
    m = ocp.control_dimension
    # dynamics
    f = ocp.dynamics
    # constraints

    control_constraints, state_constraints, mixed_constraints, boundary_conditions, control_box, state_box = nlp_constraints(ocp)
    println("control_box = ", control_box)
    dim_control_constraints = length(control_constraints[1])      # dimension of the constraints
    dim_state_constraints = length(state_constraints[1])
    dim_mixed_constraints = length(mixed_constraints[1])
    dim_boundary_conditions = length(boundary_conditions[1])
    dim_control_box = length(control_box[1])
    dim_state_box = length(state_box[1])
    has_control_constraints = !isempty(control_constraints[1])
    has_state_constraints = !isempty(state_constraints[1])
    has_mixed_constraints = !isempty(mixed_constraints[1])
    has_boundary_conditions = !isempty(boundary_conditions[1])
    has_control_box = !isempty(control_box[1])
    has_state_box = !isempty(state_box[1])

    #println("has_control_constraints = ", has_control_constraints)
    #println("has_mixed_constraints = ", has_mixed_constraints)
    #println("has_boundary_conditions = ", has_boundary_conditions)

    hasLagrangeCost = !isnothing(ocp.lagrange)
    L = ocp.lagrange

    hasMayerCost = !isnothing(ocp.mayer)
    g = ocp.mayer
    #println("hasLagrange : ", hasLagrangeCost)
    #println("Mayer = ", hasMayerCost)

    # Mayer formulation
    # use an additional state for the Lagrange cost
    #
    # remark : we pass u[1] because in our case ocp.dynamics and ocp.lagrange are defined with a scalar u
    # and we consider vectors for x and u in the discretized problem. Note that the same would apply for a scalar x.
    # question : how determine if u and x are scalar or vector ?
    # second member of the ode for the Mayer formulation

    if hasLagrangeCost
        dim_x = n_x + 1  
        nc = N*(dim_x + dim_control_constraints + dim_state_constraints + dim_mixed_constraints) +
        (dim_control_constraints +  dim_state_constraints + dim_mixed_constraints) + dim_boundary_conditions + 1       # dimension of the constraints            
    else
        dim_x = n_x  
        nc = N*(dim_x + dim_control_constraints + dim_state_constraints + dim_mixed_constraints) +
        (dim_control_constraints + dim_state_constraints + dim_mixed_constraints) + dim_boundary_conditions       # dimension of the constraints
    end

    dim_xu = (N+1)*(dim_x+m)  # dimension the the unknown xu
    has_free_final_time ? dim_xu = dim_xu + 1 : nothing

    # todo: cas vectoriel sur u a ajouter
    f_Mayer(t, x, u) = hasLagrangeCost ? [f(t, x[1:n_x], u); L(t, x[1:n_x], u)] : f(t, x, u)

    criterion = ocp.criterion

    return t0, tf, n_x, m, f, control_constraints, state_constraints, mixed_constraints, boundary_conditions, control_box, state_box, 
    dim_control_constraints, dim_state_constraints, dim_mixed_constraints, dim_boundary_conditions, dim_control_box, dim_state_box,
    has_control_constraints, has_state_constraints, has_mixed_constraints, has_boundary_conditions,
    has_control_box, has_state_box,
    hasLagrangeCost, hasMayerCost, dim_x, nc, dim_xu, 
    g, f_Mayer, has_free_final_time, criterion

end

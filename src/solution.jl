# build generic OCP solution from direct method (NLP) solution
function _OptimalControlSolution(ocp, ipopt_solution, ctd)

    # save general solution data
    ctd.NLP_stats = ipopt_solution
    if is_min(ocp)
        ctd.NLP_objective = ipopt_solution.objective
    else
        ctd.NLP_objective = - ipopt_solution.objective
    end
    ctd.NLP_constraints_violation = ipopt_solution.primal_feas
    ctd.NLP_iterations = ipopt_solution.iter
    ctd.NLP_solution = ipopt_solution.solution
    ctd.NLP_sol_constraints = ipopt_constraint(ipopt_solution.solution, ctd)

    # parse NLP variables, constraints and multipliers 
    X, U, v, P, sol_control_constraints, sol_state_constraints, sol_mixed_constraints, sol_variable_constraints, mult_control_constraints, mult_state_constraints, mult_mixed_constraints, mult_variable_constraints, mult_state_box_lower, mult_state_box_upper, mult_control_box_lower, mult_control_box_upper, mult_variable_box_lower, variable_box_upper = parse_ipopt_sol(ctd)

    # variables and misc infos
    N = ctd.dim_NLP_steps
    t0 = get_initial_time(ctd.NLP_solution, ctd)
    tf = get_final_time(ctd.NLP_solution, ctd)
    T = collect(LinRange(t0, tf, N+1))
    x = ctinterpolate(T, matrix2vec(X, 1))
    u = ctinterpolate(T, matrix2vec(U, 1))
    p = ctinterpolate(T[1:end-1], matrix2vec(P, 1))
    sol = OptimalControlSolution() # +++ constructor with ocp as argument ?
    copy!(sol,ocp)
    sol.times = T
    sol.state = t -> x(t)
    sol.costate = t -> p(t)
    sol.control = t -> u(t)
    #sol.variable = v +++no field variable?
    sol.infos[:variable] = v
    sol.objective = ctd.NLP_objective
    sol.iterations = ctd.NLP_iterations
    sol.stopping = :dummy 
    sol.message = "no message" 
    sol.success = false #

    # nonlinear constraints and multipliers
    if ctd.has_state_constraints
        cx = ctinterpolate(T, matrix2vec(sol_state_constraints, 1))
        mcx = ctinterpolate(T, matrix2vec(mult_state_constraints, 1))
        sol.infos[:dim_state_constraints] = ctd.dim_state_constraints    
        sol.infos[:state_constraints] = t -> cx(t)
        sol.infos[:mult_state_constraints] = t -> mcx(t)
    end
    if ctd.has_control_constraints
        cu = ctinterpolate(T, matrix2vec(sol_control_constraints, 1))
        mcu = ctinterpolate(T, matrix2vec(mult_control_constraints, 1))
        sol.infos[:dim_control_constraints] = ctd.dim_control_constraints  
        sol.infos[:control_constraints] = t -> cu(t)
        sol.infos[:mult_control_constraints] = t -> mcu(t)
    end
    if ctd.has_mixed_constraints
        cxu = ctinterpolate(T, matrix2vec(sol_mixed_constraints, 1))
        mcxu = ctinterpolate(T, matrix2vec(mult_mixed_constraints, 1))
        sol.infos[:dim_mixed_constraints] = ctd.dim_mixed_constraints    
        sol.infos[:mixed_constraints] = t -> cxu(t)
        sol.infos[:mult_mixed_constraints] = t -> mcxu(t)
    end
    if ctd.has_variable_constraints
        sol.infos[:dim_variable_constraints] = ctd.dim_variable_constraints
        sol.infos[:variable_constraints] = sol_variable_constraints
        sol.infos[:mult_variable_constraints] = mult_variable_constraints
    end

    # box constraints multipliers
    if ctd.has_state_box
        mbox_x_l = ctinterpolate(T, matrix2vec(mult_state_box_lower, 1))
        mbox_x_u = ctinterpolate(T, matrix2vec(mult_state_box_upper, 1))
        sol.infos[:mult_state_box_lower] = t -> mbox_x_l(t)
        sol.infos[:mult_state_box_upper] = t -> mbox_x_u(t)    
    end
    if ctd.has_control_box
        mbox_u_l = ctinterpolate(T, matrix2vec(mult_control_box_lower, 1))
        mbox_u_u = ctinterpolate(T, matrix2vec(mult_control_box_upper, 1))
        sol.infos[:mult_control_box_lower] = t -> mbox_u_l(t)
        sol.infos[:mult_control_box_upper] = t -> mbox_u_u(t)
    end
    if ctd.has_variable_box
        sol.infos[:mult_variable_box_lower] = mult_variable_box_lower
        sol.infos[:mult_variable_box_upper] = mult_variable_box_upper 
    end

    return sol

end


# parse NLP solution from ipopt into OCP variables, constraints and multipliers
function parse_ipopt_sol(ctd)
    
    N = ctd.dim_NLP_steps

    # states and controls variables, with box multipliers
    xu = ctd.NLP_solution
    mult_L = ctd.NLP_stats.multipliers_L
    mult_U = ctd.NLP_stats.multipliers_U
    X = zeros(N+1,ctd.dim_NLP_state)
    U = zeros(N+1,ctd.control_dimension)
    v = get variables(xu, ctd)
    mult_state_box_lower = zeros(N+1,ctd.dim_NLP_state)
    mult_state_box_upper = zeros(N+1,ctd.dim_NLP_state)
    mult_control_box_lower = zeros(N+1,ctd.control_dimension)
    mult_control_box_upper = zeros(N+1,ctd.control_dimension)
    mult_variable_box_lower = zeros(N+1,ctd.variable_dimension)
    mult_variable_box_upper = zeros(N+1,ctd.variable_dimension)

    # +++ v box + get
    for i in 1:N+1
        # variables
        X[i,:] = vget_state_at_time_step(xu, ctd, i-1)
        U[i,:] = vget_control_at_time_step(xu, ctd, i-1)
        # box multipliers (same layout as variables !)
        # +++ fix mult vectors instead to always be full size (todo in return from solve)
        if length(mult_L) > 0
            mult_state_box_lower[i,:] = vget_state_at_time_step(mult_L, ctd, i-1)
            mult_control_box_lower[i,:] = vget_control_at_time_step(mult_L, ctd, i-1)
        end
        if length(mult_U) > 0
            mult_state_box_upper[i,:] = vget_state_at_time_step(mult_U, ctd, i-1)
            mult_control_box_upper[i,:] = vget_control_at_time_step(mult_U, ctd, i-1)
        end
    end
    if ctd.has_variable_box
        if length(mult_L) > 0
            mult_variable_box_lower = get_variable(mult_L, ctd)
        end
        if length(mult_U) > 0
            mult_variable_box_upper = get_variable(mult_U, ctd)
        end
    end

    # constraints, costate and constraints multipliers
    P = zeros(N, ctd.dim_NLP_state)
    lambda = ctd.NLP_stats.multipliers
    c = ctd.NLP_sol_constraints
    sol_control_constraints = zeros(N+1,ctd.dim_control_constraints)
    sol_state_constraints = zeros(N+1,ctd.dim_state_constraints)
    sol_mixed_constraints = zeros(N+1,ctd.dim_mixed_constraints) 
    sol_variable_constraints = zeros(ctd.dim_variable_constraints)
    mult_control_constraints = zeros(N+1,ctd.dim_control_constraints)
    mult_state_constraints = zeros(N+1,ctd.dim_state_constraints)
    mult_mixed_constraints = zeros(N+1,ctd.dim_mixed_constraints)
    mult_variable_constraints = zeros(ctd.dim_variable_constraints)

    index = 1
    for i in 1:N
        # state equation
        P[i,:] = lambda[index:index+ctd.dim_NLP_state-1]
        index = index + ctd.dim_NLP_state
        # path constraints
        if ctd.has_control_constraints
            sol_control_constraints[i,:] = c[index:index+ctd.dim_control_constraints-1]
            mult_control_constraints[i,:] = lambda[index:index+ctd.dim_control_constraints-1]
            index = index + ctd.dim_control_constraints
        end
        if ctd.has_state_constraints
            sol_state_constraints[i,:] = c[index:index+ctd.dim_state_constraints-1]
            mult_state_constraints[i,:] = lambda[index:index+ctd.dim_state_constraints-1]
            index = index + ctd.dim_state_constraints
        end
        if ctd.has_mixed_constraints
            sol_mixed_constraints[i,:] = c[index:index+ctd.dim_mixed_constraints-1]
            mult_mixed_constraints[i,:] = lambda[index:index+ctd.dim_mixed_constraints-1]
            index = index + ctd.dim_mixed_constraints
        end
    end
    # path constraints at final time
    if ctd.has_control_constraints
        sol_control_constraints[N+1,:] = c[index:index+ctd.dim_control_constraints-1]
        mult_control_constraints[N+1,:] = lambda[index:index+ctd.dim_control_constraints-1]
        index = index + ctd.dim_control_constraints
    end
    if ctd.has_state_constraints
        sol_state_constraints[N+1,:] = c[index:index+ctd.dim_state_constraints-1] 
        mult_state_constraints[N+1,:] = lambda[index:index+ctd.dim_state_constraints-1]
        index = index + ctd.dim_state_constraints
    end
    if ctd.has_mixed_constraints
        sol_mixed_constraints[N+1,:] = c[index:index+ctd.dim_mixed_constraints-1]        
        mult_mixed_constraints[N+1,:] =  lambda[index:index+ctd.dim_mixed_constraints-1]
        index = index + ctd.dim_mixed_constraints
    end

    # boundary conditions and multipliers
    if ctd.has_boundary_conditions
        sol_boundary_conditions = c[index:index+ctd.dim_boundary_conditions-1]
        mult_boundary_conditions = lambda[index:index+ctd.dim_boundary_conditions-1]
        index = index + ctd.dim_boundary_conditions
    end

    # variable constraints and multipliers
    if ctd.has_variable_constraints
        sol_variable_constraints = c[index:index+ctd.dim_variable_constraints-1]
        mult_variable_constraints = lambda[index:index+ctd.dim_variable_constraints-1]
        index = index + ctd.dim_variable_constraints
    end


    return X, U, v, P, sol_control_constraints, sol_state_constraints, sol_mixed_constraints, sol_variable_constraints, mult_control_constraints, mult_state_constraints, mult_mixed_constraints, mult_variable_constraints, mult_state_box_lower, mult_state_box_upper, mult_control_box_lower, mult_control_box_upper, mult_variable_box_lower, mult_variable_box_upper
end

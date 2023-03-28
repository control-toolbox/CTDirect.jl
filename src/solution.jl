function parse_ipopt_sol(stats, ctd)
    
    N = ctd.dim_NLP_steps

    # states and controls
    xu = stats.solution
    X = zeros(N+1,ctd.dim_NLP_state)
    U = zeros(N+1,ctd.control_dimension)
    for i in 1:N+1
        X[i,:] = get_state_at_time_step(xu, i-1, ctd.dim_NLP_state, N)
        U[i,:] = get_control_at_time_step(xu, i-1, ctd.dim_NLP_state, N, ctd.control_dimension)
    end

    # adjoints
    P = zeros(N, ctd.dim_NLP_state)
    lambda = stats.multipliers
    P_control_constraints = zeros(N+1,ctd.dim_control_constraints)
    P_state_constraints = zeros(N+1,ctd.dim_state_constraints)
    P_mixed_constraints = zeros(N+1,ctd.dim_mixed_constraints)
    index = 1
    for i in 1:N
        # state equation
        P[i,:] = lambda[index:index+ctd.dim_NLP_state-1]
        index = index + ctd.dim_NLP_state
        # path constraints
        if ctd.has_control_constraints
            P_control_constraints[i,:] = lambda[index:index+ctd.dim_control_constraints-1]
            index = index + ctd.dim_control_constraints
        end
        if ctd.has_state_constraints
            P_state_constraints[i,:] = lambda[index:index+ctd.dim_state_constraints-1]
            index = index + ctd.dim_state_constraints
        end
        if ctd.has_mixed_constraints
            P_mixed_constraints[i,:] = lambda[index:index+ctd.dim_mixed_constraints-1]
            index = index + ctd.dim_mixed_constraints
        end
    end
    # path constraints at final time
    if ctd.has_control_constraints
        P_control_constraints[N+1,:] = lambda[index:index+ctd.dim_control_constraints-1]
        index = index + ctd.dim_control_constraints
    end
    if ctd.has_state_constraints
        P_state_constraints[N+1,:] = lambda[index:index+ctd.dim_state_constraints-1]
        index = index + ctd.dim_state_constraints
    end
    if ctd.has_mixed_constraints
        P_mixed_constraints[N+1,:] =  lambda[index:index+ctd.dim_mixed_constraints-1]
        index = index + ctd.dim_mixed_constraints
    end

    return X, U, P, P_control_constraints, P_state_constraints, P_mixed_constraints
end


function DirectSolution(ocp, N, ipopt_solution, init)

    ctd = CTDirect_data(ocp, N, init)

    # state, control, adjoint
    X, U, P, P_control_constraints, P_state_constraints, P_mixed_constraints = parse_ipopt_sol(ipopt_solution, ctd)
    
    # times
    t0 = ctd.initial_time
    tf = get_final_time(ipopt_solution.solution, ctd.final_time, ctd.has_free_final_time)
    T = collect(LinRange(t0, tf, N+1))
    
    # misc info
    if ismin(ocp) # +++ remove this once obj_factor s used instead of changing the sign of obj
        objective = ipopt_solution.objective
    else
        objective = - ipopt_solution.objective
    end

    constraints_violation = ipopt_solution.primal_feas
    iterations = ipopt_solution.iter
    #status = ipopt_solution.status this is a 'Symbol' not an int...
        
    # DirectSolution
    # To do add P_state_constraints to DirectSolution
    # and quid constraints and Lagrange multiplayers of the constraints
    #dsol  = DirectSolution(T, X, U, P, P_control_constraints, P_mixed_constraints, ctd.state_dimension, ctd.control_dimension, N, objective, constraints_violation, iterations, ipopt_solution)     

    # save solution data
    ctd.T = T 
    ctd.X = X
    ctd.U = U 
    ctd.P = P 
    ctd.P_control_constraints = P_control_constraints
    ctd.P_mixed_constraints = P_mixed_constraints
    ctd.objective = objective
    ctd.constraints_violation = constraints_violation
    ctd.iterations = iterations
    ctd.stats = ipopt_solution

    return _OptimalControlSolution(ocp, ctd)
end


function _OptimalControlSolution(ocp, ctd)

    # je ne peux pas donner directement la sortie de ctinterpolate car ce n'est pas une Function. CTBase doit etre mis a jour (?)
    # matrix2vec is in CTBase/src/utils.jl
    x = ctinterpolate(ctd.T, matrix2vec(ctd.X, 1))
    u = ctinterpolate(ctd.T, matrix2vec(ctd.U, 1))
    p = ctinterpolate(ctd.T[1:end-1], matrix2vec(ctd.P, 1)) 
    sol = OptimalControlSolution()
    sol.state_dimension = ctd.state_dimension
    sol.control_dimension = ctd.control_dimension
    sol.times = ctd.T
    sol.time_label = ocp.time_label
    sol.state = t -> x(t)
    sol.state_labels = ocp.state_labels # update CTBase to have a getter ?
    sol.adjoint = t -> p(t)
    sol.control = t -> u(t)
    sol.control_labels = ocp.control_labels
    sol.objective = ctd.objective
    sol.iterations = ctd.iterations
    # sync with CTDirectShooting: :optimality, :stagnation, :iterations
    sol.stopping = :dummy 
    # see CTDirectShooting/src/solve.jl : textsStopping
    sol.message = "no message" 
    sol.success = false #
    # field "infos" in OptimalControlSolution is Dict{Symbol, Any}, for misc data
    return sol
end
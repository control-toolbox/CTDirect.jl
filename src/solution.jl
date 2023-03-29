# build generic OCP solution from direct method (NLP) solution
# +++ put ocp inside ctd ?
function _OptimalControlSolution(ocp, ipopt_solution, ctd)

    # save general solution data
    ctd.NLP_stats = ipopt_solution
    ctd.NLP_objective = ipopt_solution.objective
    ctd.NLP_constraints_violation = ipopt_solution.primal_feas
    ctd.NLP_iterations = ipopt_solution.iter
    ctd.NLP_solution = ipopt_solution.solution
    ctd.NLP_sol_constraints = ipopt_constraint(ipopt_solution.solution, ctd)

    # parse NLP variables, constraints and multipliers
    X, U, P, CU, CX, CXU, mult_CU, mult_CX, mult_CXU = parse_ipopt_sol(ctd)

    # interpolations to build functions of time
    N = ctd.dim_NLP_steps
    t0 = ctd.initial_time
    tf = get_final_time(ctd.NLP_solution, ctd.final_time, ctd.has_free_final_time)
    T = collect(LinRange(t0, tf, N+1))
    x = ctinterpolate(T, matrix2vec(X, 1))
    u = ctinterpolate(T, matrix2vec(U, 1))
    p = ctinterpolate(T[1:end-1], matrix2vec(P, 1))
    cx = ctinterpolate(T, matrix2vec(CX, 1))
    mcx = ctinterpolate(T, matrix2vec(mult_CX, 1))

    # build OptimalControlSolution struct
    sol = OptimalControlSolution()
    sol.state_dimension = ctd.state_dimension
    sol.control_dimension = ctd.control_dimension
    sol.times = T
    sol.time_label = ocp.time_label
    sol.state = t -> x(t)
    sol.state_labels = ocp.state_labels
    sol.adjoint = t -> p(t)
    sol.control = t -> u(t)
    sol.control_labels = ocp.control_labels
    sol.objective = ctd.NLP_objective
    sol.iterations = ctd.NLP_iterations
    sol.stopping = :dummy 
    sol.message = "no message" 
    sol.success = false #
    # field "infos" in OptimalControlSolution is Dict{Symbol, Any}, for misc data
    # +++ todo: labels for costraints
    sol.infos[:dim_state_constraints] = ctd.dim_state_constraints    
    sol.infos[:state_constraints] = t -> cx(t)
    sol.infos[:mult_state_constraints] = t -> mcx(t)


    return sol

end


# parse NLP solution from ipopt into OCP variables, constraints and multipliers
function parse_ipopt_sol(ctd)
    
    N = ctd.dim_NLP_steps

    # states and controls
    xu = ctd.NLP_solution
    X = zeros(N+1,ctd.dim_NLP_state)
    U = zeros(N+1,ctd.control_dimension)
    for i in 1:N+1
        X[i,:] = get_state_at_time_step(xu, i-1, ctd.dim_NLP_state, N)
        U[i,:] = get_control_at_time_step(xu, i-1, ctd.dim_NLP_state, N, ctd.control_dimension)
    end

    # constraints, costate and constraints multipliers
    P = zeros(N, ctd.dim_NLP_state)
    lambda = ctd.NLP_stats.multipliers
    c = ctd.NLP_sol_constraints
    sol_control_constraints = zeros(N+1,ctd.dim_control_constraints)
    sol_state_constraints = zeros(N+1,ctd.dim_state_constraints)
    sol_mixed_constraints = zeros(N+1,ctd.dim_mixed_constraints)    
    mult_control_constraints = zeros(N+1,ctd.dim_control_constraints)
    mult_state_constraints = zeros(N+1,ctd.dim_state_constraints)
    mult_mixed_constraints = zeros(N+1,ctd.dim_mixed_constraints)
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

    return X, U, P, sol_control_constraints, sol_state_constraints, sol_mixed_constraints, mult_control_constraints, mult_state_constraints, mult_mixed_constraints
end

mutable struct DirectSolution
    T::Vector{<:MyNumber}
    X::Matrix{<:MyNumber}
    U::Matrix{<:MyNumber}
    P::Matrix{<:MyNumber}
    P_control_constraints::Matrix{<:MyNumber}
    P_mixed_constraints::Matrix{<:MyNumber}
    n::Integer
    m::Integer
    N::Integer
    stats  #stats::SolverCore.GenericExecutionStats
end

function DirectSolution(ocp::OptimalControlModel, N::Integer, ipopt_solution)

    # direct_infos
    t0, tf_, n_x, m, f, control_constraints, mixed_constraints, boundary_conditions, dim_control_constraints, dim_mixed_constraints, dim_boundary_conditions, 
    has_control_constraints, has_mixed_constraints, has_boundary_conditions, hasLagrangeCost, hasMayerCost, 
    dim_x, nc, dim_xu, f_Mayer, has_free_final_time, criterion = direct_infos(ocp, N)

    function parse_ipopt_sol(stats)
        """
            return
            X : matrix(N+1,n+1)
            U : matrix(N,m)
            P : matrix(N,n+1)
        """
        # states and controls
        xu = stats.solution
        X = zeros(N+1,dim_x)
        U = zeros(N,m)
        for i in 1:N
            X[i,:] =  get_state_at_time_step(xu, i-1, dim_x, N)
            U[i,:] = get_control_at_time_step(xu, i-1, dim_x, N, m)
        end
        X[N+1,:] = get_state_at_time_step(xu, N, dim_x, N)

        # adjoints
        P = zeros(N, dim_x)
        lambda = stats.multipliers
        P_control_constraints = zeros(N,dim_control_constraints)
        P_mixed_constraints = zeros(N+1,dim_mixed_constraints)
        index = 1 # counter for the constraints
        for i ∈ 1:N
            # state equation
            P[i,:] = lambda[index:index+dim_x-1]
            index = index + dim_x
            if has_control_constraints
                P_control_constraints[i,:] =  lambda[index:index+dim_control_constraints-1]
                index = index + dim_control_constraints
            end
            if has_mixed_constraints
                P_mixed_constraints[i,:] =  lambda[index:index+dim_mixed_constraints-1]
                index = index + dim_mixed_constraints
            end
            P_mixed_constraints[N+1,:] = lambda[index:index+dim_mixed_constraints-1]
        end
        return X, U, P, P_control_constraints, P_mixed_constraints
    end

    # state, control, adjoint
    X, U, P, P_control_constraints, P_mixed_constraints = parse_ipopt_sol(ipopt_solution)
    
    # times
    if has_free_final_time
        tf = stats.solution[end]
    else
        tf = tf_
    end
    T = collect(t0:(tf-t0)/N:tf)

    # DirectSolution
    sol  = DirectSolution(T, X, U, P, P_control_constraints, P_mixed_constraints, n_x, m, N, ipopt_solution)

    return sol
end

# ------------------------------------------------------------------------------------
# Direct solution
#
struct DirectSolution <: AbstractOptimalControlSolution
    T::Vector{<:MyNumber}
    X::Matrix{<:MyNumber}
    U::Matrix{<:MyNumber}
    P::Matrix{<:MyNumber}
    P_control_constraints::Matrix{<:MyNumber}
    P_mixed_constraints::Matrix{<:MyNumber}
    n::Integer
    m::Integer
    N::Integer
    objective::MyNumber
    constraints_violation::MyNumber
    iterations::Integer
    stats       # remove later ?
    #type is https://juliasmoothoptimizers.github.io/SolverCore.jl/stable/reference/#SolverCore.GenericExecutionStats
end

# getters
# todo: return all variables on a common time grid
# trapeze scheme case: state and control on all time steps [t0,...,tN], as well as path constraints
# only exception is the costate, associated to the dynamics equality constraints: N values instead of N+1
# we could use the basic extension for the final costate P_N := P_(N-1)  (or linear extrapolation) 
# NB. things will get more complicated with other discretization schemes (time stages on a different grid ...)
# NEED TO CHOOSE A COMMON OUTPUT GRID FOR ALL SCHEMES (note that building the output can be scheme-dependent)
# PM: I propose to use the time steps [t0, ... , t_N]
# - states are ok, and we can evaluate boundary conditions directly
# - control are technically on a different grid (stages) that can coincide with the steps for some schemes.
# Alternately, the averaged control (in the sense of the butcher coefficients) can be computed on each step. cf bocop
# - adjoint are for each equation linking x(t_i+1) and x(t_i), so always N values regardless of the scheme

#= state_dimension(sol::DirectSolution) = sol.n
control_dimension(sol::DirectSolution) = sol.m
time_steps_length(sol::DirectSolution) = sol.N
time_steps(sol::DirectSolution) = sol.T            
state(sol::DirectSolution) = sol.X
control(sol::DirectSolution) = sol.U
function adjoint(sol::DirectSolution)
    N = sol.N
    n = sol.n
    P = zeros(N+1, n)
    P[1:N,1:n] = sol.P[1:N,1:n]
    # trivial constant extrapolation for p(t_f)
    P[N+1,1:n] = P[N,1:n]
    return P
end
objective(sol::DirectSolution) = sol.objective
constraints_violation(sol::DirectSolution) = sol.constraints_violation  
iterations(sol::DirectSolution) = sol.iterations =#

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
    objective = ipopt_solution.objective
    constraints_violation = ipopt_solution.primal_feas
    iterations = ipopt_solution.iter
    #status = ipopt_solution.status this is a 'Symbol' not an int...
        
    # DirectSolution
    # To do add P_state_constraints to DirectSolution
    # and quid constraints and Lagrange multiplayers of the constraints
    dsol  = DirectSolution(T, X, U, P, P_control_constraints, P_mixed_constraints, ctd.state_dimension, ctd.control_dimension, N, objective, constraints_violation, iterations, ipopt_solution)     

    return _OptimalControlSolution(ocp, dsol)
end


function _OptimalControlSolution(ocp::OptimalControlModel, dsol::DirectSolution)

    # je ne peux pas donner directement la sortie de ctinterpolate car ce n'est pas une Function. CTBase doit etre mis a jour (?)
    # matrix2vec is in CTBase/src/utils.jl
    x = ctinterpolate(dsol.T, matrix2vec(dsol.X, 1))
    u = ctinterpolate(dsol.T, matrix2vec(dsol.U, 1))
    p = ctinterpolate(dsol.T[1:end-1], matrix2vec(dsol.P, 1)) 
    sol = OptimalControlSolution()
    sol.state_dimension = dsol.n
    sol.control_dimension = dsol.m
    sol.times = dsol.T
    sol.time_label = ocp.time_label
    sol.state = t -> x(t)
    sol.state_labels = ocp.state_labels # update CTBase to have a getter ?
    sol.adjoint = t -> p(t)
    sol.control = t -> u(t)
    sol.control_labels = ocp.control_labels
    sol.objective = dsol.objective
    sol.iterations = dsol.iterations
    # sync with CTDirectShooting: :optimality, :stagnation, :iterations
    sol.stopping = :dummy 
    # see CTDirectShooting/src/solve.jl : textsStopping
    sol.message = "no message" 
    sol.success = false #
    # field "infos" in OptimalControlSolution is Dict{Symbol, Any}, for misc data
    return sol
end
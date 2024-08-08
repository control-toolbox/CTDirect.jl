#+++todo:
# redo constraints/multipliers parsing 
# rewrite 3rd parser for box constraints (use tuples)
# and use proper fields in solution for constraints/multipliers

"""
$(TYPEDSIGNATURES)

Recover OCP primal variables from DOCP solution
"""
function parse_DOCP_solution_primal(docp, solution)

    ocp = docp.ocp

    # recover optimization variables
    v = get_variable(solution, docp)

    # recover states and controls variables
    N = docp.dim_NLP_steps
    X = zeros(N+1,docp.dim_NLP_x)
    U = zeros(N+1,docp.dim_NLP_u)
    for i in 1:N+1
        # state and control variables
        xi, ui = get_NLP_variables_at_time_step(solution, docp, i-1)
        X[i,:] .= xi
        U[i,:] .= ui
    end

    return X, U, v
end


"""
$(TYPEDSIGNATURES)

Recover OCP costate from DOCP multipliers
"""
function parse_DOCP_solution_dual(docp, multipliers)

    # constraints, costate and constraints multipliers
    N = docp.dim_NLP_steps
    P = zeros(N, docp.dim_NLP_x)
    
    if !isnothing(multipliers)
        index = 1
        for i in 1:N
            # state equation multiplier for costate
            P[i,:] = multipliers[index:index+docp.dim_NLP_x-1]
            index = index + docp.dim_NLP_x
            # skip path constraints multipliers
            index = index + dim_path_constraints(docp.ocp)
        end
    end

    return P
end

"""
$(TYPEDSIGNATURES)

Build OCP functional solution from the DOCP discrete solution, given as a vector. Costate will be retrieved from dual variables (multipliers) if available.
"""
function CTBase.OptimalControlSolution(docp; primal, dual=nothing)

    # time grid
    N = docp.dim_NLP_steps
    T = zeros(N+1)
    for i=1:N+1
        T[i] = get_unnormalized_time(primal, docp, docp.NLP_normalized_time_grid[i])
    end

    # recover primal variables
    X, U, v = parse_DOCP_solution_primal(docp, primal)

    # recover costate
    P = parse_DOCP_solution_dual(docp, dual)

    # recompute objective
    objective = DOCP_objective(primal, docp)

    return OCPSolutionFromDOCP_raw(docp, T, X, U, v, P, objective=objective)
end

# dummy (remove ?)
function CTBase.OptimalControlSolution(docp, docp_solution::Nothing)
    return nothing
end

"""
$(TYPEDSIGNATURES)
   
Build OCP functional solution from DOCP discrete solution (given as a SolverCore.GenericExecutionStats)
"""
function CTBase.OptimalControlSolution(docp, docp_solution_ipopt)

    # could pass some status info too (get_status ?)
    solution = docp_solution_ipopt.solution

    # time grid
    N = docp.dim_NLP_steps
    T = zeros(N+1)
    for i=1:N+1
        T[i] = get_unnormalized_time(solution, docp, docp.NLP_normalized_time_grid[i])
    end

    # adjust objective sign for maximization problems
    if is_min(docp.ocp)
        objective = docp_solution_ipopt.objective
    else        
        objective = - docp_solution_ipopt.objective
    end

    # recover primal variables
    X, U, v = parse_DOCP_solution_primal(docp, solution)

    # recover costate
    P = parse_DOCP_solution_dual(docp, docp_solution_ipopt.multipliers)

    # build and return OCP solution
    return OCPSolutionFromDOCP_raw(docp, T, X, U, v, P,
    objective=objective, iterations=docp_solution_ipopt.iter,constraints_violation=docp_solution_ipopt.primal_feas, 
    message=String(docp_solution_ipopt.solver_specific[:internal_msg]))
end


# to be updated
"""
$(TYPEDSIGNATURES)
    
Build OCP functional solution from DOCP vector solution (given as raw variables and multipliers plus some optional infos)
"""
# +++ try to reuse this for the discrete json solution !
# USE SEVERAL METHODS DEPENDING ON AVAILABLE INFO !
# rename as OCS constructor also ?
function OCPSolutionFromDOCP_raw(docp, T, X, U, v, P;
    objective=0, iterations=0, constraints_violation=0,
    message="No msg", stopping=nothing, success=nothing,
    constraints_types=nothing, constraints_mult=nothing,
    box_multipliers=nothing)

    ocp = docp.ocp
    dim_x = ocp.state_dimension()
    dim_u = ocp.control_dimension()
    dim_v = ocp.variable_dimension()

    # check that time grid is strictly increasing
    # if not proceed with list of indexes as time grid
    if !issorted(T,lt=<=)
        println("WARNING: time grid at solution is not strictly increasing, replacing with list of indices...")
        println(T)
        T = LinRange(0,docp.dim_NLP_steps,docp.dim_NLP_steps+1)
    end

    # variables: remove additional state for lagrange cost
    x = ctinterpolate(T, matrix2vec(X[:,1:dim_x], 1))
    p = ctinterpolate(T[1:end-1], matrix2vec(P[:,1:dim_x], 1))
    u = ctinterpolate(T, matrix2vec(U, 1))
    
    # force scalar output when dimension is 1
    fx = (dim_x==1) ? deepcopy(t->x(t)[1]) : deepcopy(t->x(t))
    fu = (dim_u==1) ? deepcopy(t->u(t)[1]) : deepcopy(t->u(t))
    fp = (dim_x==1) ? deepcopy(t->p(t)[1]) : deepcopy(t->p(t))
    var = (dim_v==1) ? v[1] : v

    # misc infos
    infos = Dict()
    infos[:constraints_violation] = constraints_violation

    # +++ use proper fields instead of info
    # nonlinear constraints and multipliers
    set_constraints_and_multipliers!(infos, T, constraints_types, constraints_mult)
    # box constraints multipliers
    set_box_multipliers!(infos, T, box_multipliers, dim_x, dim_u)

    # build and return solution
    return OptimalControlSolution(ocp;
    state=fx, control=fu, objective=objective, costate=fp, times=T, variable=var, iterations=iterations, stopping=stopping, message=message, success=success, infos=infos)

end

"""
$(TYPEDSIGNATURES)
    
Process data related to constraints for solution building
"""
function set_constraints_and_multipliers!(infos, T, constraints_types, constraints_mult)

    # pure state contraints
    set_constraint_block!(infos, T, (constraints_types[1], constraints_mult[1]), (:state_constraints, :mult_state_constraints))
    # pure control constraints
    set_constraint_block!(infos, T, (constraints_types[2], constraints_mult[2]), (:control_constraints, :mult_control_constraints))
    # mixed constraints
    set_constraint_block!(infos, T, (constraints_types[3], constraints_mult[3]), (:mixed_constraints, :mult_mixed_constraints))
    # variable constraints
    set_variables_block!(infos, (constraints_types[4], constraints_mult[4]), (:variable_constraints, :mult_variable_constraints))

    return infos
end

"""
$(TYPEDSIGNATURES)
    
Process data related to a constraint type for solution building
"""
function set_constraint_block!(infos, T, const_mult, keys)
    constraint, multiplier = const_mult
    key_const, key_mult = keys
    if !isnothing(constraint)
        c = ctinterpolate(T, matrix2vec(constraint, 1))
        m = ctinterpolate(T, matrix2vec(multiplier, 1))    
        infos[key_const] = t -> c(t)
        infos[key_mult] = t -> m(t)
    end
return infos
end

"""
$(TYPEDSIGNATURES)
    
Process data related to box constraints for solution building
"""
function set_box_multipliers!(infos, T, box_multipliers, dim_x, dim_u)

    # state box
    set_box_block!(infos, T, box_multipliers[1], (:mult_state_box_lower, :mult_state_box_upper), dim_x)
    # control box
    set_box_block!(infos, T, box_multipliers[2], (:mult_control_box_lower, :mult_control_box_upper), dim_u)
    # variable box
    set_variables_block!(infos, box_multipliers[3], (:mult_variable_box_lower, :mult_variable_box_upper))

    return infos
end

"""
$(TYPEDSIGNATURES)
    
Process data related to a box type for solution building
"""
function set_box_block!(infos, T, mult, keys, dim)
    mult_l, mult_u = mults
    key_l, key_u = keys
    if !isnothing(mult_l) && !isnothing(mult_u) && dim > 0
        m_l = ctinterpolate(T, matrix2vec(mult_l[:,1:dim], 1))
        m_u = ctinterpolate(T, matrix2vec(mult_u[:,1:dim], 1))
        infos[key_l] = t -> m_l(t)
        infos[key_u] = t -> m_u(t)    
    end
    return infos
end

"""
$(TYPEDSIGNATURES)
    
Process data related to variables for solution building
"""
function set_variables_block!(infos, vecs, keys)
    vec1, vec2 = vecs
    key1, key2 = keys
    if !isnothing(vec1) && !isnothing(vec2)
        infos[key1] = vec1
        infos[key2] = vec2 
    end
    return infos
end

#= 
"""
$(TYPEDSIGNATURES)

Parse DOCP solution into OCP variables, constraints and multipliers
"""
function parse_DOCP_solution(docp, solution, multipliers_constraints, multipliers_LB, multipliers_UB, constraints)
    
    ocp = docp.ocp
    # states and controls variables, with box multipliers
    N = docp.dim_NLP_steps
    X = zeros(N+1,docp.dim_NLP_x)
    U = zeros(N+1,docp.dim_NLP_u)
    v = get_variable(solution, docp)
    # if box multipliers are empty, use dummy vectors for size consistency
    if !isnothing(multipliers_LB) && length(multipliers_LB) > 0
        mult_L = multipliers_LB
    else
        mult_L = zeros(docp.dim_NLP_variables)
    end
    if !isnothing(multipliers_UB) && length(multipliers_UB) > 0
        mult_U = multipliers_UB
    else
        mult_U = zeros(docp.dim_NLP_variables)
    end

    mult_state_box_lower = zeros(N+1,docp.dim_NLP_x)
    mult_state_box_upper = zeros(N+1,docp.dim_NLP_x)
    mult_control_box_lower = zeros(N+1,docp.dim_NLP_u)
    mult_control_box_upper = zeros(N+1,docp.dim_NLP_u)
    mult_variable_box_lower = zeros(N+1,docp.dim_NLP_v)
    mult_variable_box_upper = zeros(N+1,docp.dim_NLP_v)

    for i in 1:N+1
        # state and control variables
        X[i,:] = vget_state_at_time_step(solution, docp, i-1)
        U[i,:] = vget_control_at_time_step(solution, docp, i-1)
        # box multipliers (same layout as variables !)
        # NB. will return 0 if box constraints are not present
        mult_state_box_lower[i,:] = vget_state_at_time_step(mult_L, docp, i-1)
        mult_control_box_lower[i,:] = vget_control_at_time_step(mult_L, docp, i-1)
        mult_state_box_upper[i,:] = vget_state_at_time_step(mult_U, docp, i-1)
        mult_control_box_upper[i,:] = vget_control_at_time_step(mult_U, docp, i-1)
    end
    if dim_variable_range(ocp) > 0
        mult_variable_box_lower = get_variable(mult_L, docp)
        mult_variable_box_upper = get_variable(mult_U, docp)
    end

    # constraints, costate and constraints multipliers
    P = zeros(N, docp.dim_NLP_x)
    lambda = isnothing(multipliers_constraints) ? zeros(docp.dim_NLP_constraints) : multipliers_constraints
    sol_control_constraints = zeros(N+1,dim_control_constraints(ocp))
    sol_state_constraints = zeros(N+1,dim_state_constraints(ocp))
    sol_mixed_constraints = zeros(N+1,dim_mixed_constraints(ocp)) 
    sol_variable_constraints = zeros(dim_variable_constraints(ocp))
    mult_control_constraints = zeros(N+1,dim_control_constraints(ocp))
    mult_state_constraints = zeros(N+1,dim_state_constraints(ocp))
    mult_mixed_constraints = zeros(N+1,dim_mixed_constraints(ocp))
    mult_variable_constraints = zeros(dim_variable_constraints(ocp))

    index = 1
    for i in 1:N
        # state equation
        P[i,:] = lambda[index:index+docp.dim_NLP_x-1]
        index = index + docp.dim_NLP_x
        # path constraints
        # +++ use aux function for the 3 blocks, see eval c also
        if dim_control_constraints(ocp) > 0
            sol_control_constraints[i,:] = constraints[index:index+dim_control_constraints(ocp)-1]
            mult_control_constraints[i,:] = lambda[index:index+dim_control_constraints(ocp)-1]
            index = index + dim_control_constraints(ocp)
        end
        if dim_state_constraints(ocp) > 0
            sol_state_constraints[i,:] = constraints[index:index+dim_state_constraints(ocp)-1]
            mult_state_constraints[i,:] = lambda[index:index+dim_state_constraints(ocp)-1]
            index = index + dim_state_constraints(ocp)
        end
        if dim_mixed_constraints(ocp) > 0
            sol_mixed_constraints[i,:] = constraints[index:index+dim_mixed_constraints(ocp)-1]
            mult_mixed_constraints[i,:] = lambda[index:index+dim_mixed_constraints(ocp)-1]
            index = index + dim_mixed_constraints(ocp)
        end
    end
    # path constraints at final time
    # +++ use aux function for the 3 blocks, see eval c also
    if dim_control_constraints(ocp) > 0
        sol_control_constraints[N+1,:] = constraints[index:index+dim_control_constraints(ocp)-1]
        mult_control_constraints[N+1,:] = lambda[index:index+dim_control_constraints(ocp)-1]
        index = index + dim_control_constraints(ocp)
    end
    if dim_state_constraints(ocp) > 0
        sol_state_constraints[N+1,:] = constraints[index:index+dim_state_constraints(ocp)-1] 
        mult_state_constraints[N+1,:] = lambda[index:index+dim_state_constraints(ocp)-1]
        index = index + dim_state_constraints(ocp)
    end
    if dim_mixed_constraints(ocp) > 0
        sol_mixed_constraints[N+1,:] = constraints[index:index+dim_mixed_constraints(ocp)-1]        
        mult_mixed_constraints[N+1,:] =  lambda[index:index+dim_mixed_constraints(ocp)-1]
        index = index + dim_mixed_constraints(ocp)
    end

    # boundary conditions and multipliers
    if dim_boundary_constraints(ocp) > 0
        sol_boundary_constraints = constraints[index:index+dim_boundary_constraints(ocp)-1]
        mult_boundary_constraints = lambda[index:index+dim_boundary_constraints(ocp)-1]
        index = index + dim_boundary_constraints(ocp)
    end

    # variable constraints and multipliers
    if dim_variable_constraints(ocp) > 0
        sol_variable_constraints = constraints[index:index+dim_variable_constraints(ocp)-1]
        mult_variable_constraints = lambda[index:index+dim_variable_constraints(ocp)-1]
        index = index + dim_variable_constraints(ocp)
    end

    return X, U, v, P, sol_control_constraints, sol_state_constraints, sol_mixed_constraints, sol_variable_constraints, mult_control_constraints, mult_state_constraints, mult_mixed_constraints, mult_variable_constraints, mult_state_box_lower, mult_state_box_upper, mult_control_box_lower, mult_control_box_upper, mult_variable_box_lower, mult_variable_box_upper
end
=#

#= OLD

    # recompute value of constraints at solution
    # NB. the constraint formulation is LB <= C <= UB
    constraints = zeros(docp.dim_NLP_constraints)
    DOCP_constraints!(constraints, solution, docp)
    # set constraint violation if needed
    # +++ is not saved in OCP solution currently...
    if constraints_violation==nothing
        constraints_check = zeros(docp.dim_NLP_constraints)
        DOCP_constraints_check!(constraints_check, constraints, docp)
        println("Recomputed constraints violation ", norm(constraints_check, Inf))
        variables_check = zeros(docp.dim_NLP_variables)
        DOCP_variables_check!(variables_check, solution, docp)
        println("Recomputed variable bounds violation ", norm(variables_check, Inf))
        constraints_violation = norm(append!(variables_check, constraints_check), Inf)

    end
=#

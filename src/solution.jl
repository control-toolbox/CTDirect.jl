# Build functional OCP solution from discrete DOCP solution

"""
$(TYPEDSIGNATURES)
   
Build OCP functional solution from DOCP discrete solution (given as a SolverCore.GenericExecutionStats)
"""
function CTBase.OptimalControlSolution(docp, docp_solution)

    # retrieve data (could pass some status info too (get_status ?))
    if is_min(docp.ocp)
        objective = docp_solution.objective
    else        
        objective = - docp_solution.objective
    end

    # call lower level constructor
    return OptimalControlSolution(docp, primal=docp_solution.solution, dual=docp_solution.multipliers, objective=objective, iterations=docp_solution.iter, constraints_violation=docp_solution.primal_feas, message=String(docp_solution.solver_specific[:internal_msg]), mult_LB=docp_solution.multipliers_L, mult_UB=docp_solution.multipliers_U)

end


"""
$(TYPEDSIGNATURES)

Build OCP functional solution from the DOCP discrete solution, given as a vector. Costate will be retrieved from dual variables (multipliers) if available.
"""
function CTBase.OptimalControlSolution(docp; primal=Vector(), dual=nothing, objective=nothing, iterations=nothing, constraints_violation=nothing, message=nothing, mult_LB=nothing, mult_UB=nothing)

    # time grid
    N = docp.dim_NLP_steps
    T = zeros(N+1)
    for i=1:N+1
        T[i] = get_unnormalized_time(primal, docp, docp.NLP_normalized_time_grid[i])
    end

    # recover primal variables
    X, U, v, box_multipliers = parse_DOCP_solution_primal(docp, primal, mult_LB=mult_LB, mult_UB=mult_UB)

    # recompute / check objective
    objective_r = DOCP_objective(primal, docp)
    if isnothing(objective)
        objective = objective_r
    elseif abs((objective-objective_r)/objective) > 1e-2
        println("WARNING: recomputed objective mismatch ", objective, objective_r)
    end

    # recompute constraints
    constraints = zeros(docp.dim_NLP_constraints)
    DOCP_constraints!(constraints, primal, docp)

    # recover costate and constraints with their multipliers
    P, constraints_types, constraints_mult = parse_DOCP_solution_dual(docp, dual, constraints)

    # call lowest level constructor
    return OptimalControlSolution(docp, T, X, U, v, P, objective=objective, iterations=iterations, constraints_violation=constraints_violation, message=message,
    constraints_types=constraints_types, constraints_mult=constraints_mult,
    box_multipliers=box_multipliers)
    
end


"""
$(TYPEDSIGNATURES)

Recover OCP primal variables from DOCP solution
"""
function parse_DOCP_solution_primal(docp, solution; mult_LB=nothing, mult_UB=nothing)

    # state and control variables
    N = docp.dim_NLP_steps
    X = zeros(N+1,docp.dim_NLP_x)
    U = zeros(N+1,docp.dim_NLP_u)

    # multipliers for box constraints
    if isnothing(mult_LB) || length(mult_LB) == 0
        mult_LB = zeros(docp.dim_NLP_variables)
    end
    if isnothing(mult_UB) || length(mult_UB) == 0
        mult_UB = zeros(docp.dim_NLP_variables)
    end
    mult_state_box_lower = zeros(N+1,docp.dim_NLP_x)
    mult_state_box_upper = zeros(N+1,docp.dim_NLP_x)
    mult_control_box_lower = zeros(N+1,docp.dim_NLP_u)
    mult_control_box_upper = zeros(N+1,docp.dim_NLP_u)
    mult_variable_box_lower = zeros(N+1,docp.dim_NLP_v)
    mult_variable_box_upper = zeros(N+1,docp.dim_NLP_v)

    # retrieve optimization variables
    v = get_variable(solution, docp)
    mult_variable_box_lower = get_variable(mult_LB, docp)
    mult_variable_box_upper = get_variable(mult_UB, docp)

    # loop over time steps
    for i in 1:N+1
        # state and control variables at current step
        X[i,:], U[i,:] = get_NLP_variables_at_time_step(solution, docp, i-1)

        # box multipliers
        mult_state_box_lower[i,:], mult_control_box_lower[i,:] = get_NLP_variables_at_time_step(mult_LB, docp, i-1)
        mult_state_box_upper[i,:], mult_control_box_upper[i,:] = get_NLP_variables_at_time_step(mult_UB, docp, i-1)

    end

    box_multipliers=((mult_state_box_lower,mult_state_box_upper), (mult_control_box_lower,mult_control_box_upper), (mult_variable_box_lower,mult_variable_box_upper))

    return X, U, v, box_multipliers
end


"""
$(TYPEDSIGNATURES)

Recover OCP costate and constraints multipliers from DOCP multipliers
"""
function parse_DOCP_solution_dual(docp, multipliers, constraints)

    # constraints tuple: (state, control, mixed, variable, boundary)
    N = docp.dim_NLP_steps
    P = zeros(N, docp.dim_NLP_x)

    ocp = docp.ocp
    dcc = dim_control_constraints(ocp)
    dsc = dim_state_constraints(ocp)
    dmc = dim_mixed_constraints(ocp)
    dvc = dim_variable_constraints(ocp)
    dbc = dim_boundary_constraints(ocp)

    sol_state_constraints = zeros(N+1,dsc)
    sol_control_constraints = zeros(N+1,dcc)
    sol_mixed_constraints = zeros(N+1,dmc) 
    sol_variable_constraints = zeros(dvc)
    sol_boundary_constraints = zeros(dbc)

    mult_state_constraints = zeros(N+1,dsc)
    mult_control_constraints = zeros(N+1,dcc)
    mult_mixed_constraints = zeros(N+1,dmc)
    mult_variable_constraints = zeros(dvc)
    mult_boundary_constraints = zeros(dbc)

    # if called with multipliers = nothing, fill with zeros
    if isnothing(multipliers) 
        multipliers = zeros(docp.dim_NLP_constraints)
    end

    # loop over time steps
    i_c = 1
    i_m = 1
    for i in 1:N+1

        # state equation multiplier for costate (except last step)
        if i < N+1 
            P[i,:] = multipliers[i_m : i_m+docp.dim_NLP_x-1]
            i_m += docp.dim_NLP_x
        end

        # path constraints and multipliers
        # pure control constraints
        if dcc > 0
            sol_control_constraints[i,:] = constraints[i_c : i_c+dcc-1]
            mult_control_constraints[i,:] = multipliers[i_m : i_m+dcc-1]
            i_c += dcc
            i_m += dcc
        end
        # pure state constraints
        if dsc > 0
            sol_state_constraints[i,:] = constraints[i_c : i_c+dsc-1]
            mult_state_constraints[i,:] = multipliers[i_m : i_m+dsc-1]
            i_c += dsc
            i_m += dsc
        end
        # mixed constraints
        if dmc > 0
            sol_mixed_constraints[i,:] = constraints[i_c : i_c+dmc-1]
            mult_mixed_constraints[i,:] = multipliers[i_m : i_m+dmc-1]
            i_c += dmc
            i_m += dmc
        end

    end

    # pointwise constraints: boundary then variables
    if dbc > 0
        sol_boundary_constraints[:] = constraints[i_c : i_c+dbc-1]
        mult_boundary_constraints[:] = multipliers[i_m : i_m+dbc-1]
        i_c += dbc
        i_m += dbc
    end
    if dvc > 0
        sol_variable_constraints[:] = constraints[i_c : i_c+dvc-1]
        mult_variable_constraints[:] = multipliers[i_m : i_m+dvc-1]
        i_c += dvc
        i_m += dvc
    end    

    # return tuples
    constraints_types = (sol_state_constraints, sol_control_constraints,  sol_mixed_constraints, sol_variable_constraints, sol_boundary_constraints)
    constraints_mult = (mult_state_constraints, mult_control_constraints, mult_mixed_constraints, mult_variable_constraints, mult_boundary_constraints)

    return P, constraints_types, constraints_mult
end


"""
$(TYPEDSIGNATURES)
    
Build OCP functional solution from DOCP vector solution (given as raw variables and multipliers plus some optional infos)
"""
# try to reuse this for the discrete json solution ?
function CTBase.OptimalControlSolution(docp, T, X, U, v, P;
    objective=0, iterations=0, constraints_violation=0, message="No msg", stopping=nothing, success=nothing,
    constraints_types=(nothing,nothing,nothing,nothing,nothing),
    constraints_mult=(nothing,nothing,nothing,nothing,nothing),
    box_multipliers=((nothing,nothing),(nothing,nothing),(nothing,nothing)))

    ocp = docp.ocp
    dim_x = state_dimension(ocp)
    dim_u = control_dimension(ocp)
    dim_v = variable_dimension(ocp)

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
    infos = Dict{Symbol,Any}()
    infos[:constraints_violation] = constraints_violation

    # NB later use proper fields instead of info
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

    # keys list: state, control, mixed, variable, boundary
    key_list = ((:state_constraints, :mult_state_constraints),
                (:control_constraints, :mult_control_constraints),
                (:mixed_constraints, :mult_mixed_constraints),
                (:variable_constraints, :mult_variable_constraints),
                (:boundary_constraints, :mult_boundary_constraints))

    for i=1:3
        set_constraint_block!(infos, T, (constraints_types[i], constraints_mult[i]), key_list[i])
    end
    for i=4:5
        set_variables_block!(infos, (constraints_types[i], constraints_mult[i]), key_list[i])
    end

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
function set_box_block!(infos, T, mults, keys, dim)
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

#= OLD

    # recompute value of constraints at solution
    # NB. the constraint formulation is LB <= C <= UB
    constraints = zeros(docp.dim_NLP_constraints)
    DOCP_constraints!(constraints, solution, docp)
    # set constraint violation if needed
    # is not saved in OCP solution currently...
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

# Build functional OCP solution from discrete DOCP solution

"""
$(TYPEDSIGNATURES)
   
Build OCP functional solution from DOCP discrete solution (given as a SolverCore.GenericExecutionStats)
"""
function CTBase.OptimalControlSolution(docp::DOCP, docp_solution)

    # retrieve data (could pass some status info too (get_status ?))
    if docp.is_maximization
        objective = -docp_solution.objective
    else
        objective = docp_solution.objective
    end

    # call lower level constructor
    return OptimalControlSolution(
        docp,
        primal = docp_solution.solution,
        dual = docp_solution.multipliers,
        objective = objective,
        iterations = docp_solution.iter,
        constraints_violation = docp_solution.primal_feas,
        message = string(docp_solution.solver_specific[:internal_msg][1]),  #+++ use GENERIC version and test for specific fields instead ! madnlp version can then probably be removed. 
        mult_LB = docp_solution.multipliers_L,
        mult_UB = docp_solution.multipliers_U,
    )
end

"""
$(TYPEDSIGNATURES)

Build OCP functional solution from the DOCP discrete solution, given as a vector. Costate will be retrieved from dual variables (multipliers) if available.
"""
function CTBase.OptimalControlSolution(
    docp::DOCP;
    primal = Vector(),
    dual = nothing,
    objective = nothing,
    iterations = nothing,
    constraints_violation = nothing,
    message = nothing,
    mult_LB = nothing,
    mult_UB = nothing,
)

    # time grid
    T = get_time_grid(primal, docp)

    # recover primal variables
    X, U, v, box_multipliers =
        parse_DOCP_solution_primal(docp, primal, mult_LB = mult_LB, mult_UB = mult_UB)

    # recompute / check objective
    objective_r = DOCP_objective(primal, docp)
    if docp.is_maximization
        objective_r = -objective_r
    end
    if isnothing(objective)
        objective = objective_r
    elseif abs((objective - objective_r) / objective) > 1e-2
        println("WARNING: recomputed objective mismatch ", objective, objective_r)
    end

    # recompute constraints
    constraints = zeros(docp.dim_NLP_constraints)
    DOCP_constraints!(constraints, primal, docp)

    # recover costate and constraints with their multipliers
    P, constraints_types, constraints_mult = parse_DOCP_solution_dual(docp, dual, constraints)

    # call lowest level constructor
    return OptimalControlSolution(
        docp.ocp,
        T,
        X,
        U,
        v,
        P,
        objective = objective,
        iterations = iterations,
        constraints_violation = constraints_violation,
        message = message,
        constraints_types = constraints_types,
        constraints_mult = constraints_mult,
        box_multipliers = box_multipliers
    )
end

"""
$(TYPEDSIGNATURES)

Recover OCP primal variables from DOCP solution
"""
function parse_DOCP_solution_primal(docp, solution; mult_LB = nothing, mult_UB = nothing)

    # state and control variables
    N = docp.dim_NLP_steps
    X = zeros(N + 1, docp.dim_OCP_x)
    U = zeros(N + 1, docp.dim_NLP_u)
    v = Float64[]

    # multipliers for box constraints
    if isnothing(mult_LB) || length(mult_LB) == 0
        mult_LB = zeros(docp.dim_NLP_variables)
    end
    if isnothing(mult_UB) || length(mult_UB) == 0
        mult_UB = zeros(docp.dim_NLP_variables)
    end
    mult_state_box_lower = zeros(size(X))
    mult_state_box_upper = zeros(size(X))
    mult_control_box_lower = zeros(size(U))
    mult_control_box_upper = zeros(size(U))
    mult_variable_box_lower = zeros(size(v))
    mult_variable_box_upper = zeros(size(v))

    # retrieve optimization variables
    if docp.is_variable
        v = get_OCP_variable(solution, docp)
        mult_variable_box_lower = get_OCP_variable(mult_LB, docp)
        mult_variable_box_upper = get_OCP_variable(mult_UB, docp)
    end

    # state variables and box multipliers
    for i = 1:N+1
        X[i,:] .= get_OCP_state_at_time_step(solution, docp, i)
        mult_state_box_lower[i, :] .= get_OCP_state_at_time_step(mult_LB, docp, i)
        mult_state_box_upper[i, :] .= get_OCP_state_at_time_step(mult_UB, docp, i)
    end
    # control variables and box multipliers
    for i = 1:N+1
        U[i,:] .= get_OCP_control_at_time_step(solution, docp, i)
        mult_control_box_lower[i, :] .= get_OCP_control_at_time_step(mult_LB, docp, i)
        mult_control_box_upper[i, :] .= get_OCP_control_at_time_step(mult_UB, docp, i)
    end

    box_multipliers = (
        mult_state_box_lower, mult_state_box_upper,
        mult_control_box_lower, mult_control_box_upper,
        mult_variable_box_lower, mult_variable_box_upper
    )

    return X, U, v, box_multipliers
end

"""
$(TYPEDSIGNATURES)

Recover OCP costate and constraints multipliers from DOCP multipliers
"""
function parse_DOCP_solution_dual(docp, multipliers, constraints)

    # if called with multipliers = nothing, fill with zeros
    if isnothing(multipliers)
        multipliers = zeros(docp.dim_NLP_constraints)
    end

    # constraints tuple: (state, control, mixed, variable, boundary)
    N = docp.dim_NLP_steps
    P = zeros(N, docp.dim_NLP_x)

    ocp = docp.ocp
    dcc = dim_control_constraints(ocp)
    dsc = dim_state_constraints(ocp)
    dmc = dim_mixed_constraints(ocp)
    dbc = dim_boundary_constraints(ocp)
    dvc = dim_variable_constraints(ocp)

    # constraints
    sol_state_constraints = zeros(N + 1, dsc)
    sol_boundary_constraints = zeros(dbc)
    sol_variable_constraints = zeros(dvc)
    sol_control_constraints = zeros(N + 1, dcc)
    sol_mixed_constraints = zeros(N + 1, dmc)

    # constraints multipliers
    mul_control_constraints = zeros(size(sol_control_constraints))
    mul_state_constraints = zeros(size(sol_state_constraints))
    mul_mixed_constraints = zeros(size(sol_mixed_constraints))
    mul_boundary_constraints = zeros(size(sol_boundary_constraints))
    mul_variable_constraints = zeros(size(sol_variable_constraints))

    # loop over time steps
    i_c = 1
    i_m = 1
    for i = 1:N+1

        # state equation multiplier for costate
        if i <= N
            P[i, :] = multipliers[i_m:(i_m + docp.dim_NLP_x - 1)]
            # skip state / stage constraints
            i_c += docp.discretization._state_stage_eqs_block
            i_m += docp.discretization._state_stage_eqs_block
        end

        # path constraints and multipliers
        # pure control constraints
        if dcc > 0
            sol_control_constraints[i, :] = constraints[i_c:(i_c + dcc - 1)]
            mul_control_constraints[i, :] = multipliers[i_m:(i_m + dcc - 1)]
            i_c += dcc
            i_m += dcc
        end
        # pure state constraints
        if dsc > 0
            sol_state_constraints[i, :] = constraints[i_c:(i_c + dsc - 1)]
            mul_state_constraints[i, :] = multipliers[i_m:(i_m + dsc - 1)]
            i_c += dsc
            i_m += dsc
        end
        # mixed constraints
        if dmc > 0
            sol_mixed_constraints[i, :] = constraints[i_c:(i_c + dmc - 1)]
            mul_mixed_constraints[i, :] = multipliers[i_m:(i_m + dmc - 1)]
            i_c += dmc
            i_m += dmc
        end
    end

    # pointwise constraints: boundary then variables
    if dbc > 0
        sol_boundary_constraints[:] = constraints[i_c:(i_c + dbc - 1)]
        mul_boundary_constraints[:] = multipliers[i_m:(i_m + dbc - 1)]
        i_c += dbc
        i_m += dbc
    end
    if dvc > 0
        sol_variable_constraints[:] = constraints[i_c:(i_c + dvc - 1)]
        mul_variable_constraints[:] = multipliers[i_m:(i_m + dvc - 1)]
        i_c += dvc
        i_m += dvc
    end

    # return tuples
    constraints_types = (
        sol_control_constraints,
        sol_state_constraints,
        sol_mixed_constraints,
        sol_boundary_constraints,
        sol_variable_constraints,
    )
    constraints_mult = (
        mul_control_constraints,
        mul_state_constraints,
        mul_mixed_constraints,
        mul_boundary_constraints,
        mul_variable_constraints,
    )

    return P, constraints_types, constraints_mult
end

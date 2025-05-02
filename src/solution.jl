# Build functional OCP solution from discrete DOCP solution

"""
$(TYPEDSIGNATURES)
   
Build OCP functional solution from DOCP discrete solution (given as a SolverCore.GenericExecutionStats)
"""
function build_OCP_solution(docp, docp_solution)

    ocp = docp.ocp
    solution = docp_solution.solution
    iterations, constraints_violation, message, stopping, success = SolverInfos(docp_solution)

    # time grid
    T = get_time_grid(solution, docp)

    # primal variables X, U, v and box multipliers
    X, U, v, box_multipliers = parse_DOCP_solution_primal(docp, solution; mult_LB=docp_solution.multipliers_L, mult_UB=docp_solution.multipliers_U)

    # recompute / check objective
    objective_r = DOCP_objective(solution, docp)
    if docp.flags.max
        objective = -docp_solution.objective
        objective_r = -objective_r
    else
        objective = docp_solution.objective
    end
    if isnothing(objective)
        objective = objective_r
    elseif abs((objective - objective_r) / objective) > 1e-2
        println("WARNING: recomputed objective mismatch ", objective, objective_r)
    end

    # recompute constraints
    constraints = zeros(docp.dim_NLP_constraints)
    DOCP_constraints!(constraints, solution, docp)

    # costate and constraints multipliers
    P, path_constraints, boundary_constraints, path_constraints_dual, boundary_constraints_dual = parse_DOCP_solution_dual(docp, docp_solution.multipliers, constraints)

    return CTModels.build_solution(
        ocp,
        T, X, U, v, P;
        objective=objective, iterations=iterations, constraints_violation=constraints_violation,
        message=message, stopping=stopping, success=success,
        path_constraints=path_constraints,
        path_constraints_dual=path_constraints_dual,
        boundary_constraints=boundary_constraints,
        boundary_constraints_dual=boundary_constraints_dual,
        state_constraints_lb_dual=box_multipliers[1],
        state_constraints_ub_dual=box_multipliers[2],
        control_constraints_lb_dual=box_multipliers[3],
        control_constraints_ub_dual=box_multipliers[4],
        variable_constraints_lb_dual=box_multipliers[5],
        variable_constraints_ub_dual=box_multipliers[6]
    )
end


"""
$(TYPEDSIGNATURES)

Retrieve convergence information (Ipopt version)
"""
function SolverInfos()
    return 0, 0., "undefined", :undefined, true
end
function SolverInfos(docp_solution)

    # try to detect solver here for specific fields !
    iterations = docp_solution.iter
    constraints_violation = docp_solution.primal_feas
    stopping = :undefined
    success = true
    message = "undefined"
    try
        solver_specific = docp_solution.solver_specific
        if haskey(solver_specific, :internal_msg)
            # Ipopt solver
            message = string(solver_specific[:internal_msg][1])
        end
    catch e # missing field solve_specific
    end

    return iterations, constraints_violation, message, stopping, success
end


"""
$(TYPEDSIGNATURES)

Build OCP functional solution from DOCP discrete solution (given as a SolverCore.GenericExecutionStats)
"""
function build_OCP_solution(docp; primal, dual=nothing, mult_LB=nothing, mult_UB=nothing)

    ocp = docp.ocp
    solution = primal
    iterations, constraints_violation, message, stopping, success = SolverInfos()

    # time grid
    T = get_time_grid(solution, docp)

    # primal variables X, U, v and box multipliers
    X, U, v, box_multipliers = parse_DOCP_solution_primal(docp, solution; mult_LB=mult_LB, mult_UB=mult_UB)

    # recompute objective
    objective = DOCP_objective(solution, docp)
    if docp.flags.max
        objective = -objective
    end

    # recompute constraints
    constraints = zeros(docp.dim_NLP_constraints)
    DOCP_constraints!(constraints, solution, docp)

    # costate and constraints multipliers
    P, path_constraints, boundary_constraints, path_constraints_dual, boundary_constraints_dual = parse_DOCP_solution_dual(docp, dual, constraints)

    return CTModels.build_solution(
        ocp,
        T, X, U, v, P;
        objective=objective, iterations=iterations, constraints_violation=constraints_violation,
        message=message, stopping=stopping, success=success,
        path_constraints=path_constraints,
        path_constraints_dual=path_constraints_dual,
        boundary_constraints=boundary_constraints,
        boundary_constraints_dual=boundary_constraints_dual,
        state_constraints_lb_dual=box_multipliers[1],
        state_constraints_ub_dual=box_multipliers[2],
        control_constraints_lb_dual=box_multipliers[3],
        control_constraints_ub_dual=box_multipliers[4],
        variable_constraints_lb_dual=box_multipliers[5],
        variable_constraints_ub_dual=box_multipliers[6]
    )
end


"""
$(TYPEDSIGNATURES)

Recover OCP primal variables from DOCP solution
"""
function parse_DOCP_solution_primal(docp, solution; mult_LB=nothing, mult_UB=nothing)

    # state and control variables
    N = docp.time.steps
    X = zeros(N + 1, docp.dims.OCP_x)
    U = zeros(N + 1, docp.dims.NLP_u)
    v = zeros(docp.dims.NLP_v)

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
    if docp.dims.NLP_v > 0
        v .= get_OCP_variable(solution, docp)
        mult_variable_box_lower .= get_OCP_variable(mult_LB, docp)
        mult_variable_box_upper .= get_OCP_variable(mult_UB, docp)
    end

    # state variables and box multipliers
    for i = 1:(N+1)
        X[i, :] .= get_OCP_state_at_time_step(solution, docp, i)
        mult_state_box_lower[i, :] .= get_OCP_state_at_time_step(mult_LB, docp, i)
        mult_state_box_upper[i, :] .= get_OCP_state_at_time_step(mult_UB, docp, i)
    end
    # control variables and box multipliers
    for i = 1:(N+1)
        U[i, :] .= get_OCP_control_at_time_step(solution, docp, i)
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

    # costate
    N = docp.time.steps
    P = zeros(N, docp.dims.NLP_x)
    ocp = docp.ocp

    # constraints
    dpc = docp.dims.path_cons
    dbc = docp.dims.boundary_cons
    sol_path_constraints = zeros(N + 1, dpc)
    sol_boundary_constraints = zeros(dbc)

    # constraints multipliers
    mul_path_constraints = zeros(size(sol_path_constraints))
    mul_boundary_constraints = zeros(size(sol_boundary_constraints))

    # loop over time steps
    i_c = 1
    i_m = 1
    for i = 1:(N+1)

        # state equation multiplier for costate
        if i <= N
            P[i, :] = multipliers[i_m:(i_m+docp.dims.NLP_x-1)]
            # skip state / stage constraints
            i_c += docp.discretization._state_stage_eqs_block
            i_m += docp.discretization._state_stage_eqs_block
        end

        # path constraints and multipliers
        if dpc > 0
            sol_path_constraints[i, :] = constraints[i_c:(i_c+dpc-1)]
            mul_path_constraints[i, :] = multipliers[i_m:(i_m+dpc-1)]
            i_c += dpc
            i_m += dpc
        end
    end

    # pointwise constraints: boundary then variables
    if dbc > 0
        sol_boundary_constraints[:] = constraints[i_c:(i_c+dbc-1)]
        mul_boundary_constraints[:] = multipliers[i_m:(i_m+dbc-1)]
        i_c += dbc
        i_m += dbc
    end

    return P, sol_path_constraints, sol_boundary_constraints, mul_path_constraints, mul_boundary_constraints
end

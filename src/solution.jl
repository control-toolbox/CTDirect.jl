# Build functional OCP solution from discrete DOCP solution

is_empty(t) = (isnothing(t) || length(t) == 0)

"""
$(TYPEDSIGNATURES)
   
Build OCP functional solution from DOCP discrete solution 
(given as a SolverCore.GenericExecutionStats)
"""
function build_OCP_solution(docp, nlp_solution; nlp_model=ADNLPBackend(), nlp_solver=IpoptBackend())

    ocp = ocp_model(docp)

    # retrieve data from NLP solver
    objective, iterations, constraints_violation, message, status, successful = SolverInfos(nlp_solution)
    # fix objective sign for maximization problems with MadNLP
    # should be in Solverinfos but needs max info. can we retrieve it from nlp solution ?
    if docp.flags.max && nlp_solver isa MadNLPBackend
        objective = - objective
    end

    # arrays (explicit conversion for GPU case)
    solution = Array(nlp_solution.solution)
    multipliers = Array(nlp_solution.multipliers)
    multipliers_L = Array(nlp_solution.multipliers_L)    
    multipliers_U = Array(nlp_solution.multipliers_U)

    # time grid
    if nlp_model isa ADNLPBackend
        T = get_time_grid(solution, docp)
    else
        T = get_time_grid_exa(nlp_solution, docp)
    end

    # +++ todo: replace both parsing functions with series of getter calls
    # unify adnlp / exa cases

    # primal variables X, U, v and box multipliers
    X, U, v, box_multipliers = parse_DOCP_solution_primal(docp, solution; multipliers_L=multipliers_L, multipliers_U=multipliers_U, nlp_model=nlp_model, nlp_solution=nlp_solution)

    # costate and constraints multipliers
    P, path_constraints_dual, boundary_constraints_dual = parse_DOCP_solution_dual(docp, multipliers; nlp_model=nlp_model, nlp_solution=nlp_solution)

    return CTModels.build_solution(
        ocp,
        T,
        X,
        U,
        v,
        P;
        objective=objective,
        iterations=iterations,
        constraints_violation=constraints_violation,
        message=message,
        status=status,
        successful=successful,
        path_constraints_dual=path_constraints_dual,
        boundary_constraints_dual=boundary_constraints_dual,
        state_constraints_lb_dual=box_multipliers[1],
        state_constraints_ub_dual=box_multipliers[2],
        control_constraints_lb_dual=box_multipliers[3],
        control_constraints_ub_dual=box_multipliers[4],
        variable_constraints_lb_dual=box_multipliers[5],
        variable_constraints_ub_dual=box_multipliers[6],
    )
end

"""
$(TYPEDSIGNATURES)

Retrieve convergence information from NLP solution (SolverCore.ExecutionStats)
- objective [Float]: objective value at the solution
- iterations [Integer]: number of iterations
- constraints_violations [Real]: primal feasibility
- message [String]: optional solver dependent message
- status [Symbol]: termination status from the NLP solver
- successful [Boolean]: indicates successful convergence (first order)
"""
function SolverInfos()
    return 0., 0, 0., "undefined", :undefined, true
end
function SolverInfos(nlp_solution)

    objective = nlp_solution.objective
    iterations = nlp_solution.iter
    constraints_violation = nlp_solution.primal_feas
    status = nlp_solution.status
    successful = (status == :first_order) || (status == :acceptable)

    return objective, iterations, constraints_violation, "Ipopt/generic", status, successful
end


"""
$(TYPEDSIGNATURES)

Build OCP functional solution from DOCP discrete solution 
(given as array for primal variables, optionally dual variables and bounds multipliers)
"""
function build_OCP_solution(
    docp;
    primal,
    dual=nothing,
    multipliers_L=nothing,
    multipliers_U=nothing,
    nlp_model=ADNLPBackend(),
    nlp_solution,
)
    ocp = ocp_model(docp)
    solution = primal

    # dummy info
    objective, iterations, constraints_violation, message, status, successful = SolverInfos()

    # recompute objective
    objective = DOCP_objective(solution, docp)

    # time grid
    if nlp_model isa ADNLPBackend
        T = get_time_grid(solution, docp)
    else
        T = get_time_grid_exa(nlp_solution, docp)
    end

    # primal variables X, U, v and box multipliers
    X, U, v, box_multipliers = parse_DOCP_solution_primal(
        docp,
        solution;
        multipliers_L=multipliers_L,
        multipliers_U=multipliers_U,
        nlp_model=nlp_model,
        nlp_solution=nlp_solution,
    )

    # costate and constraints multipliers
    P, path_constraints_dual, boundary_constraints_dual = parse_DOCP_solution_dual(
        docp, dual; nlp_model=nlp_model, nlp_solution=nlp_solution
    )

    return CTModels.build_solution(
        ocp,
        T,
        X,
        U,
        v,
        P;
        objective=objective,
        iterations=iterations,
        constraints_violation=constraints_violation,
        message=message,
        status=status,
        successful=successful,
        path_constraints_dual=path_constraints_dual,
        boundary_constraints_dual=boundary_constraints_dual,
        state_constraints_lb_dual=box_multipliers[1],
        state_constraints_ub_dual=box_multipliers[2],
        control_constraints_lb_dual=box_multipliers[3],
        control_constraints_ub_dual=box_multipliers[4],
        variable_constraints_lb_dual=box_multipliers[5],
        variable_constraints_ub_dual=box_multipliers[6],
    )
end

"""
$(TYPEDSIGNATURES)

Recover OCP state, control and optimization variables from DOCP primal variables.
Bounds multipliers will be parsed as well if present.
"""
function parse_DOCP_solution_primal(
    docp,
    solution;
    multipliers_L,
    multipliers_U,
    nlp_model,
    nlp_solution
)

    # state and control variables
    N = docp.time.steps
    X = zeros(N + 1, docp.dims.OCP_x)
    U = zeros(N + 1, docp.dims.NLP_u)
    v = zeros(docp.dims.NLP_v)

    # multipliers for box constraints
    mult_state_box_lower = zeros(size(X))
    mult_state_box_upper = zeros(size(X))
    mult_control_box_lower = zeros(size(U))
    mult_control_box_upper = zeros(size(U))
    mult_variable_box_lower = zeros(size(v))
    mult_variable_box_upper = zeros(size(v))

    if nlp_model isa ExaBackend # Exa
        getter = docp.exa_getter
        X[:] = getter(nlp_solution; val=:state)' # transpose to match choice below for ADNLP
        U[:] = getter(nlp_solution; val=:control)'
        v[:] = getter(nlp_solution; val=:variable)
        if !is_empty(multipliers_L)
            mult_state_box_lower[:] = getter(nlp_solution; val=:state_l)'
            mult_control_box_lower[:] = getter(nlp_solution; val=:control_l)'
            mult_variable_box_lower[:] = getter(nlp_solution; val=:variable_l)
        end
        if !is_empty(multipliers_U) 
            mult_state_box_upper[:] = getter(nlp_solution; val=:state_u)'
            mult_control_box_upper[:] = getter(nlp_solution; val=:control_u)'
            mult_variable_box_upper[:] = getter(nlp_solution; val=:variable_u)
        end

    else # ADNLP

        # replace ipopt 0-length arrays with full 0 arrays
        is_empty(multipliers_L) && (multipliers_L = zeros(docp.dim_NLP_variables))
        is_empty(multipliers_U) && (multipliers_U = zeros(docp.dim_NLP_variables))

        # retrieve optimization variables
        if docp.dims.NLP_v > 0
            v .= get_OCP_variable(solution, docp)
            mult_variable_box_lower .= get_OCP_variable(multipliers_L, docp)
            mult_variable_box_upper .= get_OCP_variable(multipliers_U, docp)
        end

        # state variables and box multipliers
        for i in 1:(N + 1)
            X[i, :] .= get_OCP_state_at_time_step(solution, docp, i)
            mult_state_box_lower[i, :] .= get_OCP_state_at_time_step(multipliers_L, docp, i)
            mult_state_box_upper[i, :] .= get_OCP_state_at_time_step(multipliers_U, docp, i)
        end
        # control variables and box multipliers
        for i in 1:(N + 1)
            U[i, :] .= get_OCP_control_at_time_step(solution, docp, i)
            mult_control_box_lower[i, :] .= get_OCP_control_at_time_step(multipliers_L, docp, i)
            mult_control_box_upper[i, :] .= get_OCP_control_at_time_step(multipliers_U, docp, i)
        end
    end

    box_multipliers = (
        mult_state_box_lower,
        mult_state_box_upper,
        mult_control_box_lower,
        mult_control_box_upper,
        mult_variable_box_lower,
        mult_variable_box_upper,
    )

    return X, U, v, box_multipliers
end

"""
$(TYPEDSIGNATURES)

Recover OCP costate and constraints multipliers from DOCP dual variables.
"""
function parse_DOCP_solution_dual(
    docp, multipliers; nlp_model=ADNLPBackend(), nlp_solution
)

    # costate
    N = docp.time.steps
    P = zeros(N, docp.dims.NLP_x)

    if nlp_model isa ExaBackend # Exa
        getter = docp.exa_getter
        P[:] = getter(nlp_solution; val=:costate)' # transpose to match choice below for ADNLP
        dpc = docp.dims.path_cons
        dbc = docp.dims.boundary_cons
        mul_path_constraints = zeros(N + 1, dpc) # todo: add getters for path constraints for :exa in CTParser
        mul_boundary_constraints = zeros(dbc) # todo: add getters for boundary constraints for :exa in CTParser

    else # ADNLP

        disc = disc_model(docp)

        # if called with multipliers = nothing, fill with zeros
        isnothing(multipliers) && (multipliers = zeros(docp.dim_NLP_constraints))

        # convert GPU arrays if needed
        multipliers = Array(multipliers)

        # dimensions
        dpc = docp.dims.path_cons
        dbc = docp.dims.boundary_cons

        # constraints multipliers
        mul_path_constraints = zeros(N + 1, dpc)
        mul_boundary_constraints = zeros(dbc)

        # loop over time steps
        i_m = 1
        for i in 1:(N + 1)

            # state equation multiplier for costate
            if i <= N
                P[i, :] = multipliers[i_m:(i_m + docp.dims.NLP_x - 1)]
                # skip state / stage constraints
                i_m += disc._state_stage_eqs_block
            end

            # path constraints and multipliers
            if dpc > 0
                mul_path_constraints[i, :] = multipliers[i_m:(i_m + dpc - 1)]
                i_m += dpc
            end
        end

        # pointwise constraints: boundary then variables
        if dbc > 0
            mul_boundary_constraints[:] = multipliers[i_m:(i_m + dbc - 1)]
            i_m += dbc
        end
    end

    return P, mul_path_constraints, mul_boundary_constraints
end

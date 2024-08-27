# CTDirect interface

"""
$(TYPEDSIGNATURES)

Return the list of available methods to solve the optimal control problem.
"""
function available_methods()
    # available methods by order of preference
    algorithms = ()
    algorithms = add(algorithms, (:adnlp, :ipopt))
    algorithms = add(algorithms, (:adnlp, :madnlp))
    return algorithms
end

"""
$(TYPEDSIGNATURES)

Discretize an optimal control problem into a nonlinear optimization problem (ie direct transcription)
"""
function direct_transcription(
    ocp::OptimalControlModel,
    description...;
    init = CTBase.__ocp_init(),
    grid_size = __grid_size(),
    time_grid = __time_grid(),
)

    # build DOCP
    docp = DOCP(ocp, grid_size, time_grid)

    # set bounds in DOCP
    variables_bounds!(docp)
    constraints_bounds!(docp)

    # set initial guess
    x0 = DOCP_initial_guess(
        docp,
        OptimalControlInit(
            init,
            state_dim = ocp.state_dimension,
            control_dim = ocp.control_dimension,
            variable_dim = ocp.variable_dimension,
        ),
    )

    # call NLP problem constructor
    nlp = ADNLPModel!(
        x -> DOCP_objective(x, docp),
        x0,
        docp.var_l,
        docp.var_u,
        (c, x) -> DOCP_constraints!(c, x, docp),
        docp.con_l,
        docp.con_u,
        backend = :optimized,
    )

    return docp, nlp
end

"""
$(TYPEDSIGNATURES)

Set initial guess in the DOCP
"""
function set_initial_guess(docp::DOCP, nlp, init)
    ocp = docp.ocp
    nlp.meta.x0 .= DOCP_initial_guess(
        docp,
        OptimalControlInit(
            init,
            state_dim = ocp.state_dimension,
            control_dim = ocp.control_dimension,
            variable_dim = ocp.variable_dimension,
        ),
    )
end

"""
$(TYPEDSIGNATURES)

Solve an OCP with a direct method
"""
function direct_solve(
    ocp::OptimalControlModel,
    description::Symbol...;
    init = CTBase.__ocp_init(),
    grid_size::Int = CTDirect.__grid_size(),
    time_grid = CTDirect.__time_grid(),
    kwargs...,
)
    method = getFullDescription(description, available_methods())
    #println(method)

    # build discretized OCP, including initial guess
    docp, nlp = direct_transcription(
        ocp,
        description,
        init = init,
        grid_size = grid_size,
        time_grid = time_grid,
    )

    # solve DOCP
    if :ipopt ∈ method
        tag = CTDirect.IpoptTag()
    elseif :madnlp ∈ method
        tag = CTDirect.MadNLPTag()
    else
        error("no known solver in method", method)
    end
    docp_solution = CTDirect.solve_docp(tag, docp, nlp; kwargs...)

    # build and return OCP solution
    return OptimalControlSolution(docp, docp_solution)
end

# placeholders (see CTSolveExt*** extensions)
abstract type SolverTag end
struct IpoptTag <: SolverTag end
struct MadNLPTag <: SolverTag end

function solve_docp(solver_tag, args...; kwargs...)
    if typeof(solver_tag) == IpoptTag
        throw(ExtensionError(:NLPModelsIpopt))
    elseif typeof(solver_tag) == MadNLPTag
        throw(ExtensionError(:MadNLP))
    else
        error("Unknown solver type", typeof(solver_tag))
    end
end

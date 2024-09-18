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
    disc_method = __disc_method()
)

    # build DOCP
    docp = DOCP(ocp; grid_size=grid_size, time_grid=time_grid, disc_method=string(disc_method))

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
        #x -> DOCP_objective(x, docp),
        docp.objective,
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
    disc_method = __disc_method(),
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
        disc_method = disc_method
    )

    # solve DOCP
    if :ipopt ∈ method
        solver_backend = CTDirect.IpoptBackend()
    elseif :madnlp ∈ method
        solver_backend = CTDirect.MadNLPBackend()
    else
        error("no known solver in method", method)
    end
    docp_solution = CTDirect.solve_docp(solver_backend, docp, nlp; kwargs...)

    # build and return OCP solution
    return OptimalControlSolution(docp, docp_solution)
end

# placeholders (see CTSolveExt*** extensions)
abstract type AbstractSolverBackend end
struct IpoptBackend <: AbstractSolverBackend end
struct MadNLPBackend <: AbstractSolverBackend end

weakdeps = Dict(IpoptBackend => :NLPModelsIpopt, MadNLPBackend => :MadNLP)

function solve_docp(solver_backend::T, args...; kwargs...) where {T <: AbstractSolverBackend}
    throw(ExtensionError(weakdeps[T]))
end
# CTDirect interface

# available methods by order of preference: from top to bottom
algorithms = ()
algorithms = add(algorithms, (:adnlp, :ipopt))

"""
$(TYPEDSIGNATURES)

Return the list of available methods to solve the optimal control problem.
"""
function available_methods()::Tuple{Tuple{Vararg{Symbol}}}
    return algorithms
end


"""
$(TYPEDSIGNATURES)

Discretize an optimal control problem into a nonlinear optimization problem (ie direct transcription)
"""
function direct_transcription(ocp::OptimalControlModel,
    description...;
    init=nothing,
    grid_size::Integer=__grid_size_direct(),
    time_grid=__time_grid_direct())

    # build DOCP
    docp = DOCP(ocp, grid_size, time_grid)
    
    # set bounds in DOCP
    variables_bounds!(docp)
    constraints_bounds!(docp)

    # set initial guess
    x0 = DOCP_initial_guess(docp, OptimalControlInit(init, state_dim=ocp.state_dimension, control_dim=ocp.control_dimension, variable_dim=ocp.variable_dimension))

    # call NLP problem constructor
    nlp = ADNLPModel!(x -> DOCP_objective(x, docp), 
        x0,
        docp.var_l, docp.var_u, 
        (c, x) -> DOCP_constraints!(c, x, docp), 
        docp.con_l, docp.con_u,
        backend = :optimized)

return docp, nlp

end



"""
$(TYPEDSIGNATURES)

Set initial guess in the DOCP
"""
function set_initial_guess(docp::DOCP, nlp, init)

    ocp = docp.ocp
    nlp.meta.x0 .= DOCP_initial_guess(docp, OptimalControlInit(init, state_dim=ocp.state_dimension, control_dim=ocp.control_dimension, variable_dim=ocp.variable_dimension))

end

# placeholders (see CTSolveExt)
#=
function solve_docp(args...; kwargs...)
    error("Please execute `using NLPModelsIpopt` before calling the solve method.")
end
=#

abstract type AbstractSolver end
struct MadNLPSolverSolver <: AbstractSolver end

function solve_docp(s, args...; kwargs...)
    if typeof(s) == IpoptSolver
        error("Please execute `using NLPModelsIpopt` before calling the solve method.")
    elseif typeof(s) == MadNLPSolver
        error("Please execute `using MadNLP` before calling the solve method.")
    else
        error("Unknown solver type", typeof(s))
    end
end


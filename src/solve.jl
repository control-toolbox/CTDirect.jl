# CTDirect interface

# available methods by order of preference: from top to bottom
algorithms = ()
algorithms = add(algorithms, (:adnlp, :ipopt))
algorithms = add(algorithms, (:adnlp, :madnlp))

"""
$(TYPEDSIGNATURES)

Return the list of available methods to solve the optimal control problem.
"""
function available_methods()
    return algorithms
end


"""
$(TYPEDSIGNATURES)

Discretize an optimal control problem into a nonlinear optimization problem (ie direct transcription)
"""
function direct_transcription(ocp::OptimalControlModel,
    description...;
    init=CTBase.__ocp_init(),
    grid_size=__grid_size(),
    time_grid=__time_grid())

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

# placeholders (see CTSolveExt*** extensions)
abstract type SolverTag end
struct IpoptTag <: SolverTag end
struct MadNLPTag <: SolverTag end

function solve_docp(solver_tag, args...; kwargs...)
    if typeof(solver_tag) == IpoptTag
        #error("Please execute `using NLPModelsIpopt` before calling the solve method.")
        throw(ExtensionError("Please execute `using NLPModelsIpopt`"))
    elseif typeof(solver_tag) == MadNLPTag
        #error("Please execute `using MadNLP` before calling the solve method.")
        throw(ExtensionError("Please execute `using MadNLP`"))
    else
        error("Unknown solver type", typeof(solver_tag))
    end
end

# CTDirect interface

# available methods by order of preference: from top to bottom
algorithmes = ()
algorithmes = add(algorithmes, (:adnlp, :ipopt))

"""
$(TYPEDSIGNATURES)

Return the list of available methods to solve the optimal control problem.
"""
function available_methods()::Tuple{Tuple{Vararg{Symbol}}}
    return algorithmes
end


"""
$(TYPEDSIGNATURES)

Discretize an optimal control problem into a nonlinear optimization problem (ie direct transcription)
"""
function directTranscription(ocp::OptimalControlModel,
    description...;
    init=nothing,
    grid_size::Integer=__grid_size_direct(),
    time_grid=nothing)

    # build DOCP
    docp = DOCP(ocp, grid_size, time_grid)

    # set initial guess
    x0 = DOCP_initial_guess(docp, _OptimalControlInit(init, state_dim=ocp.state_dimension, control_dim=ocp.control_dimension, variable_dim=ocp.variable_dimension))

    # set bounds
    docp.var_l, docp.var_u = variables_bounds(docp)
    docp.con_l, docp.con_u = constraints_bounds(docp)

    # call NLP problem constructor
    docp.nlp = ADNLPModel!(x -> DOCP_objective(x, docp), 
                    x0,
                    docp.var_l, docp.var_u, 
                    (c, x) -> DOCP_constraints!(c, x, docp), 
                    docp.con_l, docp.con_u,
                    backend = :optimized)

return docp

end


"""
$(TYPEDSIGNATURES)

Extract the NLP problem from the DOCP
"""
function getNLP(docp::DOCP)
    return docp.nlp
end


"""
$(TYPEDSIGNATURES)

Extract the NLP problem from the DOCP
"""
function setInitialGuess(docp::DOCP, init)

    nlp = getNLP(docp)
    ocp = docp.ocp
    nlp.meta.x0 .= DOCP_initial_guess(docp, _OptimalControlInit(init, state_dim=ocp.state_dimension, control_dim=ocp.control_dimension, variable_dim=ocp.variable_dimension))

end

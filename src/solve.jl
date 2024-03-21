# TODO
# add function to set the intial guess for a docp: need to rebuild the nlp model completely ? 

# availble methods by order of preference: from top to bottom
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
function DirectTranscription(ocp::OptimalControlModel,
    description...;
    init::OptimalControlInit=OptimalControlInit(),
    grid_size::Integer=__grid_size_direct())
    
    # initialization is optional
    docp = DOCP(ocp, grid_size)
    x0 = initial_guess(docp, init)
    l_var, u_var = variables_bounds(docp)
    lb, ub = constraints_bounds(docp)
    docp.nlp = ADNLPModel!(x -> ipopt_objective(x, docp), 
                    x0,
                    l_var, u_var, 
                    (c, x) -> ipopt_constraint!(c, x, docp), 
                    lb, ub, 
                    backend = :optimized)

return docp

end

function getNLP(docp::DOCP)
return docp.nlp
end

"""
$(TYPEDSIGNATURES)

Solve a discretized optimal control problem
"""
function solveDOCP(docp::DOCP;
    display::Bool=__display(),
    #init::OptimalControlInit=OptimalControlInit(),
    print_level::Integer=__print_level_ipopt(),
    mu_strategy::String=__mu_strategy_ipopt(),
    kwargs...)

    # +++ set initial guess if provided
    # setDOCPInitialGuess(+++)

    # solve DOCP with NLP solver
    # sb="yes": remove ipopt header +++ make that default
    print_level = display ?  print_level : 0
    ipopt_solution = ipopt(docp.nlp, print_level=print_level, mu_strategy=mu_strategy, sb="yes"; kwargs...)

    # build OCP solution from DOCP result
    sol = _OptimalControlSolution(ipopt_solution, docp)

    return sol
end

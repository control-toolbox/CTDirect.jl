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
function directTranscription(ocp::OptimalControlModel,
    description...;
    init::OptimalControlInit=OptimalControlInit(),
    grid_size::Integer=__grid_size_direct())
    
    # initialization is optional
    docp = DOCP(ocp, grid_size)
    x0 = initial_guess(docp, init)
    l_var, u_var = variables_bounds(docp)
    lb, ub = constraints_bounds(docp)
    docp.nlp = ADNLPModel!(x -> DOCP_objective(x, docp), 
                    x0,
                    l_var, u_var, 
                    (c, x) -> DOCP_constraint!(c, x, docp), 
                    lb, ub, 
                    backend = :optimized)

return docp

end

function getNLP(docp::DOCP)
    return docp.nlp
end

function setDOCPInit(docp::DOCP, init::OptimalControlInit)
    nlp = getNLP(docp)
    nlp.meta.x0 .= initial_guess(docp, init)
end

"""
$(TYPEDSIGNATURES)

Solve a discretized optimal control problem
"""
function solveDOCP(docp::DOCP;
    init=nothing,
    display::Bool=__display(),
    print_level::Integer=__print_level_ipopt(),
    mu_strategy::String=__mu_strategy_ipopt(),
    kwargs...)

    # solve DOCP with NLP solver
    print_level = display ?  print_level : 0
    if init == nothing
        docp_solution = ipopt(getNLP(docp), print_level=print_level, mu_strategy=mu_strategy, sb="yes"; kwargs...)
    else
        docp_solution = ipopt(getNLP(docp),x0=initial_guess(docp, init), print_level=print_level, mu_strategy=mu_strategy, sb="yes"; kwargs...)
    end

    # return solution for original OCP
    return OCPSolutionFromDOCP(docp, docp_solution)
end


function solveDirect(ocp::OptimalControlModel,
    description...;
    init::OptimalControlInit=OptimalControlInit(),
    grid_size::Integer=__grid_size_direct(),
    display::Bool=__display(),
    print_level::Integer=__print_level_ipopt(),
    mu_strategy::String=__mu_strategy_ipopt(),
    kwargs...)

    # build discretized OCP
    docp = directTranscription(ocp, description, init=init, grid_size=grid_size)
    # solve DOCP and retrieve OCP solution
    ocp_solution = solveDOCP(docp; display=display, print_level=print_level, mu_strategy=mu_strategy, kwargs...)

    return ocp_solution
end
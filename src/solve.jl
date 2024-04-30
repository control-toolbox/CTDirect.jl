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
    init::OCPInit=OCPInit(),
    grid_size::Integer=__grid_size_direct())
    
    # initialization is optional
    docp = DOCP(ocp, grid_size)
    x0 = initial_guess(docp, init)
    docp.var_l, docp.var_u = variables_bounds(docp)
    docp.con_l, docp.con_u = constraints_bounds(docp)
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
function setDOCPInit(docp::DOCP, init::OCPInit)
    nlp = getNLP(docp)
    nlp.meta.x0 .= initial_guess(docp, init)
end


"""
$(TYPEDSIGNATURES)

Solve a discretized optimal control problem DOCP
"""
function solve(docp::DOCP;
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


"""
$(TYPEDSIGNATURES)

Solve an optimal control problem OCP by direct method
"""
function solve(ocp::OptimalControlModel,
    description...;
    init::Union{OCPInit, OptimalControlSolution}=OCPInit(),
    grid_size::Integer=__grid_size_direct(),
    display::Bool=__display(),
    print_level::Integer=__print_level_ipopt(),
    mu_strategy::String=__mu_strategy_ipopt(),
    kwargs...)

    # build init if needed (+++add dictionary as option ?)
    if init isa OptimalControlSolution
        init = OCPInit(init)
    end     

    # build discretized OCP
    docp = directTranscription(ocp, description, init=init, grid_size=grid_size)

    # solve DOCP and retrieve OCP solution
    ocp_solution = solve(docp; display=display, print_level=print_level, mu_strategy=mu_strategy, kwargs...)

    return ocp_solution
end


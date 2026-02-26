# CT packages
using CTBase
using CTParser: CTParser, @def, @init
using CTModels: CTModels, objective, time_grid, iterations, state, control, variable, costate
using CTDirect
using CTSolvers

# other
using CommonSolve

# NLP modelers
using ADNLPModels
using ExaModels

# NLP solvers
using NLPModels
using NLPModelsIpopt
using MadNLPMumps

# misc
using SplitApplyCombine # for flatten in some tests

# solve a given OCP problem (given as struct)
# +++ make more similar to high level solve from Optimal Control (options ...)
function solve_problem(prob;
    max_iter=1000,
    tol=1e-6,
    discretizer=:collocation,
    modeler=:adnlp, #+++exa
    solver=:ipopt, #+++madnlp
    display=false,
    graph=false,
    init=nothing,
    adnlp_backend=:optimized,
    exa_backend=nothing,
    linear_solver=nothing,
    bound_relax_factor=nothing,
    kwargs...)

    # discretized problem (model and solution builders)
    if discretizer == :collocation
        my_discretizer = CTDirect.Collocation(; kwargs...)
    elseif discretizer == :direct_shooting
        my_discretizer = CTDirect.DirectShooting(; kwargs...)
    else
        error("Unknown discretizer: ", discretizer)
    end
    docp = CTDirect.discretize(prob.ocp, my_discretizer)
    
    # initial guess for model builders
    if isnothing(init)
        my_init = CTModels.build_initial_guess(prob.ocp, prob.init)
    else
        my_init = CTModels.build_initial_guess(prob.ocp, init)
    end

    # NLP modeler
    if modeler == :adnlp
        my_modeler = CTSolvers.ADNLP(; backend=adnlp_backend) # kwargs
    elseif modeler == :exa
        my_modeler = CTSolvers.Exa(; backend=exa_backend) # kwargs
    else
        error("Unknown modeler: ", modeler)
    end

    # NLP solver
    if solver == :ipopt
        ipopt_options = Dict(
            :max_iter => max_iter,
            :tol => tol,
            :print_level => 5,
            :mu_strategy => "adaptive",
            :linear_solver => "mumps",
            :sb => "yes",
        )
        my_solver = CTSolvers.Ipopt(; ipopt_options...)
    elseif solver == :madnlp
        if isnothing(linear_solver)
            solver = MumpsSolver
        else
            solver = linear_solver
        end
        madnlp_options = Dict(
            :max_iter => max_iter,
            :tol => tol,
            :linear_solver => solver,
        )
        if !isnothing(bound_relax_factor)
            madnlp_options[:bound_relax_factor] = bound_relax_factor
            madnlp_options[:mode] = :permissive
        end
        my_solver = CTSolvers.MadNLP(; madnlp_options...)
    else
        error("Unknown solver: ", solver)
    end

    # solve DOCP and return OCP solution
    sol = CommonSolve.solve(docp, my_init, my_modeler, my_solver; display=display)

    return sol
end

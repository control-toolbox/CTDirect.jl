# CTDirect interface
using CTBase

"""
$(TYPEDSIGNATURES)

Return the list of available methods to solve the optimal control problem.
"""
function available_methods()
    # available methods by order of preference
    algorithms = ()
    algorithms = CTBase.add(algorithms, (:adnlp, :ipopt))
    algorithms = CTBase.add(algorithms, (:adnlp, :madnlp))
    algorithms = CTBase.add(algorithms, (:adnlp, :knitro))
    algorithms = CTBase.add(algorithms, (:exa, :ipopt))
    algorithms = CTBase.add(algorithms, (:exa, :madnlp))
    algorithms = CTBase.add(algorithms, (:exa, :knitro))
    return algorithms
end

# NLP solver extensions (see ext/CTSolveExt***)
abstract type AbstractNLPSolverBackend end
struct IpoptBackend <: AbstractNLPSolverBackend end
struct MadNLPBackend <: AbstractNLPSolverBackend end
struct KnitroBackend <: AbstractNLPSolverBackend end

weakdeps = Dict(IpoptBackend => :NLPModelsIpopt, MadNLPBackend => :MadNLP, KnitroBackend => :NLPModelsKnitro)

function solve_docp(nlp_solver::T, args...; kwargs...) where {T<:AbstractNLPSolverBackend}
    throw(CTBase.ExtensionError(weakdeps[T]))
end

# NLP model (future) extensions
abstract type AbstractNLPModelBackend end
struct ADNLPBackend <: AbstractNLPModelBackend end
struct ExaBackend <: AbstractNLPModelBackend end

function parse_description(description)

    # default: Ipopt, ADNLPModels
    method = CTBase.complete(description; descriptions=available_methods())

    # get NLP solver choice
    if :ipopt ∈ method
        nlp_solver = CTDirect.IpoptBackend()
    elseif :madnlp ∈ method
        nlp_solver = CTDirect.MadNLPBackend()
    elseif :knitro ∈ method
        nlp_solver = CTDirect.KnitroBackend()
    else
        error("no known solver (:ipopt, :madnlp, :knitro) in method", method)
    end

    # get NLP model choice
    if :adnlp ∈ method
        nlp_model = CTDirect.ADNLPBackend()
    elseif :exa ∈ method
        nlp_model = CTDirect.ExaBackend()
    else
        error("no known model (:adnlp, :exa) in method", method)
    end 

    return nlp_solver, nlp_model
end


"""
$(TYPEDSIGNATURES)

Solve an OCP with a direct method

# Arguments
* ocp: optimal control problem as defined in `CTBase`
* [description]: can specifiy for instance the NLP model and / or solver (:ipopt, :madnlp or :knitro)

# Keyword arguments (optional)
* `display`: ([true], false) will disable output if set to false
* `grid_size`: number of time steps for the discretized problem ([250])
* `disc_method`: discretization method ([`:trapeze`], `:midpoint`, `gauss_legendre_2`)
* `time_grid`: explicit time grid (can be non uniform)
* `init`: info for the starting guess (values or existing solution)
* `nlp_model`: modeller used for optimisation ([`:adnlp`], `:exa`)
* `adnlp_backend`: backend for automatic differentiation in ADNLPModels ([`:optimized`], `:manual`, `:default`)
* `exa_backend`: backend for ExaModels ([`nothing`])

All further keywords are passed to the inner call of `solve_docp`
"""
function solve(
    ocp::CTModels.Model,
    description::Symbol...;
    display::Bool=__display(),
    grid_size::Int=__grid_size(),
    disc_method=__disc_method(),
    time_grid=__time_grid(),
    init=__ocp_init(),
    adnlp_backend=__adnlp_backend(),
    kwargs...,
)

    # build discretized optimal control problem (DOCP)
    # NB. this includes the initial guess for the resulting NLP
    docp = direct_transcription(
        ocp,
        description...;
        init=init,
        grid_size=grid_size,
        time_grid=time_grid,
        disc_method=disc_method,
        kwargs...,
    )

    # get NLP solver choice and solve DOCP
    nlp_solver, nlp_model = parse_description(description)
    docp_solution = CTDirect.solve_docp(nlp_solver, docp; display=display, kwargs...)

    # build and return OCP solution
    return build_OCP_solution(docp, docp_solution)
end


"""
$(TYPEDSIGNATURES)

Discretize an optimal control problem into a nonlinear optimization problem (ie direct transcription)

# Arguments
* ocp: optimal control problem as defined in `CTModels`
* [description]: can specifiy for instance the NLP model and / or solver (:ipopt, :madnlp or :knitro)

# Keyword arguments (optional)
* `grid_size`: number of time steps for the discretized problem ([250])
* `disc_method`: discretization method ([`:trapeze`], `:euler`, `:euler_implicit`, `:midpoint`, `gauss_legendre_2`, `gauss_legendre_3`)
* `time_grid`: explicit time grid (can be non uniform)
* `init`: info for the starting guess (values as named tuple or existing solution)
* `nlp_model`: modeller used for optimisation ([`:adnlp`], `:exa`)
* `adnlp_backend`: backend for automatic differentiation in ADNLPModels ([`:optimized`], `:manual`, `:default`)
* `exa_backend`: tries to solve on GPU whenever possible ([`false`], `true`)
* `show_time`: (:true, [:false]) show timing details from ADNLPModels

"""
function direct_transcription(
    ocp::CTModels.Model,
    description...;
    grid_size=__grid_size(),
    disc_method=__disc_method(),
    time_grid=__time_grid(),
    init=__ocp_init(),
    adnlp_backend=__adnlp_backend(),
    kwargs...,
)

    nlp_solver, nlp_model = parse_description(description)

    # build DOCP
    if nlp_model isa ExaBackend
        docp = DOCP(ocp; grid_size=grid_size, time_grid=time_grid, disc_method=disc_method, lagrange_to_mayer=false)
    else
        docp = DOCP(ocp; grid_size=grid_size, time_grid=time_grid, disc_method=disc_method)
    end

    # set bounds in DOCP
    variables_bounds!(docp)
    constraints_bounds!(docp)

    # build and set initial guess in DOCP
    docp_init = CTModels.Init(
        init;
        state_dim=CTModels.state_dimension(ocp),
        control_dim=CTModels.control_dimension(ocp),
        variable_dim=CTModels.variable_dimension(ocp),
    )
    x0 = DOCP_initial_guess(docp, docp_init)

    # build nlp
    docp.nlp = build_nlp(nlp_model, 
    docp, 
    x0; 
    nlp_solver=nlp_solver, 
    adnlp_backend=adnlp_backend, # for adnlpmodel
    grid_size=grid_size, disc_method=disc_method, # for examodel
    kwargs...)

    return docp
end


function build_nlp(
    nlp_model::ADNLPBackend,
    docp::CTDirect.DOCP,
    x0;
    adnlp_backend=__adnlp_backend(),
    show_time=false, #+default
    matrix_free=false, #+default
    nlp_solver=nothing,
    kwargs...
)

    # redeclare objective and constraints functions
    f = x -> DOCP_objective(x, docp)
    c! = (c, x) -> DOCP_constraints!(c, x, docp)

    # call NLP problem constructor
    if adnlp_backend == :manual

        # build sparsity pattern
        J_backend = ADNLPModels.SparseADJacobian(docp.dim_NLP_variables, f, docp.dim_NLP_constraints, c!, DOCP_Jacobian_pattern(docp))
        H_backend = ADNLPModels.SparseReverseADHessian(docp.dim_NLP_variables, f, docp.dim_NLP_constraints, c!, DOCP_Hessian_pattern(docp))

        # build NLP with given patterns; disable unused backends according to solver info
        if (nlp_solver isa IpoptBackend || nlp_solver isa MadNLPBackend || nlp_solver isa KnitroBackend)
            nlp = ADNLPModel!(
                f,
                x0,
                docp.bounds.var_l, docp.bounds.var_u,
                c!,
                docp.bounds.con_l, docp.bounds.con_u,
                gradient_backend=ADNLPModels.ReverseDiffADGradient,
                jacobian_backend=J_backend,
                hessian_backend=H_backend,
                hprod_backend=ADNLPModels.EmptyADbackend,
                jtprod_backend=ADNLPModels.EmptyADbackend,
                jprod_backend=ADNLPModels.EmptyADbackend,
                ghjvprod_backend=ADNLPModels.EmptyADbackend,
                show_time=show_time,
                #excluded_backend = [:jprod_backend, :jtprod_backend, :hprod_backend, :ghjvprod_backend]
            )
        else
            nlp = ADNLPModel!(
                f,
                x0,
                docp.bounds.var_l, docp.bounds.var_u,
                c!,
                docp.bounds.con_l, docp.bounds.con_u,
                gradient_backend=ADNLPModels.ReverseDiffADGradient,
                jacobian_backend=J_backend,
                hessian_backend=H_backend,
                show_time=show_time,
            )
        end
    else
        # build NLP; disable unused backends according to solver info
        if (nlp_solver isa IpoptBackend || nlp_solver isa MadNLPBackend || nlp_solver isa KnitroBackend)
            nlp = ADNLPModel!(
                f,
                x0,
                docp.bounds.var_l, docp.bounds.var_u,
                c!,
                docp.bounds.con_l, docp.bounds.con_u,
                backend=adnlp_backend,
                hprod_backend=ADNLPModels.EmptyADbackend,
                jtprod_backend=ADNLPModels.EmptyADbackend,
                jprod_backend=ADNLPModels.EmptyADbackend,
                ghjvprod_backend=ADNLPModels.EmptyADbackend,
                show_time=show_time,
            )
        else
            # use manual settings including matrix_free
            nlp = ADNLPModel!(
                f,
                x0,
                docp.bounds.var_l, docp.bounds.var_u,
                c!,
                docp.bounds.con_l, docp.bounds.con_u,
                backend=adnlp_backend,
                show_time=show_time,
                matrix_free=matrix_free
            )
        end
    end

    return nlp
end


function build_nlp(
    nlp_model::CTDirect.ExaBackend,
    docp::CTDirect.DOCP,
    x0;
    grid_size=__grid_size(),
    disc_method=__disc_method(),
    exa_backend=__exa_backend(),
    kwargs...,
)

    # build nlp
    # debug: (time_grid != __time_grid()) || throw("non uniform time grid not available for nlp_model = :exa") # todo: remove when implemented in CTParser
    build_exa = CTModels.get_build_examodel(docp.ocp)
    nlp = build_exa(; grid_size = grid_size, backend = exa_backend, scheme = disc_method) 
    
    # set initial guess
    nlp.meta.x0[1:docp.dim_NLP_variables] .= x0 # NB we currently have an unused final control in examodel

    return nlp
end


"""
$(TYPEDSIGNATURES)

Set initial guess in the DOCP
"""
function set_initial_guess(docp::DOCP, init)
    ocp = docp.ocp
    docp_init = CTModels.Init(
        init;
        state_dim=CTModels.state_dimension(ocp),
        control_dim=CTModels.control_dimension(ocp),
        variable_dim=CTModels.variable_dimension(ocp),
    )
    docp.nlp.meta.x0 .= DOCP_initial_guess(docp, docp_init)
end

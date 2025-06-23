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

# ----------------------------------------------------------------------
# Packages associated to Symbols: used for display
const PACKAGES = Dict(
    # NLP solver
    :ipopt  => :NLPModelsIpopt,
    :madnlp => :MadNLP,
    :knitro => :NLPModelsKnitro,
    # NLP modeller
    :adnlp  => :ADNLPModels,
    :exa    => :ExaModels,
)

# ----------------------------------------------------------------------
# EXTENSIONS

# NLP solver backend extensions 
abstract type AbstractNLPSolverBackend end
struct IpoptBackend <: AbstractNLPSolverBackend end
struct MadNLPBackend <: AbstractNLPSolverBackend end
struct KnitroBackend <: AbstractNLPSolverBackend end

# NLP model (future) backend extensions
abstract type AbstractNLPModelBackend end
struct ADNLPBackend <: AbstractNLPModelBackend end
struct ExaBackend <: AbstractNLPModelBackend end

## Extensions and weak dependencies (see ext/CTDirectExt***)
const WEAKDEPS = Dict{Type, Any}(
    # NLP solver
    IpoptBackend  => [:NLPModelsIpopt],
    MadNLPBackend => [:MadNLP],
    KnitroBackend => [:NLPModelsKnitro],
    # NLP modeller
    ADNLPBackend  => [:ADNLPModels],
    ExaBackend    => [:ExaModels, :MadNLPGPU, :CUDA],
)

# solver
function solve_docp(nlp_solver::T, args...; kwargs...) where {T<:AbstractNLPSolverBackend}
    throw(CTBase.ExtensionError(WEAKDEPS[T]...))
end

# modeller
function build_nlp(nlp_model::T, args...; kwargs...) where {T<:AbstractNLPModelBackend}
    throw(CTBase.ExtensionError(WEAKDEPS[T]...))
end
# ----------------------------------------------------------------------

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

    # display infos about the chosen method
    display && display_method(ocp, description...; 
        grid_size=grid_size,
        time_grid=time_grid,
        disc_method=disc_method,
        kwargs...)

    # build discretized optimal control problem (DOCP)
    # NB. this includes the initial guess for the resulting NLP
    docp = direct_transcription(
        ocp,
        description...;
        init=init,
        grid_size=grid_size,
        time_grid=time_grid,
        disc_method=disc_method,
        adnlp_backend=adnlp_backend,
        kwargs...,
    )

    # get NLP solver choice and solve DOCP
    nlp_solver, nlp_model = parse_description(description)
    docp_solution = CTDirect.solve_docp(nlp_solver, docp; display=display, kwargs...)

    # build and return OCP solution
    return build_OCP_solution(docp, docp_solution)
end

function display_method(ocp, description::Symbol...; grid_size, disc_method, time_grid, kwargs...,)

    # complete description
    method = CTBase.complete(description; descriptions=available_methods())

    #
    print("▫ The optimal control problem is solved with ")
    printstyled("CTDirect", color = :black, bold = true)
    print(" version v$(version()).", "\n\n", "   ┌─ The discretised problem is modelled with ")
    printstyled(PACKAGES[method[1]], color = :black, bold = true)
    print(" and solved with ")
    printstyled(PACKAGES[method[2]], color = :black, bold = true)
    println(".")
    println("   │")

    #
    time = DOCPtime(ocp, grid_size, time_grid)
    N = time.steps

    println("   ├─ Number of time steps⋅: ", N)
    println("   └─ Discretisation scheme: ", disc_method)
	println("")

    # for ipopt
    if !(:print_level ∈ keys(kwargs) && kwargs[:print_level] != 5)
        print("▫ ")
    end

    return nothing
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

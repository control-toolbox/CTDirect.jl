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
    return algorithms
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
* `adnlp_backend`: backend for automatic differentiation in ADNLPModels ([`:optimized`], `:manual`, `:default`)

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

    # get solver choice
    method = CTBase.getFullDescription(description, available_methods())
    if :ipopt ∈ method
        solver_backend = CTDirect.IpoptBackend()
    elseif :madnlp ∈ method
        solver_backend = CTDirect.MadNLPBackend()
    elseif :knitro ∈ method
        solver_backend = CTDirect.KnitroBackend()
    else
        error("no known solver in method", method)
    end

    # build discretized OCP, including initial guess
    docp, nlp = direct_transcription(
        ocp,
        description;
        init=init,
        grid_size=grid_size,
        time_grid=time_grid,
        disc_method=disc_method,
        adnlp_backend=adnlp_backend,
        solver_backend=solver_backend
    )

    # solve DOCP
    docp_solution = CTDirect.solve_docp(solver_backend, docp, nlp; display=display, kwargs...)

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
* `adnlp_backend`: backend for automatic differentiation in ADNLPModels ([`:optimized`], `:manual`, `:default`)
* show_time: (:true, [:false]) show timing details from ADNLPModels

"""
function direct_transcription(
    ocp::CTModels.Model,
    description...;
    grid_size=__grid_size(),
    disc_method=__disc_method(),
    time_grid=__time_grid(),
    init=__ocp_init(),
    adnlp_backend=__adnlp_backend(),
    solver_backend=nothing,
    show_time=false,
    matrix_free=false
)

    # build DOCP
    docp = DOCP(ocp; grid_size=grid_size, time_grid=time_grid, disc_method=disc_method)

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

    # redeclare objective and constraints functions
    f = x -> DOCP_objective(x, docp)
    c! = (c, x) -> DOCP_constraints!(c, x, docp)

    # call NLP problem constructor
    if adnlp_backend == :manual

        # build sparsity pattern
        J_backend = ADNLPModels.SparseADJacobian(docp.dim_NLP_variables, f, docp.dim_NLP_constraints, c!, DOCP_Jacobian_pattern(docp))
        H_backend = ADNLPModels.SparseReverseADHessian(docp.dim_NLP_variables, f, docp.dim_NLP_constraints, c!, DOCP_Hessian_pattern(docp))

        # build NLP with given patterns; disable unused backends according to solver info
        if (solver_backend isa IpoptBackend || solver_backend isa MadNLPBackend || solver_backend isa KnitroBackend)
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
        if (solver_backend isa IpoptBackend || solver_backend isa MadNLPBackend || solver_backend isa KnitroBackend)
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

    return docp, nlp
end


"""
$(TYPEDSIGNATURES)

Set initial guess in the DOCP
"""
function set_initial_guess(docp::DOCP, nlp, init)
    ocp = docp.ocp
    docp_init = CTModels.Init(
        init;
        state_dim=CTModels.state_dimension(ocp),
        control_dim=CTModels.control_dimension(ocp),
        variable_dim=CTModels.variable_dimension(ocp),
    )
    nlp.meta.x0 .= DOCP_initial_guess(docp, docp_init)
end


# placeholders (see CTSolveExt*** extensions)
abstract type AbstractSolverBackend end
struct IpoptBackend <: AbstractSolverBackend end
struct MadNLPBackend <: AbstractSolverBackend end
struct KnitroBackend <: AbstractSolverBackend end

weakdeps = Dict(IpoptBackend => :NLPModelsIpopt, MadNLPBackend => :MadNLP, KnitroBackend => :NLPModelsKnitro)

function solve_docp(solver_backend::T, args...; kwargs...) where {T<:AbstractSolverBackend}
    throw(CTBase.ExtensionError(weakdeps[T]))
end

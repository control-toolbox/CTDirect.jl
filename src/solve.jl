# CTDirect interface

"""
$(TYPEDSIGNATURES)

Return a tuple of available NLP model and solver combinations for solving optimal control problems.

# Returns

- `algorithms::Tuple`: A tuple of symbol pairs representing the available methods.

# Example

```julia-repl
julia> available_methods()
((:adnlp, :ipopt), (:adnlp, :madnlp), (:adnlp, :knitro), (:exa, :ipopt), (:exa, :madnlp), (:exa, :knitro))
```
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
    :ipopt => :NLPModelsIpopt,
    :madnlp => :MadNLPMumps,
    :knitro => :NLPModelsKnitro,
    # NLP modeller
    :adnlp => :ADNLPModels,
    :exa => :ExaModels,
)

# solver
"""
$(TYPEDSIGNATURES)

Solve a discretized optimal control problem using the specified solver backend.

# Arguments

- `solver_backend::T`: An NLP solver backend (subtype of `AbstractNLPSolverBackend`).
- `docp::CTDirect.DOCP`: The discretized optimal control problem.

# Returns

- Throws `CTBase.ExtensionError` if the solver backend is unavailable.

# Example

```julia-repl
julia> solve_docp(IpoptBackend(), docp)
ERROR: ExtensionError(...)
```
"""
function solve_docp(
    solver_backend::T, docp::CTDirect.DOCP; kwargs...
) where {T<:AbstractNLPSolverBackend}
    throw(CTBase.ExtensionError(WEAKDEPS[T]...))
end

# modeller
"""
$(TYPEDSIGNATURES)

Build the NLP model for a discretized optimal control problem using the specified NLP backend.

# Arguments

- `docp::CTDirect.DOCP`: The discretized optimal control problem.
- `nlp_model::T`: The NLP model backend (subtype of `AbstractNLPModelBackend`).
- `x0`: Initial guess for decision variables.

# Returns

- Throws `CTBase.ExtensionError` if the NLP model backend is unavailable.

# Example

```julia-repl
julia> build_nlp!(docp, ADNLPBackend(), x0)
ERROR: ExtensionError(...)
```
"""
function build_nlp!(
    docp::CTDirect.DOCP, nlp_model::T, x0; kwargs...
) where {T<:AbstractNLPModelBackend}
    throw(CTBase.ExtensionError(WEAKDEPS[T]...))
end
# ----------------------------------------------------------------------

"""
$(TYPEDSIGNATURES)

Parse the method description to determine the NLP solver or model.

# Arguments

- `description`: A tuple of symbols representing the desired solver and/or model.
    - NLP solver: `ipopt`, `madnlp` or `knitro` 
    - NLP model: `:adnlp` or `:exa`
- `info::Symbol`: Either `:solver` to return the solver backend or `:model` to return the NLP model backend.

# Returns

- `nlp_solver` or `nlp_model`: The corresponding backend instance.

# Example

```julia-repl
julia> parse_description((:adnlp, :ipopt), :solver)
CTDirect.IpoptBackend()

julia> parse_description((:exa, :madnlp), :model)
CTDirect.ExaBackend()
```
"""
function parse_description(description, info)

    # default: Ipopt, ADNLPModels
    method = CTBase.complete(description; descriptions=available_methods())

    if info == :solver
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

        # patch: replaces ipopt by madnlp for :exa as long as the issue with getters for a posteriori treatment is not fixed
        #=if (:exa ∈ method) && (:ipopt ∈ method)
            nlp_solver = CTDirect.MadNLPBackend()
            @warn "currently replacing Ipopt with MadNLP for :exa"
        end=#
        return nlp_solver

    elseif info == :model
        # get NLP model choice
        if :adnlp ∈ method
            nlp_model = CTDirect.ADNLPBackend()
        elseif :exa ∈ method
            nlp_model = CTDirect.ExaBackend()
        else
            error("no known model (:adnlp, :exa) in method", method)
        end
        return nlp_model
    else
        error("parse_description info should be either :solver or :model, got ", info)
        return nothing
    end
end

"""
$(TYPEDSIGNATURES)

Solve an optimal control problem using a direct transcription method.

# Arguments

- `ocp::CTModels.Model`: The continuous-time optimal control problem.
- `description::Symbol...`: Symbols specifying the NLP model (`:adnlp` or `:exa`) and/or solver (`:ipopt`, `:madnlp`, `:knitro`).

# Keyword Arguments (optional)

- `display::Bool`: Display solver output ([`true`], `false`).
- `grid_size::Int`: Number of discretization steps ([`250`]).
- `disc_method`: Discretization scheme (`:trapeze`, [`:midpoint`], `:gauss_legendre_2`, etc.).
- `time_grid`: Explicit time grid, uniform or not.
- `init`: Initial guess for states, controls, or variables.
- `adnlp_backend`, `exa_backend`: Backend options for NLP modelers.
- `lagrange_to_mayer`: Convert Lagrange cost to Mayer cost
- `kwargs...`: Additional parameters passed to NLP modelers and solvers.

# Returns

- `solution::CTModels.Solution`: The continuous-time solution with objective, state/control trajectories, solver status, and convergence information. Main features:
    - `objective(sol)`: value of the objective
    - `state(sol)`, `control(sol)`: functions for state and control variables (trajectory)
    - `variable(sol)`: optimization variables if any (e.g. free final time)
    - `successful(sol)`: boolean indicating successful convergence of the NLP solver
    - `status(sol)`: symbol for the return code of the NLP solver
    - `message(sol)`: string with specific info from the NLP solver, if any
    - `constraints_violation(sol)`: primal feasibility at the solution
    - `iterations(sol)`: number of iterations 

# Example

```julia-repl
julia> sol = solve(ocp, :adnlp, :ipopt; grid_size=100)
CTModels.Solution(...)
```
"""
function solve(
    ocp::CTModels.Model,
    description::Symbol...;
    display::Bool=__display(),
    grid_size::Int=__grid_size(),
    disc_method=__disc_method(),
    time_grid=__time_grid(),
    init=__ocp_init(),
    lagrange_to_mayer=__lagrange_to_mayer(),
    adnlp_backend=__adnlp_backend(),
    exa_backend=__exa_backend(),
    kwargs...,
)

    # display infos about the chosen method
    display && display_method(
        ocp,
        description...;
        grid_size=grid_size,
        time_grid=time_grid,
        disc_method=disc_method,
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
        lagrange_to_mayer=lagrange_to_mayer,
        adnlp_backend=adnlp_backend,
        exa_backend=exa_backend,
        kwargs...,
    )

    # get NLP solver choice and solve DOCP
    nlp_solver = parse_description(description, :solver)
    nlp_model = parse_description(description, :model)
    docp_solution = CTDirect.solve_docp(nlp_solver, docp; display=display, kwargs...)

    # build and return OCP solution
    return build_OCP_solution(
        docp, docp_solution; nlp_model=nlp_model, nlp_solver=nlp_solver
    )
end

"""
$(TYPEDSIGNATURES)

Display information about the chosen NLP model, solver, discretization scheme, and number of steps.

# Arguments

- `ocp`: The continuous-time optimal control problem.
- `description::Symbol...`: Symbols specifying the solver and model.
- `grid_size::Int`: Number of time steps.
- `disc_method`: Discretization scheme.
- `time_grid`: Optional explicit time grid.
- `kwargs...`: Additional keyword arguments.

# Returns

- `nothing`

# Example

```julia-repl
julia> display_method(ocp, :adnlp, :ipopt; grid_size=100, disc_method=:trapeze)
▫ The optimal control problem is solved with CTDirect version vX.Y.Z.
...
```
"""
function display_method(
    ocp, description::Symbol...; grid_size, disc_method, time_grid, kwargs...
)

    # complete description
    method = CTBase.complete(description; descriptions=available_methods())

    #
    print("▫ The optimal control problem is solved with ")
    printstyled("CTDirect"; color=:black, bold=true)
    print(" version v$(version()).", "\n\n", "   ┌─ The NLP is modelled with ")
    printstyled(PACKAGES[method[1]]; color=:black, bold=true)
    print(" and solved with ")
    printstyled(PACKAGES[method[2]]; color=:black, bold=true)
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

Convert a continuous-time optimal control problem into a discretized nonlinear programming problem.

# Arguments

- `ocp::CTModels.Model`: Continuous-time optimal control problem.
- `description...`: Symbols specifying the NLP model ([`:adnlp`] or `:exa`) and/or solver ([`:ipopt`], :madnlp or :knitro).

# Keyword Arguments (optional)

- `grid_size::Int`: Number of discretization steps ([`250`]).
- `disc_method`: Discretization scheme (`:trapeze`, `:euler`, `:euler_implicit`, [`:midpoint`], `gauss_legendre_2`, `gauss_legendre_3`).
- `time_grid`: Explicit time grid (can be non uniform).
- `init`: Initial guess values or existing solution.
- `lagrange_to_mayer::Bool`: Convert Lagrange cost to Mayer cost (`true` or `false`).
- `kwargs...`: Additional arguments passed to the NLP modeler.

# Returns

- `docp::CTDirect.DOCP`: Discretized optimal control problem ready for NLP solving.

# Example

```julia-repl
julia> docp = direct_transcription(ocp, :adnlp, :ipopt; grid_size=100, disc_method=:trapeze)
CTDirect.DOCP(...)
```
"""
function direct_transcription(
    ocp::CTModels.Model,
    description...;
    grid_size=__grid_size(),
    disc_method=__disc_method(),
    time_grid=__time_grid(),
    init=__ocp_init(),
    lagrange_to_mayer=__lagrange_to_mayer(),
    kwargs...,
)
    nlp_model = parse_description(description, :model)

    # build DOCP
    if nlp_model isa ExaBackend
        docp = DOCP(
            ocp,
            nlp_model;
            grid_size=grid_size,
            time_grid=time_grid,
            disc_method=disc_method,
            lagrange_to_mayer=false,
        )
    else
        docp = DOCP(
            ocp,
            nlp_model;
            grid_size=grid_size,
            time_grid=time_grid,
            disc_method=disc_method,
            lagrange_to_mayer=lagrange_to_mayer,
        )
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
    build_nlp!(
        docp,
        nlp_model, # +++is now in docp, can be removed
        x0;
        grid_size=grid_size,
        disc_method=disc_method,
        kwargs...,
    )

    return docp
end

"""
$(TYPEDSIGNATURES)

Set the initial guess for the decision variables in a discretized optimal control problem.

# Arguments

- `docp::DOCP`: The discretized optimal control problem.
- `init`: Initial guess values as a named tuple or existing solution.

# Returns

- `nothing`

# Example

```julia-repl
julia> set_initial_guess(docp, init)
```
"""
function set_initial_guess(docp::DOCP, init)
    ocp = ocp_model(docp)
    docp_init = CTModels.Init(
        init;
        state_dim=CTModels.state_dimension(ocp),
        control_dim=CTModels.control_dimension(ocp),
        variable_dim=CTModels.variable_dimension(ocp),
    )
    docp.nlp.meta.x0 .= DOCP_initial_guess(docp, docp_init)
end

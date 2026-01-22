# CTDirect interface

# # modeller
# """
# $(TYPEDSIGNATURES)

# Build the NLP model for a discretized optimal control problem using the specified NLP backend.

# # Arguments

# - `docp::CTDirect.DOCP`: The discretized optimal control problem.
# - `nlp_model_backend::T`: The NLP model backend (subtype of `AbstractNLPModelBackend`).
# - `x0`: Initial guess for decision variables.

# # Returns

# - Throws `CTBase.ExtensionError` if the NLP model backend is unavailable.

# # Example

# ```julia-repl
# julia> build_nlp!(docp, ADNLPBackend(), x0)
# ERROR: ExtensionError(...)
# ```
# """
# function build_nlp!(
#     docp::CTDirect.DOCP{<:CTDirect.Discretization,<:CTModels.Model,T}, x0; kwargs...
# ) where {T<:AbstractNLPModelBackend}
#     throw(CTBase.ExtensionError(WEAKDEPS[T]...))
# end
# # ----------------------------------------------------------------------



# """
# $(TYPEDSIGNATURES)

# Solve an optimal control problem using a direct transcription method.

# # Arguments

# - `ocp::CTModels.Model`: The continuous-time optimal control problem.
# - `description::Symbol...`: Symbols specifying the NLP model (`:adnlp` or `:exa`) and/or solver (`:ipopt`, `:madnlp`, `:knitro`).

# # Keyword Arguments (optional)

# - `display::Bool`: Display solver output ([`true`], `false`).
# - `grid_size::Int`: Number of discretization steps ([`250`]).
# - `disc_method`: Discretization scheme (`:trapeze`, [`:midpoint`], `:gauss_legendre_2`, etc.).
# - `time_grid`: Explicit time grid, uniform or not.
# - `init`: Initial guess for states, controls, or variables.
# - `adnlp_backend`, `exa_backend`: Backend options for NLP modelers.
# - `kwargs...`: Additional parameters passed to NLP modelers and solvers.

# # Returns

# - `solution::CTModels.Solution`: The continuous-time solution with objective, state/control trajectories, solver status, and convergence information. Main features:
#     - `objective(sol)`: value of the objective
#     - `state(sol)`, `control(sol)`: functions for state and control variables (trajectory)
#     - `variable(sol)`: optimization variables if any (e.g. free final time)
#     - `successful(sol)`: boolean indicating successful convergence of the NLP solver
#     - `status(sol)`: symbol for the return code of the NLP solver
#     - `message(sol)`: string with specific info from the NLP solver, if any
#     - `constraints_violation(sol)`: primal feasibility at the solution
#     - `iterations(sol)`: number of iterations 

# # Example

# ```julia-repl
# julia> sol = solve(ocp, :adnlp, :ipopt; grid_size=100)
# CTModels.Solution(...)
# ```
# """
# function solve(
#     ocp::CTModels.Model,
#     description::Symbol...;
#     display::Bool=__display(),
#     grid_size::Int=__grid_size(),
#     disc_method=__disc_method(),
#     time_grid=__time_grid(),
#     init=__ocp_init(),
#     adnlp_backend=__adnlp_backend(),
#     exa_backend=__exa_backend(),
#     kwargs...,
# )

#     # display infos about the chosen method
#     display && display_method(
#         ocp,
#         description...;
#         grid_size=grid_size,
#         time_grid=time_grid,
#         disc_method=disc_method,
#         kwargs...,
#     )

#     # build discretized optimal control problem (DOCP)
#     # NB. this includes the initial guess for the resulting NLP
#     docp = direct_transcription(
#         ocp,
#         description...;
#         init=init,
#         grid_size=grid_size,
#         time_grid=time_grid,
#         disc_method=disc_method,
#         adnlp_backend=adnlp_backend,
#         exa_backend=exa_backend,
#         kwargs...,
#     )

#     # get NLP solver choice and solve DOCP
#     nlp_solver_backend = parse_description(description, :solver)
#     nlp_solution = CTDirect.solve_docp(nlp_solver_backend, docp; display=display, kwargs...)

#     # build and return OCP solution
#     return build_OCP_solution(docp, nlp_solution)
# end

# """
# $(TYPEDSIGNATURES)

# Display information about the chosen NLP model, solver, discretization scheme, and number of steps.

# # Arguments

# - `ocp`: The continuous-time optimal control problem.
# - `description::Symbol...`: Symbols specifying the solver and model.
# - `grid_size::Int`: Number of time steps.
# - `disc_method`: Discretization scheme.
# - `time_grid`: Optional explicit time grid.
# - `kwargs...`: Additional keyword arguments.

# # Returns

# - `nothing`

# # Example

# ```julia-repl
# julia> display_method(ocp, :adnlp, :ipopt; grid_size=100, disc_method=:trapeze)
# ▫ The optimal control problem is solved with CTDirect version vX.Y.Z.
# ...
# ```
# """
# function display_method(
#     ocp, description::Symbol...; grid_size, disc_method, time_grid, kwargs...
# )

#     # ----------------------------------------------------------------------
#     # Packages associated to Symbols: used for display
#     PACKAGES = Dict(
#         # NLP solver
#         :ipopt => "NLPModelsIpopt",
#         :madnlp => "MadNLP suite",
#         :knitro => "NLPModelsKnitro",
#         # NLP modeller
#         :adnlp => "ADNLPModels",
#         :exa => "ExaModels",
#     )

#     # complete description
#     method = CTBase.complete(description; descriptions=available_methods())

#     #
#     print("▫ The optimal control problem is solved with ")
#     printstyled("CTDirect"; color=:black, bold=true)
#     print(" version v$(version()).", "\n\n", "   ┌─ The NLP is modelled with ")
#     printstyled(PACKAGES[method[1]]; color=:black, bold=true)
#     print(" and solved with ")
#     printstyled(PACKAGES[method[2]]; color=:black, bold=true)
#     println(".")
#     println("   │")

#     #
#     time = DOCPtime(ocp, grid_size, time_grid)
#     N = time.steps

#     println("   ├─ Number of time steps⋅: ", N)
#     println("   └─ Discretisation scheme: ", disc_method)
#     println("")

#     # for ipopt
#     if !(:print_level ∈ keys(kwargs) && kwargs[:print_level] != 5)
#         print("▫ ")
#     end

#     return nothing
# end

#=
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
    modeler;
    grid_size=__grid_size(),
    disc_method=__disc_method(),
    time_grid=__time_grid(),
    init=__ocp_init(),
    kwargs...,
)

    # build DOCP
    docp = DOCP(
        ocp;
        grid_size=grid_size,
        time_grid=time_grid,
        disc_method=disc_method,
    )

    # set bounds in DOCP
    variables_bounds!(docp)
    constraints_bounds!(docp)

    # build and set initial guess in DOCP
    docp_init = CTModels.build_initial_guess(ocp, init)
    x0 = DOCP_initial_guess(docp, docp_init)

    # build nlp
    nlp = build_nlp_adnlp(docp, x0; grid_size=grid_size, disc_method=disc_method, kwargs...)

    return docp
end
=#

# """
# $(TYPEDSIGNATURES)

# Set the initial guess for the decision variables in a discretized optimal control problem.

# # Arguments

# - `docp::DOCP`: The discretized optimal control problem.
# - `init`: Initial guess values as a named tuple or existing solution.

# # Returns

# - `nothing`

# # Example

# ```julia-repl
# julia> set_initial_guess(docp, init)
# ```
# """
# function set_initial_guess(docp::DOCP, init)
#     ocp = ocp_model(docp)
#     docp_init = CTModels.Init(
#         init;
#         state_dim=CTModels.state_dimension(ocp),
#         control_dim=CTModels.control_dimension(ocp),
#         variable_dim=CTModels.variable_dimension(ocp),
#     )
#     docp.nlp.meta.x0 .= DOCP_initial_guess(docp, docp_init)
# end

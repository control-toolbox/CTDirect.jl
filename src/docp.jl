# Discretized Optimal Control Problem DOCP

"""
$(TYPEDEF)

Abstract type representing a discretization strategy for an optimal
control problem.  

Concrete subtypes of `Discretization` define specific schemes for
transforming a continuous-time problem into a discrete-time
representation suitable for numerical solution.

# Example

```julia-repl
julia> struct MyDiscretization <: Discretization end
MyDiscretization
```
"""
abstract type Discretization end

"""
$(TYPEDEF)

Internal struct holding boolean flags that characterize properties of the 
discretized optimal control problem (DOCP).

# Fields

- `freet0::Bool`: Whether the OCP has a free initial time.
- `freetf::Bool`: Whether the OCP has a free final time.
- `lagrange::Bool`: Whether the OCP includes a Lagrange cost.
- `mayer::Bool`: Whether the OCP includes a Mayer cost.
- `lagrange_to_mayer::Bool`: Whether the Lagrange cost is reformulated as a Mayer term.
- `max::Bool`: Whether the OCP is a maximization problem.

# Example

```julia-repl
julia> DOCPFlags(true, false, true, true, false, false)
DOCPFlags(true, false, true, true, false, false)
```
"""
struct DOCPFlags
    freet0::Bool
    freetf::Bool
    lagrange::Bool
    mayer::Bool
    lagrange_to_mayer::Bool
    max::Bool
end

"""
$(TYPEDSIGNATURES)

Construct a [`DOCPFlags`](@ref) struct from an OCP model.

# Arguments

- `ocp::CTModels.Model`: The optimal control problem model.
- `lagrange_to_mayer::Bool`: Whether to reformulate the Lagrange cost as a Mayer term.

# Returns

- `DOCPFlags`: A struct encoding the problem’s boolean properties.

# Example

```julia-repl
julia> DOCPFlags(ocp, true)
DOCPFlags(false, true, true, false, true, false)
```
"""
function DOCPFlags(ocp::CTModels.Model, lagrange_to_mayer::Bool)
    has_free_initial_time = CTModels.has_free_initial_time(ocp)
    has_free_final_time = CTModels.has_free_final_time(ocp)
    has_lagrange = CTModels.has_lagrange_cost(ocp)
    has_mayer = CTModels.has_mayer_cost(ocp)
    is_maximization = CTModels.criterion(ocp) == :max

    return DOCPFlags(
        has_free_initial_time,
        has_free_final_time,
        has_lagrange,
        has_mayer,
        lagrange_to_mayer,
        is_maximization,
    )
end

"""
$(TYPEDEF)

Internal struct holding problem dimensions for a DOCP.

# Fields

- `NLP_x::Int`: State dimension, possibly including an extra variable for Lagrange cost.
- `NLP_u::Int`: Control dimension.
- `NLP_v::Int`: Variable dimension.
- `OCP_x::Int`: State dimension of the original OCP.
- `path_cons::Int`: Path constraints dimension.
- `boundary_cons::Int`: Boundary constraints dimension.

# Example

```julia-repl
julia> DOCPdims(4, 2, 1, 3, 2, 1)
DOCPdims(4, 2, 1, 3, 2, 1)
```
"""
struct DOCPdims
    NLP_x::Int
    NLP_u::Int
    NLP_v::Int
    OCP_x::Int
    path_cons::Int
    boundary_cons::Int
end

"""
$(TYPEDSIGNATURES)

Construct a [`DOCPdims`](@ref) struct from an OCP model.

# Arguments

- `ocp::CTModels.Model`: The optimal control problem model.
- `lagrange_to_mayer::Bool`: Whether the Lagrange cost is reformulated as Mayer.

# Returns

- `DOCPdims`: A struct containing the problem dimensions.

# Example

```julia-repl
julia> DOCPdims(ocp, true)
DOCPdims(5, 2, 1, 4, 2, 1)
```
"""
function DOCPdims(ocp::CTModels.Model, lagrange_to_mayer::Bool)
    if CTModels.has_lagrange_cost(ocp) && lagrange_to_mayer
        dim_NLP_x = CTModels.state_dimension(ocp) + 1
    else
        dim_NLP_x = CTModels.state_dimension(ocp)
    end
    dim_NLP_u = CTModels.control_dimension(ocp)
    dim_NLP_v = CTModels.variable_dimension(ocp)
    dim_OCP_x = CTModels.state_dimension(ocp)
    dim_path_cons = CTModels.dim_path_constraints_nl(ocp)
    dim_boundary_cons = CTModels.dim_boundary_constraints_nl(ocp)

    return DOCPdims(
        dim_NLP_x, dim_NLP_u, dim_NLP_v, dim_OCP_x, dim_path_cons, dim_boundary_cons
    )
end

"""
$(TYPEDEF)

Internal struct holding time grid information for a DOCP.

# Fields

- `steps::Int`: Number of time steps.
- `normalized_grid::Vector{Float64}`: Normalized time grid in `[0,1]`.
- `fixed_grid::Vector{Float64}`: Fixed time grid in `[t0, tf]`.

# Example

```julia-repl
julia> DOCPtime(10, collect(0:0.1:1), collect(0.0:0.1:1.0))
DOCPtime(10, [0.0, 0.1, …, 1.0], [0.0, 0.1, …, 1.0])
```
"""
struct DOCPtime
    steps::Int
    normalized_grid::Vector{Float64}
    fixed_grid::Vector{Float64}
end

"""
$(TYPEDSIGNATURES)

Construct a [`DOCPtime`](@ref) struct from an OCP model.

# Arguments

- `ocp::CTModels.Model`: The optimal control problem model.
- `grid_size::Int`: Number of grid steps if no grid is provided.
- `time_grid`: Custom time grid (or `nothing` to auto-generate).

# Returns

- `DOCPtime`: A struct encoding the time discretization.

# Example

```julia-repl
julia> DOCPtime(ocp, 10, nothing)
DOCPtime(10, [0.0, 0.1, …, 1.0], [0.0, 0.1, …, 1.0])
```
"""
function DOCPtime(ocp::CTModels.Model, grid_size::Int, time_grid)

    # 1. build/recover normalized time grid
    if time_grid === nothing
        NLP_normalized_time_grid = convert(
            Vector{Float64}, collect(LinRange(0, 1, grid_size + 1))
        )
        dim_NLP_steps = grid_size
    else
        # check strictly increasing
        if !issorted(time_grid; lt=<=)
            throw(ArgumentError("given time grid is not strictly increasing. Aborting..."))
            return nothing
        end
        # normalize input grid if needed
        if (time_grid[1] != 0) || (time_grid[end] != 1)
            #println("INFO: normalizing given time grid...")
            t0 = time_grid[1]
            tf = time_grid[end]
            NLP_normalized_time_grid = (time_grid .- t0) ./ (tf - t0)
        else
            NLP_normalized_time_grid = time_grid
        end
        dim_NLP_steps = length(time_grid) - 1
    end

    # 2. build fixed time grid if needed
    if CTModels.has_free_initial_time(ocp) || CTModels.has_free_final_time(ocp)
        # time grid will be recomputed at each NLP iteration
        NLP_fixed_time_grid = Vector{Float64}(undef, dim_NLP_steps+1)
    else
        # compute time grid once for all
        t0 = CTModels.initial_time(ocp)
        tf = CTModels.final_time(ocp)
        NLP_fixed_time_grid = @. t0 + (NLP_normalized_time_grid * (tf - t0))
    end

    return DOCPtime(dim_NLP_steps, NLP_normalized_time_grid, NLP_fixed_time_grid)
end

"""
$(TYPEDEF)

Internal struct holding variable and constraint bounds for a DOCP.

# Fields

- `var_l::Vector{Float64}`: Lower bounds for NLP variables.
- `var_u::Vector{Float64}`: Upper bounds for NLP variables.
- `con_l::Vector{Float64}`: Lower bounds for NLP constraints.
- `con_u::Vector{Float64}`: Upper bounds for NLP constraints.

# Example

```julia-repl
julia> DOCPbounds([-1.0, -2.0], [1.0, 2.0], [0.0], [0.0])
DOCPbounds([-1.0, -2.0], [1.0, 2.0], [0.0], [0.0])
```
"""
struct DOCPbounds
    var_l::Vector{Float64}
    var_u::Vector{Float64}
    con_l::Vector{Float64}
    con_u::Vector{Float64}
end

"""
$(TYPEDEF)

Struct representing a discretized optimal control problem (DOCP).

# Fields

- `discretization::D`: The discretization scheme.
- `ocp::O`: The original OCP model.
- `nlp_model_backend::N`: The NLP model backend.
- `nlp`: The constructed NLP instance.
- `exa_getter::Union{Nothing,Function}`: Getter for ExaModels if used.
- `flags::DOCPFlags`: Boolean flags describing problem structure.
- `dims::DOCPdims`: Problem dimensions.
- `time::DOCPtime`: Time discretization.
- `bounds::DOCPbounds`: Variable and constraint bounds.
- `dim_NLP_variables::Int`: Number of NLP variables.
- `dim_NLP_constraints::Int`: Number of NLP constraints.

# Example

```julia-repl
julia> DOCP(ocp, nlp_model_backend)
DOCP{...}(...)
```
"""
mutable struct DOCP{
    D<:CTDirect.Discretization,O<:CTModels.Model,N<:CTDirect.AbstractNLPModelBackend
}

    # discretization scheme
    discretization::D

    # OCP
    ocp::O # parametric instead of just qualifying reduces allocations (but not time). Specialization ?

    # NLP
    nlp_model_backend::N
    nlp
    exa_getter::Union{Nothing,Function} # getter for ExaModels (if used)

    # boolean flags
    flags::DOCPFlags

    # dimensions
    dims::DOCPdims

    # time grid
    time::DOCPtime

    # lower and upper bounds for variables and constraints
    bounds::DOCPbounds

    # NLP variables and constraints
    dim_NLP_variables::Int
    dim_NLP_constraints::Int

    # constructor
    function DOCP(
        ocp::CTModels.Model,
        nlp_model_backend::CTDirect.AbstractNLPModelBackend;
        grid_size=__grid_size(),
        time_grid=__time_grid(),
        disc_method=__disc_method(),
        lagrange_to_mayer=__lagrange_to_mayer(),
    )

        # boolean flags
        flags = DOCPFlags(ocp, lagrange_to_mayer)

        # dimensions
        dims = DOCPdims(ocp, lagrange_to_mayer)

        # time grid
        time = DOCPtime(ocp, grid_size, time_grid)

        # discretization method (+++ try to unify this if possible)
        if disc_method == :trapeze
            discretization, dim_NLP_variables, dim_NLP_constraints = CTDirect.Trapeze(
                time.steps,
                dims.NLP_x,
                dims.NLP_u,
                dims.NLP_v,
                dims.path_cons,
                dims.boundary_cons,
            )
        elseif disc_method == :midpoint
            discretization, dim_NLP_variables, dim_NLP_constraints = CTDirect.Midpoint(
                time.steps,
                dims.NLP_x,
                dims.NLP_u,
                dims.NLP_v,
                dims.path_cons,
                dims.boundary_cons,
            )
        elseif disc_method == :euler ||
            disc_method == :euler_explicit ||
            disc_method == :euler_forward
            discretization, dim_NLP_variables, dim_NLP_constraints = CTDirect.Euler(
                time.steps,
                dims.NLP_x,
                dims.NLP_u,
                dims.NLP_v,
                dims.path_cons,
                dims.boundary_cons,
            )
        elseif disc_method == :euler_implicit || disc_method == :euler_backward
            discretization, dim_NLP_variables, dim_NLP_constraints = CTDirect.Euler(
                time.steps,
                dims.NLP_x,
                dims.NLP_u,
                dims.NLP_v,
                dims.path_cons,
                dims.boundary_cons;
                explicit=false,
            )
        elseif disc_method == :gauss_legendre_2
            discretization, dim_NLP_variables, dim_NLP_constraints = CTDirect.Gauss_Legendre_2(
                time.steps,
                dims.NLP_x,
                dims.NLP_u,
                dims.NLP_v,
                dims.path_cons,
                dims.boundary_cons,
            )
        elseif disc_method == :gauss_legendre_3
            discretization, dim_NLP_variables, dim_NLP_constraints = CTDirect.Gauss_Legendre_3(
                time.steps,
                dims.NLP_x,
                dims.NLP_u,
                dims.NLP_v,
                dims.path_cons,
                dims.boundary_cons,
            )
        else
            error(
                "Unknown discretization method: ",
                disc_method,
                "\nValid options are disc_method={:trapeze, :midpoint, :euler | :euler_explicit | :euler_forward, :euler_implicit | :euler_backward, :gauss_legendre_2, :gauss_legendre_3}\n",
                typeof(disc_method),
            )
        end

        # add initial condition for lagrange state
        if flags.lagrange && flags.lagrange_to_mayer
            dim_NLP_constraints += 1
        end

        # lower and upper bounds for variables and constraints
        bounds = DOCPbounds(
            -Inf * ones(dim_NLP_variables),
            Inf * ones(dim_NLP_variables),
            zeros(dim_NLP_constraints),
            zeros(dim_NLP_constraints),
        )

        # call constructor with const fields
        docp = new{typeof(discretization),typeof(ocp),typeof(nlp_model_backend)}(
            discretization,
            ocp,
            nlp_model_backend,
            nothing, # nlp
            nothing, # exa_getter
            flags,
            dims,
            time,
            bounds,
            dim_NLP_variables,
            dim_NLP_constraints,
        )

        return docp
    end
end

"""
$(TYPEDSIGNATURES)

Check if an OCP is solvable by [`solve`](@ref).

# Arguments

- `ocp`: The OCP model.

# Returns

- `solvable::Bool`: Always returns `true` in the current implementation.

# Example

```julia-repl
julia> is_solvable(ocp)
true
```
"""
function is_solvable(ocp)
    solvable = true
    return solvable
end

"""
$(TYPEDSIGNATURES)

Build lower and upper bounds vectors for the nonlinear constraints of a DOCP.

# Arguments

- `docp::DOCP`: The discretized OCP.

# Returns

- `(lb, ub)::Tuple{Vector{Float64},Vector{Float64}}`: Lower and upper bounds.

# Example

```julia-repl
julia> constraints_bounds!(docp)
([-1.0, …], [1.0, …])
```
"""
function constraints_bounds!(docp::DOCP)
    lb = docp.bounds.con_l
    ub = docp.bounds.con_u
    disc = disc_model(docp)

    offset = 0
    for i in 1:(docp.time.steps + 1)
        if i <= docp.time.steps
            # skip (ie leave 0) for state / stage equations 
            offset = offset + disc._state_stage_eqs_block
        end
        # path constraints
        if docp.dims.path_cons > 0
            lb[(offset + 1):(offset + docp.dims.path_cons)] = CTModels.path_constraints_nl(
                ocp_model(docp)
            )[1]
            ub[(offset + 1):(offset + docp.dims.path_cons)] = CTModels.path_constraints_nl(
                ocp_model(docp)
            )[3]
            offset = offset + docp.dims.path_cons
        end
    end

    # boundary constraints
    if docp.dims.boundary_cons > 0
        lb[(offset + 1):(offset + docp.dims.boundary_cons)] = CTModels.boundary_constraints_nl(
            ocp_model(docp)
        )[1]
        ub[(offset + 1):(offset + docp.dims.boundary_cons)] = CTModels.boundary_constraints_nl(
            ocp_model(docp)
        )[3]
        offset = offset + docp.dims.boundary_cons
    end

    # null initial condition for lagrangian cost state
    if docp.flags.lagrange && docp.flags.lagrange_to_mayer
        lb[offset + 1] = 0.0
        ub[offset + 1] = 0.0
        offset = offset + 1
    end

    return lb, ub
end

"""
$(TYPEDSIGNATURES)

Build lower and upper bounds vectors for the variable box constraints of a DOCP.

# Arguments

- `docp::DOCP`: The discretized OCP.

# Returns

- `(var_l, var_u)::Tuple{Vector{Float64},Vector{Float64}}`: Lower and upper bounds for variables.

# Example

```julia-repl
julia> variables_bounds!(docp)
([-Inf, …], [Inf, …])
```
"""
function variables_bounds!(docp::DOCP)
    N = docp.time.steps
    var_l = docp.bounds.var_l
    var_u = docp.bounds.var_u
    ocp = ocp_model(docp)

    # build full ordered sets of bounds
    x_lb, x_ub = build_bounds(
        docp.dims.OCP_x,
        CTModels.dim_state_constraints_box(ocp),
        CTModels.state_constraints_box(ocp),
    )
    u_lb, u_ub = build_bounds(
        docp.dims.NLP_u,
        CTModels.dim_control_constraints_box(ocp),
        CTModels.control_constraints_box(ocp),
    )

    # set state / control box along time steps
    for i in 1:(N + 1)
        set_state_at_time_step!(var_l, x_lb, docp, i)
        set_state_at_time_step!(var_u, x_ub, docp, i)
        set_control_at_time_step!(var_l, u_lb, docp, i)
        set_control_at_time_step!(var_u, u_ub, docp, i)
    end

    # variable box
    if docp.dims.NLP_v > 0
        v_lb, v_ub = build_bounds(
            docp.dims.NLP_v,
            CTModels.dim_variable_constraints_box(ocp),
            CTModels.variable_constraints_box(ocp),
        )
        set_optim_variable!(var_l, v_lb, docp)
        set_optim_variable!(var_u, v_ub, docp)
    end

    return var_l, var_u
end

"""
$(TYPEDSIGNATURES)

Compute the objective value of a discretized OCP.

# Arguments

- `xu`: Vector of NLP decision variables.
- `docp::DOCP`: The discretized OCP.

# Returns

- `obj::Float64`: Objective function value.

# Example

```julia-repl
julia> DOCP_objective(xu, docp)
12.34
```
"""
function DOCP_objective(xu, docp::DOCP)

    # initialization
    if docp.flags.freet0 || docp.flags.freetf
        time_grid = get_time_grid(xu, docp)
    else
        time_grid = docp.time.fixed_grid
    end
    v = get_OCP_variable(xu, docp)
    ocp = ocp_model(docp)

    # mayer cost
    if docp.flags.mayer
        x0 = get_OCP_state_at_time_step(xu, docp, 1)
        xf = get_OCP_state_at_time_step(xu, docp, docp.time.steps+1)
        obj_mayer = CTModels.mayer(ocp)(x0, xf, v)
    else
        obj_mayer = 0.0
    end

    # lagrange cost
    if docp.flags.lagrange
        if docp.flags.lagrange_to_mayer
            obj_lagrange = get_lagrange_state_at_time_step(xu, docp, docp.time.steps+1)
        else
            obj_lagrange = runningCost(docp, xu, v, time_grid)
        end
    else
        obj_lagrange = 0.0
    end

    # total cost
    obj = obj_mayer + obj_lagrange

    #= maximization problem
    # +++ add a max_to_min flag in DOCP
    # option minimize[=true] for adnlpmodels
    if docp.flags.max
        obj = -obj
    end=#

    return obj
end

"""
$(TYPEDSIGNATURES)

Compute the nonlinear constraints of a DOCP.

The constraints are modeled as `lb <= C(x) <= ub`.

# Arguments

- `c`: Preallocated constraint vector.
- `xu`: Vector of NLP decision variables.
- `docp::DOCP`: The discretized OCP.

# Returns

- `c`: The filled constraint vector.

# Example

```julia-repl
julia> DOCP_constraints!(zeros(docp.dim_NLP_constraints), xu, docp)
[0.0, 0.1, …]
```
"""
function DOCP_constraints!(c, xu, docp::DOCP)

    # initialization
    if docp.flags.freet0 || docp.flags.freetf
        time_grid = get_time_grid(xu, docp)
    else
        time_grid = docp.time.fixed_grid
    end
    v = get_OCP_variable(xu, docp)
    work = setWorkArray(docp, xu, time_grid, v)

    # main loop on time steps
    for i in 1:(docp.time.steps + 1)
        setStepConstraints!(docp, c, xu, v, time_grid, i, work)
    end

    # point constraints
    setPointConstraints!(docp, c, xu, v)

    # NB. the function *needs* to return c for AD...
    return c
end

"""
$(TYPEDSIGNATURES)

Set boundary and variable point constraints for a DOCP.

# Arguments

- `docp::DOCP`: The discretized OCP.
- `c`: Constraint vector to modify.
- `xu`: Vector of NLP decision variables.
- `v`: Additional OCP variables.

# Returns

- `nothing`: Modifies `c` in place.

# Example

```julia-repl
julia> setPointConstraints!(docp, c, xu, v)
```
"""
function setPointConstraints!(docp::DOCP, c, xu, v)

    # add lagrange state null initial condition at the end
    if docp.flags.lagrange && docp.flags.lagrange_to_mayer
        c[docp.dim_NLP_constraints] = get_lagrange_state_at_time_step(xu, docp, 1)
        offset = docp.dim_NLP_constraints - docp.dims.boundary_cons - 1
    else
        offset = docp.dim_NLP_constraints - docp.dims.boundary_cons
    end

    # boundary constraints
    if docp.dims.boundary_cons > 0
        ocp = ocp_model(docp)
        x0 = get_OCP_state_at_time_step(xu, docp, 1)
        xf = get_OCP_state_at_time_step(xu, docp, docp.time.steps+1)
        CTModels.boundary_constraints_nl(ocp)[2](
            (@view c[(offset + 1):(offset + docp.dims.boundary_cons)]), x0, xf, v
        )
    end
end

"""
$(TYPEDSIGNATURES)

Build an initial guess vector for the discretized OCP.

# Arguments

- `docp::DOCP`: The discretized OCP.
- `init::CTModels.Init`: Initialization settings (default: `CTModels.Init()`).

# Returns

- `NLP_X::Vector{Float64}`: Initial guess vector.

# Example

```julia-repl
julia> DOCP_initial_guess(docp)
[0.1, 0.1, …]
```
"""
function DOCP_initial_guess(docp::DOCP, init::CTModels.Init=CTModels.Init())

    # default initialization (internal variables such as lagrange cost, k_i for RK schemes) will keep these default values 
    NLP_X = 0.1 * ones(docp.dim_NLP_variables)

    # set variables if provided (needed first in case of free times !)
    if !isnothing(init.variable_init)
        set_optim_variable!(NLP_X, init.variable_init, docp)
    end

    # set state / control variables if provided (final control case handled by setter)
    time_grid = get_time_grid(NLP_X, docp)
    for i in 1:(docp.time.steps + 1)
        ti = time_grid[i]
        set_state_at_time_step!(NLP_X, init.state_init(ti), docp, i)
        set_control_at_time_step!(NLP_X, init.control_init(ti), docp, i)
    end

    return NLP_X
end

"""
$(TYPEDSIGNATURES)

Return the time grid for problems with free initial or final times.

# Arguments

- `xu`: Vector of NLP decision variables.
- `docp::DOCP`: The discretized OCP.

# Returns

- `grid::Vector{Float64}`: Time grid corresponding to current NLP variables.

# Example

```julia-repl
julia> get_time_grid(xu, docp)
[0.0, 0.1, …, 1.0]
```
"""
function get_time_grid(xu, docp::DOCP)
    grid = similar(xu, docp.time.steps+1)
    ocp = ocp_model(docp)

    if docp.flags.freet0 || docp.flags.freetf
        v = get_OCP_variable(xu, docp)
    end

    if docp.flags.freet0
        t0 = CTModels.initial_time(ocp, v)
    else
        t0 = CTModels.initial_time(ocp)
    end
    if docp.flags.freetf
        tf = CTModels.final_time(ocp, v)
    else
        tf = CTModels.final_time(ocp)
    end

    @. grid = t0 + docp.time.normalized_grid * (tf - t0)
    return grid
end

"""
$(TYPEDSIGNATURES)

Helper for invalid execution paths.

Always throws an error.

# Example

```julia-repl
julia> get_time_grid_exa()
ERROR: you should not be here
```
"""
function get_time_grid_exa()
    error("you should not be here")
    return nothing
end

"""
$(TYPEDSIGNATURES)

Build lower and upper bound vectors for state, control, or optimization variables.

# Arguments

- `dim_var::Int`: Variable dimension.
- `dim_box::Int`: Number of box constraints.
- `box_triplet`: Triplet defining box constraints.

# Returns

- `(x_lb, x_ub)::Tuple{Vector{Float64},Vector{Float64}}`: Lower and upper bounds.

# Example

```julia-repl
julia> build_bounds(3, 1, ([0.0], [2], [1.0]))
([-Inf, 0.0, -Inf], [Inf, 1.0, Inf])
```
"""
function build_bounds(dim_var, dim_box, box_triplet)
    x_lb = -Inf * ones(dim_var)
    x_ub = Inf * ones(dim_var)
    for j in 1:(dim_box)
        indice = box_triplet[2][j]
        x_lb[indice] = box_triplet[1][j]
        x_ub[indice] = box_triplet[3][j]
    end

    return x_lb, x_ub
end

# getters for high level structs
"""
$(TYPEDSIGNATURES)

Return the nonlinear programming (NLP) model associated with a given
discretized optimal control problem (`DOCP`).

# Arguments

- `docp::DOCP`: The discretized optimal control problem.

# Returns

- `nlp::Any`: The underlying NLP model stored in `docp`.

# Example

```julia-repl
julia> nlp_model(docp)
NLPModel(...)
```
"""
nlp_model(docp::DOCP) = docp.nlp

"""
$(TYPEDSIGNATURES)

Return the continuous-time optimal control problem (OCP) model
associated with a given discretized optimal control problem (`DOCP`).

# Arguments

- `docp::DOCP`: The discretized optimal control problem.

# Returns

- `ocp::Any`: The underlying OCP model stored in `docp`.

# Example

```julia-repl
julia> ocp_model(docp)
OCPModel(...)
```
"""
ocp_model(docp::DOCP) = docp.ocp

"""
$(TYPEDSIGNATURES)

Return the discretization model associated with a given discretized
optimal control problem (`DOCP`).

# Arguments

- `docp::DOCP`: The discretized optimal control problem.

# Returns

- `discretization::Any`: The discretization model stored in `docp`.

# Example

```julia-repl
julia> disc_model(docp)
DiscretizationModel(...)
```
"""
disc_model(docp::DOCP) = docp.discretization

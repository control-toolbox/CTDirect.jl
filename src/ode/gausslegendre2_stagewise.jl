# ============================================================================
# src/ode/irk_gl2_stagewise.jl
#
# Minimal specialized stagewise version of Gauss-Legendre 2 collocation.
#
# Internal layout for NLP variables:
# [X_0, U_0^1, U_0^2, K_0^1, K_0^2,
#  X_1, U_1^1, U_1^2, K_1^1, K_1^2,
#  ...
#  X_N-1, U_N-1^1, U_N-1^2, K_N-1^1, K_N-1^2,
#  X_N, V]
# ============================================================================

abstract type GenericIRKStagewise <: Scheme end

"""
Implicit Gauss-Legendre collocation for s=2 with stagewise controls.
"""
struct Gauss_Legendre_2_Stagewise <: GenericIRKStagewise
    info::String
    stage::Int
    butcher_a::Matrix{Float64}
    butcher_b::Vector{Float64}
    butcher_c::Vector{Float64}
    _control_block::Int
    _step_variables_block::Int
    _state_stage_eqs_block::Int
    _step_pathcons_block::Int
    _final_control::Bool

    function Gauss_Legendre_2_Stagewise(dims::DOCPdims, time::DOCPtime)
        control_block,
        step_variables_block,
        state_stage_eqs_block,
        step_pathcons_block,
        dim_NLP_variables,
        dim_NLP_constraints = GL2_stagewise_dims(dims, time)

        disc = new(
            "Implicit Gauss-Legendre collocation for s=2, 4th order, stagewise controls",
            2,
            [0.25 (0.25 - sqrt(3) / 6);
             (0.25 + sqrt(3) / 6) 0.25],
            [0.5, 0.5],
            [0.5 - sqrt(3) / 6, 0.5 + sqrt(3) / 6],
            control_block,
            step_variables_block,
            state_stage_eqs_block,
            step_pathcons_block,
            false,
        )

        return disc, dim_NLP_variables, dim_NLP_constraints
    end
end

"""
Dimensions for the specialized Gauss-Legendre-2 stagewise scheme.
"""
function GL2_stagewise_dims(dims::DOCPdims, time::DOCPtime)
    stage = 2

    control_block = dims.NLP_u * time.control_steps

    # per step:
    # x_i + U_i^1 + U_i^2 + K_i^1 + K_i^2
    step_variables_block = dims.NLP_x + 2 * control_block + stage * dims.NLP_x

    # state equation + 2 stage equations
    state_stage_eqs_block = dims.NLP_x * (1 + stage)

    # same convention as irk.jl: path constraints at time steps
    step_pathcons_block = dims.path_cons

    dim_NLP_variables = time.steps * step_variables_block + dims.NLP_x + dims.NLP_v

    dim_NLP_constraints =
        time.steps * (state_stage_eqs_block + step_pathcons_block) +
        step_pathcons_block +
        dims.boundary_cons

    return control_block,
           step_variables_block,
           state_stage_eqs_block,
           step_pathcons_block,
           dim_NLP_variables,
           dim_NLP_constraints
end

# ----------------------------------------------------------------------------
# Indexing helpers
# ----------------------------------------------------------------------------

"""
Return control at stage j on time step i.
Convention: 1 <= i <= time.steps, j in {1,2}.
"""
function get_stagecontrol_at_time_step(xu, docp::DOCP{<:Gauss_Legendre_2_Stagewise}, i::Int, j::Int)
    disc = disc_model(docp)
    dims = docp.dims

    @assert 1 <= j <= 2

    var_offset = (i - 1) * disc._step_variables_block
    x_block_end = var_offset + dims.NLP_x
    ctrl0 = x_block_end + (j - 1) * disc._control_block + 1
    ctrl1 = x_block_end + j * disc._control_block

    return @view xu[ctrl0:ctrl1]
end

"""
Compatibility accessor: return a "control at time step".
By default we return U_i^1.
"""
function get_OCP_control_at_time_step(
    xu,
    docp::DOCP{<:Gauss_Legendre_2_Stagewise},
    i::Int,
    j::Int,
)
    if i == docp.time.steps + 1
        i = docp.time.steps
    end

    return get_stagecontrol_at_time_step(xu, docp, i, j)
end

function get_OCP_control_at_time_step(xu, docp::DOCP{<:Gauss_Legendre_2_Stagewise}, i::Int)
    return get_OCP_control_at_time_step(xu, docp, i, 1)
end

"""
Return stage variables K_i^j on time step i.
Convention: 1 <= i <= time.steps, j in {1,2}.
"""
function get_stagevars_at_time_step(xu, docp::DOCP{<:Gauss_Legendre_2_Stagewise}, i::Int, j::Int)
    disc = disc_model(docp)
    dims = docp.dims

    @assert 1 <= j <= 2

    var_offset = (i - 1) * disc._step_variables_block
    stage_block_start = var_offset + dims.NLP_x + 2 * disc._control_block
    k0 = stage_block_start + (j - 1) * dims.NLP_x + 1
    k1 = stage_block_start + j * dims.NLP_x

    return @view xu[k0:k1]
end

# ----------------------------------------------------------------------------
# Work array
# ----------------------------------------------------------------------------

"""
Set work array for dynamics and lagrange cost evaluations.
Layout:
- [x_i^j ; sum_bk]
"""
function setWorkArray(docp::DOCP{<:Gauss_Legendre_2_Stagewise}, xu, time_grid, v)
    dims = docp.dims
    work = similar(xu, dims.NLP_x + dims.NLP_x)
    return work
end

# ----------------------------------------------------------------------------
# Bounds + initial guess
# ----------------------------------------------------------------------------

function __variables_bounds!(docp::DOCP{<:Gauss_Legendre_2_Stagewise})

    var_l = docp.bounds.var_l
    var_u = docp.bounds.var_u
    ocp = ocp_model(docp)

    x_lb, x_ub = build_bounds_block(
        docp.dims.NLP_x,
        CTModels.dim_state_constraints_box(ocp),
        CTModels.state_constraints_box(ocp),
    )
    u_lb, u_ub = build_bounds_block(
        docp.dims.NLP_u,
        CTModels.dim_control_constraints_box(ocp),
        CTModels.control_constraints_box(ocp),
    )

    for i in 1:(docp.time.steps + 1)
        set_state_at_time_step!(var_l, x_lb, docp, i)
        set_state_at_time_step!(var_u, x_ub, docp, i)
    end

    if docp.dims.NLP_u > 0
        for i in 1:docp.time.steps
            for j in 1:2
                get_stagecontrol_at_time_step(var_l, docp, i, j) .= u_lb
                get_stagecontrol_at_time_step(var_u, docp, i, j) .= u_ub
            end
        end
    end

    if docp.dims.NLP_v > 0
        v_lb, v_ub = build_bounds_block(
            docp.dims.NLP_v,
            CTModels.dim_variable_constraints_box(ocp),
            CTModels.variable_constraints_box(ocp),
        )
        set_optim_variable!(var_l, v_lb, docp)
        set_optim_variable!(var_u, v_ub, docp)
    end

    return var_l, var_u
end

function __initial_guess(
    docp::DOCP{<:Gauss_Legendre_2_Stagewise},
    init::CTModels.InitialGuess,
)
    NLP_X = 0.1 * ones(docp.dim_NLP_variables)
    disc = disc_model(docp)

    set_optim_variable!(NLP_X, init.variable, docp)

    time_grid = get_time_grid(NLP_X, docp)
    for i in 1:(docp.time.steps + 1)
        ti = time_grid[i]
        set_state_at_time_step!(NLP_X, init.state(ti), docp, i)
    end

    if docp.dims.NLP_u > 0
        for i in 1:docp.time.steps
            ti = time_grid[i]
            hi = time_grid[i + 1] - ti

            for j in 1:disc.stage
                tij = ti + disc.butcher_c[j] * hi
                uij = init.control(tij)

                if !isnothing(uij)
                    get_stagecontrol_at_time_step(NLP_X, docp, i, j) .= uij
                end
            end
        end
    end

    return NLP_X
end

# ----------------------------------------------------------------------------
# Running cost
# ----------------------------------------------------------------------------

"""
Compute running cost using Gauss quadrature and stagewise controls.
"""
function integral(docp::DOCP{<:Gauss_Legendre_2_Stagewise}, xu, v, time_grid, f)
    disc = disc_model(docp)
    dims = docp.dims

    value = 0.0
    work_xij = similar(xu, dims.NLP_x)

    for i in 1:docp.time.steps
        ti = time_grid[i]
        xi = get_OCP_state_at_time_step(xu, docp, i)
        hi = time_grid[i + 1] - ti

        local_sum = 0.0

        for j in 1:2
            tij = ti + disc.butcher_c[j] * hi
            uij = get_stagecontrol_at_time_step(xu, docp, i, j)

            @views @. work_xij = xi
            for l in 1:2
                kil = get_stagevars_at_time_step(xu, docp, i, l)
                @views @. work_xij = work_xij + hi * disc.butcher_a[j, l] * kil
            end

            local_sum += disc.butcher_b[j] * f(tij, work_xij, uij, v)
        end

        value += hi * local_sum
    end

    return value
end

# ----------------------------------------------------------------------------
# Dynamics + stage constraints
# ----------------------------------------------------------------------------

"""
Set state and stage constraints.
Convention: 1 <= i <= docp.time.steps
"""
function stepStateConstraints!(
    docp::DOCP{<:Gauss_Legendre_2_Stagewise},
    c,
    xu,
    v,
    time_grid,
    i,
    work,
)
    ocp = ocp_model(docp)
    disc = disc_model(docp)
    dims = docp.dims

    # work layout: [x_ij ; sum_bk]
    work_xij = @view work[1:dims.NLP_x]
    work_sumbk = @view work[(dims.NLP_x + 1):(2 * dims.NLP_x)]

    offset = (i - 1) * (disc._state_stage_eqs_block + disc._step_pathcons_block)

    ti = time_grid[i]
    tip1 = time_grid[i + 1]
    hi = tip1 - ti

    xi = get_OCP_state_at_time_step(xu, docp, i)
    xip1 = get_OCP_state_at_time_step(xu, docp, i + 1)

    offset_stage_eqs = dims.NLP_x

    for j in 1:2
        tij = ti + disc.butcher_c[j] * hi
        kij = get_stagevars_at_time_step(xu, docp, i, j)
        uij = get_stagecontrol_at_time_step(xu, docp, i, j)

        if j == 1
            @views @. work_sumbk = disc.butcher_b[j] * kij
        else
            @views @. work_sumbk = work_sumbk + disc.butcher_b[j] * kij
        end

        @views @. work_xij = xi
        for l in 1:2
            kil = get_stagevars_at_time_step(xu, docp, i, l)
            @views @. work_xij = work_xij + hi * disc.butcher_a[j, l] * kil
        end

        CTModels.dynamics(ocp)(
            @view(c[(offset + offset_stage_eqs + 1):(offset + offset_stage_eqs + dims.NLP_x)]),
            tij,
            work_xij,
            uij,
            v,
        )

        @views @. c[(offset + offset_stage_eqs + 1):(offset + offset_stage_eqs + dims.NLP_x)] =
            kij - c[(offset + offset_stage_eqs + 1):(offset + offset_stage_eqs + dims.NLP_x)]

        offset_stage_eqs += dims.NLP_x
    end

    @views @. c[(offset + 1):(offset + dims.NLP_x)] =
        xip1 - (xi + hi * work_sumbk)

    return nothing
end

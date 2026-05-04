# ============================================================================
# src/ode/irk_stagewise.jl
#
# Generic stagewise version of implicit Runge-Kutta collocation.
#
# Internal layout for NLP variables (vector xu):
# [X_0, U_0^1, ..., U_0^s, K_0^1, ..., K_0^s,
#  X_1, U_1^1, ..., U_1^s, K_1^1, ..., K_1^s,
#  ...
#  X_N-1, U_N-1^1, ..., U_N-1^s, K_N-1^1, ..., K_N-1^s,
#  X_N, V]
#
# Internal layout for NLP constraints (vector c):
# [C_0^x, C_0^{k,1}, ..., C_0^{k,s}, G_0,
#  C_1^x, C_1^{k,1}, ..., C_1^{k,s}, G_1,
#  ...
#  C_{N-1}^x, C_{N-1}^{k,1}, ..., C_{N-1}^{k,s}, G_{N-1},
#  G_N,
#  B]
#
# with:
# C_i^x      = state constraint on step i
#              X_{i+1} - (X_i + h_i * sum_j b_j K_i^j)
#
# C_i^{k,j}  = stage-j constraint on step i
#              K_i^j - f(t_i^j, X_i^j, U_i^j, V)
#
# G_i        = path constraints at time t_i   (if path_cons > 0)
# G_N        = path constraints at final time (if path_cons > 0)
# B          = boundary / terminal constraints

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
        stage = 2
        control_block,
        step_variables_block,
        state_stage_eqs_block,
        step_pathcons_block,
        dim_NLP_variables,
        dim_NLP_constraints = IRK_stagewise_dims(dims, time, stage)

        disc = new(
            "Implicit Gauss-Legendre collocation for s=2, 4th order, stagewise controls",
            stage,
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
Implicit Gauss-Legendre collocation for s=3 with stagewise controls.
"""
struct Gauss_Legendre_3_Stagewise <: GenericIRKStagewise
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

    function Gauss_Legendre_3_Stagewise(dims::DOCPdims, time::DOCPtime)
        stage = 3
        control_block,
        step_variables_block,
        state_stage_eqs_block,
        step_pathcons_block,
        dim_NLP_variables,
        dim_NLP_constraints = IRK_stagewise_dims(dims, time, stage)

        disc = new(
            "Implicit Gauss-Legendre collocation for s=3, 6th order, stagewise controls",
            stage,
            [
                (5.0 / 36.0) (2 / 9 - sqrt(15) / 15) (5 / 36 - sqrt(15) / 30);
                (5.0 / 36.0 + sqrt(15) / 24) (2.0 / 9.0) (5.0 / 36.0 - sqrt(15) / 24);
                (5 / 36 + sqrt(15) / 30) (2 / 9 + sqrt(15) / 15) (5.0 / 36.0)
            ],
            [5.0 / 18.0, 4.0 / 9.0, 5.0 / 18.0],
            [0.5 - 0.1 * sqrt(15), 0.5, 0.5 + 0.1 * sqrt(15)],
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
Dimensions for the specialized IRK stagewise scheme.
    reminder:
    dims.NLP_x : number of state variables
    dims.NLP_u : number of control variables
    time.control_steps : not fully clear yet, assumed to be 1 here
    dims.NLP_v : number of additional optimization variables (free final or initial time for example)

    returns:
    step_variables_block : number of variables stored for a single time step
    state_stage_eqs_block : how many state-vector-sized blocks are stored for the equations, i.e. the state equation and the stages
    dim_NLP_variables : total number of optimization variables in the problem, all steps plus the final step
    dim_NLP_constraints : for each time step, ((state + stage constraints) + path constraints) + path constraints at final time + boundary constraints

"""
function IRK_stagewise_dims(dims::DOCPdims, time::DOCPtime, stage::Int)

    control_block = dims.NLP_u * stage

    # per step:
    # x_i + U_i^1 + ... + U_i^s + K_i^1 + ... + K_i^s
    step_variables_block = dims.NLP_x + control_block + stage * dims.NLP_x

    # state equation + stage equations
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
Convention: 1 <= i <= time.steps, j in {1,...,s}.
"""
function get_stagecontrol_at_time_step(xu, docp::DOCP{<:GenericIRKStagewise}, i::Int, j::Int)
    disc = disc_model(docp)
    dims = docp.dims

    @assert 1 <= j <= disc.stage
    # offset to move to the start of the block
    # X_i, U_i^1, ..., U_i^s, K_i^1, ..., K_i^s,
    var_offset = (i - 1) * disc._step_variables_block
    # we now are here: U_i^1, ..., U_i^s, K_i^1, ..., K_i^s
    x_block_end = var_offset + dims.NLP_x
    # define start and end as a function of the control size
    ctrl0 = x_block_end + (j - 1) * dims.NLP_u + 1
    ctrl1 = x_block_end + j * dims.NLP_u
    # return the view
    return @view xu[ctrl0:ctrl1]
end

"""
Compatibility accessor: return an averaged control at time step.
Note that j is here ignored. See generic control getter in common.jl:90
+++ update: return actual stage controls for solution ?
+++ nb incompatibility between control.steps ans stage for j...
+++ would require a new getter for control block ?
"""
function get_OCP_control_at_time_step(xu, docp::DOCP{<:GenericIRKStagewise}, i::Int, j::Int)
    (i == docp.time.steps +  1) && (i = docp.time.steps) 
    disc = disc_model(docp)
    ui = disc.butcher_b[1] .* get_stagecontrol_at_time_step(xu, docp, i, 1)
    for j=2:disc.stage
        ui += disc.butcher_b[j] .* get_stagecontrol_at_time_step(xu, docp, i, j)
    end
    return ui
end

"""
Return stage variables K_i^j on time step i.
Convention: 1 <= i <= time.steps, j in {1,...,s}.
We extract K_i^j.
"""
function get_stagevars_at_time_step(xu, docp::DOCP{<:GenericIRKStagewise}, i::Int, j::Int)
    disc = disc_model(docp)
    dims = docp.dims

    @assert 1 <= j <= disc.stage

    var_offset = (i - 1) * disc._step_variables_block
    stage_block_start = var_offset + dims.NLP_x + disc._control_block
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
function setWorkArray(docp::DOCP{<:GenericIRKStagewise}, xu, time_grid, v)
    dims = docp.dims
    work = similar(xu, dims.NLP_x + dims.NLP_x)
    return work
end

# ----------------------------------------------------------------------------
# Bounds + initial guess
# ----------------------------------------------------------------------------
"""
Give the bounds for each NLP variable.
var_l <= xu <= var_u
var_l and var_u have the same structure as xu, so we can
reuse the same functions to access the correct location.
"""
function __variables_bounds!(docp::DOCP{<:GenericIRKStagewise})
    # retrieve the bounds and the continuous model
    var_l = docp.bounds.var_l
    var_u = docp.bounds.var_u
    ocp = ocp_model(docp)
    disc = disc_model(docp)
    #
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
    # we change var_l and var_u at every time step with the
    # right bound
    for i in 1:(docp.time.steps + 1)
        set_state_at_time_step!(var_l, x_lb, docp, i)
        set_state_at_time_step!(var_u, x_ub, docp, i)
    end
    # if there is a control variable (change in the new zero-control version)
    if docp.dims.NLP_u > 0
        # for each timestep
        for i in 1:docp.time.steps
            # for each control stage
            for j in 1:disc.stage
                # access the correct index and use the view to set it
                # should later be replaced with a proper setter
                get_stagecontrol_at_time_step(var_l, docp, i, j) .= u_lb
                get_stagecontrol_at_time_step(var_u, docp, i, j) .= u_ub
            end
        end
    end
    # same thing for supplementary variables
    if docp.dims.NLP_v > 0
        v_lb, v_ub = build_bounds_block(
            docp.dims.NLP_v,
            CTModels.dim_variable_constraints_box(ocp),
            CTModels.variable_constraints_box(ocp),
        )
        # stored at the end of xu, therefore also at the end of var_l and var_u
        # so the old setters still work
        set_optim_variable!(var_l, v_lb, docp)
        set_optim_variable!(var_u, v_ub, docp)
    end

    return var_l, var_u
end

function __initial_guess(
    docp::DOCP{<:GenericIRKStagewise},
    init::CTModels.InitialGuess,
)
    # initialize everything to 0.1
    NLP_X = 0.1 * ones(docp.dim_NLP_variables)
    disc = disc_model(docp)
    # if we have initial variable values, place the available values in the correct location
    set_optim_variable!(NLP_X, init.variable, docp)
    # fill the t_i points with the initial state
    time_grid = get_time_grid(NLP_X, docp)
    for i in 1:(docp.time.steps + 1)
        ti = time_grid[i]
        set_state_at_time_step!(NLP_X, init.state(ti), docp, i)
    end
    # if there is a control
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
function integral(docp::DOCP{<:GenericIRKStagewise}, xu, v, time_grid, f)
    disc = disc_model(docp)
    dims = docp.dims

    value = 0.0
    work_xij = similar(xu, dims.NLP_x)
    # for each time step i
    for i in 1:docp.time.steps
        # starting time ti
        ti = time_grid[i]
        # get the state at the beginning xi
        xi = get_OCP_state_at_time_step(xu, docp, i)
        # compute the step size hi
        hi = time_grid[i + 1] - ti

        local_sum = 0.0
        # for each stage
        for j in 1:disc.stage
            # compute the stage time
            tij = ti + disc.butcher_c[j] * hi
            # retrieve the stage control uij
            uij = get_stagecontrol_at_time_step(xu, docp, i, j)
            # reconstruct the stage state xij
            # copy x_i into work_xij
            @views @. work_xij = xi
            # for each stage (implicit)
            for l in 1:disc.stage
                # retrieve K_i^l
                kil = get_stagevars_at_time_step(xu, docp, i, l)
                # update work_xij
                @views @. work_xij = work_xij + hi * disc.butcher_a[j, l] * kil
            end
            # evaluate the cost f(tij,xij,uij,v) and add it weighted by b_j
            local_sum += disc.butcher_b[j] * f(tij, work_xij, uij, v)
        end
        # add hi * local_sum to the total value
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
    docp::DOCP{<:GenericIRKStagewise},
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
    # allocate from the work vector defined in setWorkArray
    work_xij = @view work[1:dims.NLP_x]
    work_sumbk = @view work[(dims.NLP_x + 1):(2 * dims.NLP_x)]
    # compute where the constraints for step i start
    offset = (i - 1) * (disc._state_stage_eqs_block + disc._step_pathcons_block)
    # start, end, and size of the step
    ti = time_grid[i]
    tip1 = time_grid[i + 1]
    hi = tip1 - ti

    # state at the ends of the step
    xi = get_OCP_state_at_time_step(xu, docp, i)
    xip1 = get_OCP_state_at_time_step(xu, docp, i + 1)
    # first dimensions reserved for the state equation
    offset_stage_eqs = dims.NLP_x
    # for the stages
    for j in 1:disc.stage
        tij = ti + disc.butcher_c[j] * hi
        kij = get_stagevars_at_time_step(xu, docp, i, j)
        uij = get_stagecontrol_at_time_step(xu, docp, i, j)
        # build Σ b_j K_i^j
        if j == 1
            @views @. work_sumbk = disc.butcher_b[j] * kij
        else
            @views @. work_sumbk = work_sumbk + disc.butcher_b[j] * kij
        end
        # x_i^j = x_i + h_i Σ_l a_{j,l} K_i^l
        @views @. work_xij = xi
        for l in 1:disc.stage
            kil = get_stagevars_at_time_step(xu, docp, i, l)
            @views @. work_xij = work_xij + hi * disc.butcher_a[j, l] * kil
        end
        # compute f(t_i^j, x_i^j, u_i^j, v)
        CTModels.dynamics(ocp)(
            @view(c[(offset + offset_stage_eqs + 1):(offset + offset_stage_eqs + dims.NLP_x)]),
            tij,
            work_xij,
            uij,
            v,
        )
        # store K_i^j - f(t_i^j, x_i^j, u_i^j, v)
        # reminder: @views @. -> without reallocating, element by element
        @views @. c[(offset + offset_stage_eqs + 1):(offset + offset_stage_eqs + dims.NLP_x)] =
            kij - c[(offset + offset_stage_eqs + 1):(offset + offset_stage_eqs + dims.NLP_x)]
        # offset for the next constraint
        offset_stage_eqs += dims.NLP_x
    end
    # state evolution constraint x_{i+1} - (x_i + h_i Σ_j b_j K_i^j)
    @views @. c[(offset + 1):(offset + dims.NLP_x)] =
        xip1 - (xi + hi * work_sumbk)

    return nothing
end


"""
$(TYPEDSIGNATURES)

Build sparsity pattern for Jacobian of constraints
"""
function DOCP_Jacobian_pattern(docp::DOCP{<: GenericIRKStagewise})
    disc = disc_model(docp)
    dims = docp.dims

    # vector format for sparse matrix
    Is = Vector{Int}(undef, 0)
    Js = Vector{Int}(undef, 0)

    s = disc.stage

    # index alias for v
    v_start = docp.dim_NLP_variables - dims.NLP_v + 1
    v_end = docp.dim_NLP_variables

    # 1. main loop over steps
    for i in 1:docp.time.steps

        # constraints block and offset: state equation, stage equations, path constraints
        c_block = disc._state_stage_eqs_block + disc._step_pathcons_block
        c_offset = (i-1)*c_block
        dyn_start = c_offset + 1
        dyn_end = c_offset + dims.NLP_x
        dyn_lag = c_offset + dims.NLP_x
        stage_start = c_offset + dims.NLP_x + 1
        stage_end = c_offset + (s+1) * dims.NLP_x
        path_start = c_offset + (s+1) * dims.NLP_x + 1
        path_end = c_offset + c_block

        # contiguous variables blocks will be used when possible
        # x_i u_ij k_ij x_i+1 
        var_offset = (i-1)*disc._step_variables_block
        xi_start = var_offset + 1
        xi_end = var_offset + dims.NLP_x
        ui_start = var_offset + dims.NLP_x + 1
        ui_end = var_offset + dims.NLP_x + dims.NLP_u*s
        ki_start = var_offset + dims.NLP_x + dims.NLP_u*s + 1
        ki_end = var_offset + disc._step_variables_block
        xip1_end = var_offset + disc._step_variables_block + dims.NLP_x

        # 1.1 state eq 0 = x_i+1 - (x_i + h_i sum_j b_j k_ij)
        # depends on x_i, k_ij, x_i+1, and v for h_i in variable times case !
        # skip u_ij; should skip k_i[n+1] also but annoying...
        add_nonzero_block!(Is, Js, dyn_start, dyn_end, xi_start, xi_end)
        add_nonzero_block!(Is, Js, dyn_start, dyn_end, ki_start, xip1_end)
        add_nonzero_block!(Is, Js, dyn_start, dyn_end, v_start, v_end)

        # 1.3 stage equations k_ij = f(t_ij, x_ij, u_ij, v)
        # with x_ij = x_i + sum_l a_il k_jl
        # depends on x_i, u_ij, k_ij, and v; (we could distinguish each j...)
        add_nonzero_block!(Is, Js, stage_start, stage_end, xi_start, ki_end)
        add_nonzero_block!(Is, Js, stage_start, stage_end, v_start, v_end)

        # 1.4 path constraint g(t_i, x_i, u_i, v)
        # depends on x_i, u_i (check actual value used), v;
        add_nonzero_block!(Is, Js, path_start, path_end, xi_start, ui_end)
        add_nonzero_block!(Is, Js, path_start, path_end, v_start, v_end)
    end

    # 2. final path constraints (xf, uf, v)
    c_offset = docp.time.steps * (disc._state_stage_eqs_block + disc._step_pathcons_block)
    c_block = disc._step_pathcons_block
    var_offset = docp.time.steps*disc._step_variables_block
    xf_start = var_offset + 1
    xf_end = var_offset + dims.NLP_x
    # NB convention u(tf) = U_N-1
    uf_start = var_offset - disc._step_variables_block + dims.NLP_x + 1
    uf_end = var_offset - disc._step_variables_block + dims.NLP_x + dims.NLP_u*s
    add_nonzero_block!(Is, Js, c_offset+1, c_offset+c_block, xf_start, xf_end)
    add_nonzero_block!(Is, Js, c_offset+1, c_offset+c_block, uf_start, uf_end)
    add_nonzero_block!(Is, Js, c_offset+1, c_offset+c_block, v_start, v_end)

    # 3. boundary constraints (x0, xf, v)
    c_offset =
        docp.time.steps * (disc._state_stage_eqs_block + disc._step_pathcons_block) +
        disc._step_pathcons_block
    c_block = dims.boundary_cons
    x0_start = 1
    x0_end = dims.NLP_x
    add_nonzero_block!(Is, Js, c_offset+1, c_offset+c_block, x0_start, x0_end)
    add_nonzero_block!(Is, Js, c_offset+1, c_offset+c_block, xf_start, xf_end)
    add_nonzero_block!(Is, Js, c_offset+1, c_offset+c_block, v_start, v_end)
    # 3.4 null initial condition for lagrangian cost state l0
    if docp.flags.lagrange
        add_nonzero_block!(Is, Js, docp.dim_NLP_constraints, dims.NLP_x)
    end

    # build and return sparse matrix
    nnzj = length(Is)
    Vs = ones(Bool, nnzj)
    return SparseArrays.sparse(Is, Js, Vs, docp.dim_NLP_constraints, docp.dim_NLP_variables)
end

"""
$(TYPEDSIGNATURES)

Build sparsity pattern for Hessian of Lagrangian
"""
function DOCP_Hessian_pattern(docp::DOCP{<: GenericIRKStagewise})
    disc = disc_model(docp)
    dims = docp.dims

    # NB. need to provide full pattern for coloring, not just upper/lower part
    Is = Vector{Int}(undef, 0)
    Js = Vector{Int}(undef, 0)

    s = disc.stage

    # index alias for v
    v_start = docp.dim_NLP_variables - dims.NLP_v + 1
    v_end = docp.dim_NLP_variables

    # 0. objective
    # 0.1 mayer cost (x0, xf, v) 
    # -> grouped with term 3. for boundary conditions
    # 0.2 lagrange case sum h_i l(ti, xi, ui, v)
    # -> included in stage equations terms see 1.2

    # 1. main loop over steps
    # 1.0 v / v term
    add_nonzero_block!(Is, Js, v_start, v_end, v_start, v_end)

    for i in 1:docp.time.steps

        # contiguous variables blocks will be used when possible
        # x_i u_ij k_ij x_i+1
        var_offset = (i-1)*disc._step_variables_block
        xi_start = var_offset + 1
        xi_end = var_offset + dims.NLP_x
        ui_start = var_offset + dims.NLP_x + 1
        ui_end = var_offset + dims.NLP_x + dims.NLP_u*s
        ki_start = var_offset + dims.NLP_x + dims.NLP_u*s + 1
        ki_end = var_offset + disc._step_variables_block

        # 1.1 state eq 0 = x_i+1 - (x_i + h_i sum_j b_j k_ij)
        # -> 2nd order terms are zero

        # 1.2 stage equations 0 = k_ij - f(t_ij, x_ij, u_ij, v)
        # with x_ij = x_i + sum_l a_il k_jl
        # 2nd order terms depend on x_i, u_ij, k_ij, and v; (we could distinguish each j...)
        add_nonzero_block!(Is, Js, xi_start, ki_end, xi_start, ki_end)
        add_nonzero_block!(Is, Js, xi_start, ki_end, v_start, v_end; sym=true)

        # 1.3 path constraint g(t_i, x_i, u_i, v)
        # -> included in 1.3
    end

    # 2. final path constraints (xf, uf, v) (assume present) +++ done in 1.4 above ?
    var_offset = docp.time.steps * disc._step_variables_block
    xf_start = var_offset + 1
    xf_end = var_offset + dims.NLP_x
    # NB convention u(tf) = U_N-1
    uf_start = var_offset - disc._step_variables_block + dims.NLP_x + 1
    uf_end = var_offset - disc._step_variables_block + dims.NLP_x + dims.NLP_u*s
    add_nonzero_block!(Is, Js, xf_start, xf_end, xf_start, xf_end)
    add_nonzero_block!(Is, Js, uf_start, uf_end, uf_start, uf_end)
    add_nonzero_block!(Is, Js, xf_start, xf_end, uf_start, uf_end; sym=true)
    add_nonzero_block!(Is, Js, xf_start, uf_end, v_start, v_end; sym=true)
    add_nonzero_block!(Is, Js, uf_start, uf_end, v_start, v_end; sym=true)

    # 3. boundary constraints (x0, xf, v) or mayer cost g0(x0, xf, v) (assume present)
    # -> x0 / x0, x0 / v terms included in first loop iteration
    # -> xf / xf, xf / v terms included in 2.
    x0_start = 1
    x0_end = dims.NLP_x
    add_nonzero_block!(Is, Js, x0_start, x0_end, xf_start, xf_end; sym=true)

    # build and return sparse matrix
    nnzj = length(Is)
    Vs = ones(Bool, nnzj)
    return SparseArrays.sparse(Is, Js, Vs, docp.dim_NLP_variables, docp.dim_NLP_variables)
end

#= Functions for implicit midpoint discretization scheme
Internal layout for NLP variables: 
[X_1,U_1, .., X_N,U_N, X_N+1, V]
with the convention u([t_i,t_i+1[) = U_i and u(tf) = U_N
NB. This version is much faster than the one using stage variables
Version without work array gives identical performance but is less readable
=#

struct Midpoint <: Scheme
    info::String
    _step_variables_block::Int
    _state_stage_eqs_block::Int
    _step_pathcons_block::Int
    _final_control::Bool

    # constructor
    function Midpoint(dims::DOCPdims, time::DOCPtime)

        step_variables_block = dims.NLP_x + dims.NLP_u * time.control_steps
        state_stage_eqs_block = dims.NLP_x
        step_pathcons_block = dims.path_cons

        # NLP variables size ([state, controls]_1..N, final state, variable)
        dim_NLP_variables = time.steps * step_variables_block + dims.NLP_x + dims.NLP_v

        # NLP constraints size ([state eq, path]_1..N, final_path, boundary)
        dim_NLP_constraints = time.steps * (dims.NLP_x + step_pathcons_block) + 
        step_pathcons_block + dims.boundary_cons

        disc = new(
            "Implicit Midpoint aka Gauss-Legendre collocation for s=1, 2nd order, symplectic",
            step_variables_block,
            state_stage_eqs_block,
            step_pathcons_block,
            false,
        )

        return disc, dim_NLP_variables, dim_NLP_constraints
    end
end

"""
$(TYPEDSIGNATURES)

Set work array for all dynamics cost evaluations
"""
function setWorkArray(docp::DOCP{Midpoint}, xu, time_grid, v)
    work = similar(xu, docp.dims.NLP_x * docp.time.steps * docp.time.control_steps)
    ocp = ocp_model(docp)

    # loop over time steps
    offset = 0
    for i in 1:docp.time.steps
        ts = 0.5 * (time_grid[i] + time_grid[i + 1])
        xs =
            0.5 * (
                get_OCP_state_at_time_step(xu, docp, i) +
                get_OCP_state_at_time_step(xu, docp, i+1)
            )
        # loop over control steps
        for j in 1:docp.time.control_steps
            uij = get_OCP_control_at_time_step(xu, docp, i; j=j)
            # OCP dynamics
            CTModels.dynamics(ocp)(
            (@view work[(offset + 1):(offset + docp.dims.NLP_x)]), ts, xs, uij, v
            )
            offset += docp.dims.NLP_x
        end
    end

    return work
end

"""
$(TYPEDSIGNATURES)

Compute the running cost
"""
function integral(docp::DOCP{Midpoint}, xu, v, time_grid, f)
    value = 0.0

    # Collocation: 1 control per step
    # Direct shooting: >=1 controls per step

    if docp.time.control_steps == 1
        # loop over time steps
        for i in 1:docp.time.steps
            hi = time_grid[i + 1] - time_grid[i]
            ts = 0.5 * (time_grid[i] + time_grid[i + 1])
            xs =
                0.5 * (
                    get_OCP_state_at_time_step(xu, docp, i) +
                    get_OCP_state_at_time_step(xu, docp, i+1)
                )
            ui = get_OCP_control_at_time_step(xu, docp, i)
            value +=  hi * f(ts, xs, ui, v)
        end
    else
        # loop over time steps
        for i in 1:docp.time.steps
            hi = (time_grid[i + 1] - time_grid[i]) / docp.time.control_steps
            xs = 0.5 * (
                    get_OCP_state_at_time_step(xu, docp, i) +
                    get_OCP_state_at_time_step(xu, docp, i+1)
                )
            # loop over control steps
            for j in 1:docp.time.control_steps
                tij = time_grid[i] + (j - 0.5) * hi
                uij = get_OCP_control_at_time_step(xu, docp, i; j=j)
                value +=  hi * f(tij, xs, uij, v)
            end
        end
    end

    return value
end

"""
$(TYPEDSIGNATURES)

Set the constraints corresponding to the state equation
Convention: 1 <= i <= dim_NLP_steps+1
"""
function stepStateConstraints!(docp::DOCP{Midpoint}, c, xu, v, time_grid, i, work)
    ocp = ocp_model(docp)
    disc = disc_model(docp)

    # compute state variables at next step: midpoint rule
    ti = time_grid[i]
    xi = get_OCP_state_at_time_step(xu, docp, i)
    tip1 = time_grid[i + 1]
    xip1 = get_OCP_state_at_time_step(xu, docp, i+1)
    hi = (tip1 - ti) / docp.time.control_steps
    offset_dyn_i = (i-1) * docp.dims.NLP_x * docp.time.control_steps
    # +++ allocations here ?
    x_next = xi
    for j in 1:docp.time.control_steps
        x_next += hi * work[(offset_dyn_i + 1):(offset_dyn_i + docp.dims.NLP_x)]
        offset_dyn_i += docp.dims.NLP_x
    end

    # set state equation as constraints (equal to 0)
    offset = (i-1)*(disc._state_stage_eqs_block + disc._step_pathcons_block)
    @views @. c[(offset + 1):(offset + docp.dims.NLP_x)] = xip1 - x_next

end

"""
$(TYPEDSIGNATURES)

Build sparsity pattern for Jacobian of constraints
"""
function DOCP_Jacobian_pattern(docp::DOCP{Midpoint})
    disc = disc_model(docp)

    # vector format for sparse matrix
    Is = Vector{Int}(undef, 0)
    Js = Vector{Int}(undef, 0)

    # index alias for v
    v_start = docp.dim_NLP_variables - docp.dims.NLP_v + 1
    v_end = docp.dim_NLP_variables

    # 1. main loop over steps
    for i in 1:docp.time.steps

        # constraints block and offset: state equation, path constraints
        c_block = disc._state_stage_eqs_block + disc._step_pathcons_block
        c_offset = (i-1)*c_block

        # contiguous variables blocks will be used when possible
        # x_i (l_i) u_i x_i+1 (l_i+1)
        var_offset = (i-1)*disc._step_variables_block
        xi_start = var_offset + 1
        xi_end = var_offset + docp.dims.NLP_x
        ui_start = var_offset + docp.dims.NLP_x + 1
        ui_end = var_offset + docp.dims.NLP_x + docp.dims.NLP_u
        xip1_end = var_offset + disc._step_variables_block + docp.dims.NLP_x
        var_end = var_offset + disc._step_variables_block + docp.dims.NLP_x

        # dynamics constraint: depends on x_i, u_i, x_i+1 [, v]
        add_nonzero_block!(
            Is, Js, c_offset+1, c_offset+docp.dims.NLP_x, xi_start, xip1_end
        )

        # path constraint: depends on x_i, u_i [, v]
        add_nonzero_block!(
            Is, Js, c_offset+docp.dims.NLP_x+1, c_offset+c_block, xi_start, ui_end
        )

        # dependency wrt v
        add_nonzero_block!(Is, Js, c_offset+1, c_offset+c_block, v_start, v_end)

    end

    # 2. final path constraints (xf, uf, v)
    c_offset = docp.time.steps * (disc._state_stage_eqs_block + disc._step_pathcons_block)
    c_block = disc._step_pathcons_block
    var_offset = docp.time.steps*disc._step_variables_block
    xf_start = var_offset + 1
    xf_end = var_offset + docp.dims.NLP_x
    uf_start = var_offset-disc._step_variables_block + docp.dims.NLP_x + 1
    uf_end = var_offset-disc._step_variables_block + docp.dims.NLP_x + docp.dims.NLP_u
    add_nonzero_block!(Is, Js, c_offset+1, c_offset+c_block, xf_start, xf_end)
    add_nonzero_block!(Is, Js, c_offset+1, c_offset+c_block, uf_start, uf_end)
    add_nonzero_block!(Is, Js, c_offset+1, c_offset+c_block, v_start, v_end)

    # 3. boundary constraints (x0, xf, v)
    c_offset =
        docp.time.steps * (disc._state_stage_eqs_block + disc._step_pathcons_block) +
        disc._step_pathcons_block
    c_block = docp.dims.boundary_cons
    x0_start = 1
    x0_end = docp.dims.NLP_x
    add_nonzero_block!(Is, Js, c_offset+1, c_offset+c_block, x0_start, x0_end)
    add_nonzero_block!(Is, Js, c_offset+1, c_offset+c_block, xf_start, xf_end)
    add_nonzero_block!(Is, Js, c_offset+1, c_offset+c_block, v_start, v_end)

    # build and return sparse matrix
    nnzj = length(Is)
    Vs = ones(Bool, nnzj)
    return sparse(Is, Js, Vs, docp.dim_NLP_constraints, docp.dim_NLP_variables)
end

"""
$(TYPEDSIGNATURES)

Build sparsity pattern for Hessian of Lagrangian
"""
function DOCP_Hessian_pattern(docp::DOCP{Midpoint})
    disc = disc_model(docp)

    # NB. need to provide full pattern for coloring, not just upper/lower part
    Is = Vector{Int}(undef, 0)
    Js = Vector{Int}(undef, 0)

    # index alias for v
    v_start = docp.dim_NLP_variables - docp.dims.NLP_v + 1
    v_end = docp.dim_NLP_variables

    # 0. objective
    # 0.1 mayer cost (x0, xf, v) 
    # -> grouped with term 3. for boundary conditions
    # +++ 0.2 lagrange cost ?

    # 1. main loop over steps
    # 1.0 v / v term
    add_nonzero_block!(Is, Js, v_start, v_end, v_start, v_end)

    for i in 1:docp.time.steps

        # contiguous variables blocks will be used when possible
        # x_i (l_i) u_i x_i+1 (l_i+1)
        var_offset = (i-1)*disc._step_variables_block
        xi_start = var_offset + 1
        xi_end = var_offset + docp.dims.NLP_x
        xip1_end = var_offset + disc._step_variables_block + docp.dims.NLP_x
        ui_start = var_offset + docp.dims.NLP_x + 1

        # 1.1 state eq 0 = x_i+1 - (x_i + h_i * f(t_s, x_s, u_i, v))
        # with t_s = (t_i + t_i+1)/2    x_s = (x_i + x_i+1)/2
        # 2nd order terms depend on x_i, u_i, x_i+1, and v; 
        add_nonzero_block!(Is, Js, xi_start, xip1_end, xi_start, xip1_end)
        add_nonzero_block!(Is, Js, xi_start, xip1_end, v_start, v_end; sym=true)

        # 1.3 path constraint g(t_i, x_i, u_i, v)
        # -> included in 1.1
    end

    # 2. final path constraints (xf, uf, v)
    # -> included in last loop iteration (with x_i+1 as x_f and u_i as u_f)

    # 3. boundary constraints (x0, xf, v) or mayer cost g0(x0, xf, v) (assume present)
    # -> x0 / x0, x0 / v terms included in first loop iteration
    # -> xf / xf, xf / v terms included in last loop iteration (with x_i+1 as x_f)
    x0_start = 1
    x0_end = docp.dims.NLP_x
    var_offset = docp.time.steps*disc._step_variables_block
    xf_start = var_offset + 1
    xf_end = var_offset + docp.dims.NLP_x
    add_nonzero_block!(Is, Js, x0_start, x0_end, xf_start, xf_end; sym=true)

    # 3.1 null initial condition for lagrangian cost state l0
    # -> 2nd order term is zero

    # build and return sparse matrix
    nnzj = length(Is)
    Vs = ones(Bool, nnzj)
    return sparse(Is, Js, Vs, docp.dim_NLP_variables, docp.dim_NLP_variables)
end

#= Functions for trapeze discretization scheme
Internal layout for NLP variables: 
[X_1,U_1, .., X_N+1,U_N+1, V]
=#

struct Trapeze <: Discretization
    info::String
    _step_variables_block::Int
    _state_stage_eqs_block::Int
    _step_pathcons_block::Int
    _final_control::Bool

    # constructor
    function Trapeze(
        dim_NLP_steps, dim_NLP_x, dim_NLP_u, dim_NLP_v, dim_path_cons, dim_boundary_cons
    )

        # Trapeze is better with final control (used in final dynamics)
        final_control = true #false about 10% slower

        # aux variables
        step_variables_block = dim_NLP_x + dim_NLP_u
        state_stage_eqs_block = dim_NLP_x
        step_pathcons_block = dim_path_cons

        # NLP variables size ([state, control]_1..N+1, variable)
        dim_NLP_variables = dim_NLP_steps * step_variables_block + dim_NLP_x + dim_NLP_v
        final_control && (dim_NLP_variables += dim_NLP_u)

        # NLP constraints size ([dynamics, stage, path]_1..N, final path, boundary, variable)
        dim_NLP_constraints =
            dim_NLP_steps * (state_stage_eqs_block + step_pathcons_block) +
            step_pathcons_block +
            dim_boundary_cons

        disc = new(
            "Implicit Trapeze aka Crank-Nicolson, 2nd order, A-stable",
            step_variables_block,
            state_stage_eqs_block,
            step_pathcons_block,
            final_control,
        )

        return disc, dim_NLP_variables, dim_NLP_constraints
    end
end

"""
$(TYPEDSIGNATURES)

Set work array for all dynamics evaluations
"""
function setWorkArray(docp::DOCP{Trapeze}, xu, time_grid, v)

    #= use work array to store all dynamics
    NB. using a smaller work array to store a single dynamics between steps
    appears slower, maybe due to the copy involved ? =#
    dims = docp.dims
    work = similar(xu, dims.NLP_x * (docp.time.steps+1))

    # loop over time steps
    for i in 1:(docp.time.steps + 1)
        offset = (i-1) * dims.NLP_x
        ti = time_grid[i]
        xi = get_OCP_state_at_time_step(xu, docp, i)
        ui = get_OCP_control_at_time_step(xu, docp, i)
        # OCP dynamics
        CTModels.dynamics(docp.ocp)(
            (@view work[(offset + 1):(offset + dims.NLP_x)]), ti, xi, ui, v
        )
    end
    return work
end

"""
$(TYPEDSIGNATURES)

Compute the running cost
"""
function runningCost(docp::DOCP{Trapeze}, xu, v, time_grid)
    dims = docp.dims
    obj_lagrange = 0.0

    # sum_i=1..N h_i * (l_i + l_i+1) / 2 = (h_1 / 2) l_1 + sum_i=2..N (h_i-1+h_i)/2 * l_i + (h_N / 2) l_N+1 

    # first term
    i = 1
    ti = time_grid[i]
    xi = get_OCP_state_at_time_step(xu, docp, i)
    ui = get_OCP_control_at_time_step(xu, docp, i)
    h = time_grid[i + 1] - time_grid[i]
    obj_lagrange = obj_lagrange + h / 2.0 * CTModels.lagrange(docp.ocp)(ti, xi, ui, v)

    # loop over time steps
    for i in 2:docp.time.steps
        offset = (i-1) * dims.NLP_x
        ti = time_grid[i]
        xi = get_OCP_state_at_time_step(xu, docp, i)
        ui = get_OCP_control_at_time_step(xu, docp, i)
        h2 = time_grid[i + 1] - time_grid[i - 1]
        obj_lagrange = obj_lagrange + h2 / 2.0 * CTModels.lagrange(docp.ocp)(ti, xi, ui, v)
    end

    # last term
    i = docp.time.steps+1
    ti = time_grid[i]
    xi = get_OCP_state_at_time_step(xu, docp, i)
    ui = get_OCP_control_at_time_step(xu, docp, i)
    h = time_grid[i] - time_grid[i - 1]
    obj_lagrange = obj_lagrange + h / 2.0 * CTModels.lagrange(docp.ocp)(ti, xi, ui, v)

    return obj_lagrange
end

"""
$(TYPEDSIGNATURES)

Set the constraints corresponding to the state equation
Convention: 1 <= i <= dim_NLP_steps+1
"""
function setStepConstraints!(docp::DOCP{Trapeze}, c, xu, v, time_grid, i, work)
    disc = disc_model(docp)
    dims = docp.dims

    # offset for previous steps
    offset = (i-1)*(disc._state_stage_eqs_block + disc._step_pathcons_block)
    # c block version 
    # offset = 0

    # 0. variables
    ti = time_grid[i]
    xi = get_OCP_state_at_time_step(xu, docp, i)

    # 1. state equation
    if i <= docp.time.steps
        # more variables
        tip1 = time_grid[i + 1]
        xip1 = get_OCP_state_at_time_step(xu, docp, i+1)
        half_hi = 0.5 * (tip1 - ti)
        offset_dyn_i = (i-1)*dims.NLP_x
        offset_dyn_ip1 = i*dims.NLP_x

        # trapeze rule (no allocations ^^)
        @views @. c[(offset + 1):(offset + dims.NLP_x)] =
            xip1 - (
                xi +
                half_hi * (
                    work[(offset_dyn_i + 1):(offset_dyn_i + dims.NLP_x)] +
                    work[(offset_dyn_ip1 + 1):(offset_dyn_ip1 + dims.NLP_x)]
                )
            )

        offset += dims.NLP_x
    end

    # 2. path constraints
    if dims.path_cons > 0
        ui = get_OCP_control_at_time_step(xu, docp, i)
        CTModels.path_constraints_nl(docp.ocp)[2](
            (@view c[(offset + 1):(offset + dims.path_cons)]), ti, xi, ui, v
        )
    end
end

"""
$(TYPEDSIGNATURES)

Build sparsity pattern for Jacobian of constraints
"""
function DOCP_Jacobian_pattern(docp::DOCP{Trapeze})

    # +++ possible improvements (besides getting actual pattern of OCP functions !)
    # - handle l_i separately ie dont skip them but handle them at the end
    # ie put zeros everywhere then re add the few nonzeros
    # NB. requires a new function remove_nnz_block and remove_nnz
    # - handle variable time ie dependency to v via time step h_i ?
    disc = disc_model(docp)
    dims = docp.dims

    # vector format for sparse matrix
    Is = Vector{Int}(undef, 0)
    Js = Vector{Int}(undef, 0)

    # index alias for v
    v_start = docp.dim_NLP_variables - dims.NLP_v + 1
    v_end = docp.dim_NLP_variables

    # 1. main loop over steps
    for i in 1:docp.time.steps

        # constraints block and offset: state equation, path constraints
        c_block = disc._state_stage_eqs_block + disc._step_pathcons_block
        c_offset = (i-1)*c_block
        dyn_start = c_offset + 1
        dyn_end = c_offset + dims.NLP_x
        dyn_lag = c_offset + dims.NLP_x
        path_start = c_offset + dims.NLP_x + 1
        path_end = c_offset + c_block

        # contiguous variables blocks will be used when possible 
        # x_i (l_i) u_i x_i+1 (l_i+1) u_i+1
        var_block = disc._step_variables_block * 2
        var_offset = (i-1)*disc._step_variables_block
        xi_start = var_offset + 1
        xi_end = var_offset + dims.NLP_x
        ui_start = var_offset + dims.NLP_x + 1
        ui_end = var_offset + dims.NLP_x + dims.NLP_u
        xip1_end = var_offset + dims.NLP_x + dims.NLP_u + dims.NLP_x
        uip1_start = var_offset + dims.NLP_x*2 + dims.NLP_u + 1
        uip1_end = var_offset + dims.NLP_x*2 + dims.NLP_u*2

        # 1.1 state eq 0 = x_i+1 - (x_i + h_i/2 (f(t_i,x_i,u_i,v) + f(t_i+1,x_i+1,u_i+1,v)))
        # depends on x_i, u_i, x_i+1, u_i+1; skip l_i, l_i+1; v cf 1.4
        add_nonzero_block!(Is, Js, dyn_start, dyn_end, xi_start, xi_end)
        add_nonzero_block!(Is, Js, dyn_start, dyn_end, ui_start, xip1_end)
        add_nonzero_block!(Is, Js, dyn_start, dyn_end, uip1_start, uip1_end)

        # 1.3 path constraint g(t_i, x_i, u_i, v) 
        # depends on x_i, u_i; skip l_i; v cf 1.4
        add_nonzero_block!(Is, Js, path_start, path_end, xi_start, xi_end)
        add_nonzero_block!(Is, Js, path_start, path_end, ui_start, ui_end)

        # 1.4 whole constraint block depends on v
        add_nonzero_block!(Is, Js, path_start, path_end, v_start, v_end)
    end

    # 2. final path constraints (xf, uf, v)
    c_offset = docp.time.steps*(disc._state_stage_eqs_block + disc._step_pathcons_block)
    c_block = disc._step_pathcons_block
    var_offset = docp.time.steps*disc._step_variables_block
    xf_start = var_offset + 1
    xf_end = var_offset + dims.NLP_x
    uf_start = var_offset + dims.NLP_x + 1
    uf_end = var_offset + dims.NLP_x + dims.NLP_u
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

    # build and return sparse matrix
    nnzj = length(Is)
    Vs = ones(Bool, nnzj)
    return sparse(Is, Js, Vs, docp.dim_NLP_constraints, docp.dim_NLP_variables)
end

"""
$(TYPEDSIGNATURES)

Build sparsity pattern for Hessian of Lagrangian
"""
function DOCP_Hessian_pattern(docp::DOCP{Trapeze})
    disc = disc_model(docp)
    dims = docp.dims

    # NB. need to provide full pattern for coloring, not just upper/lower part
    Is = Vector{Int}(undef, 0)
    Js = Vector{Int}(undef, 0)

    # index alias for v
    v_start = docp.dim_NLP_variables - dims.NLP_v + 1
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
        # x_i (l_i) u_i x_i+1 (l_i+1) u_i+1
        var_block = disc._step_variables_block * 2
        var_offset = (i-1)*disc._step_variables_block

        # 1.1 state eq 0 = x_i+1 - (x_i + h_i/2 (f(t_i,x_i,u_i,v) + f(t_i+1,x_i+1,u_i+1,v)))
        # 2nd order terms depend on x_i, u_i, x_i+1, u_i+1, and v -> included in 1.2
        # -> use single block for all step variables
        add_nonzero_block!(
            Is, Js, var_offset+1, var_offset+var_block, var_offset+1, var_offset+var_block
        )
        add_nonzero_block!(
            Is, Js, var_offset+1, var_offset+var_block, v_start, v_end; sym=true
        )

        # 1.3 path constraint g(t_i, x_i, u_i, v)
        # -> included in 1.2
    end

    # 2. final path constraints (xf, uf, v)
    # -> included in last loop iteration

    # 3. boundary constraints (x0, xf, v)
    # -> (x0, v) terms included in first loop iteration
    # -> (xf, v) terms included in last loop iteration
    if docp.flags.mayer || dims.boundary_cons > 0
        var_offset = docp.time.steps*disc._step_variables_block
        x0_start = 1
        x0_end = dims.NLP_x
        xf_start = var_offset + 1
        xf_end = var_offset + dims.NLP_x
        add_nonzero_block!(Is, Js, x0_start, x0_end, xf_start, xf_end; sym=true)
    end
    # 3.1 null initial condition for lagrangian cost state l0
    # -> 2nd order term is zero

    # build and return sparse matrix
    nnzj = length(Is)
    Vs = ones(Bool, nnzj)
    return sparse(Is, Js, Vs, docp.dim_NLP_variables, docp.dim_NLP_variables)
end

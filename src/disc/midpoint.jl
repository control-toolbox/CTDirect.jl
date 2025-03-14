#= Functions for implicit midpoint discretization scheme
Internal layout for NLP variables: 
[X_1,U_1, .., X_N,U_N, X_N+1, V]
with the convention u([t_i,t_i+1[) = U_i and u(tf) = U_N
NB. This version is much faster than the one using stage variables
Version without work array gives identical performance but is less readable
=#

struct Midpoint <: Discretization

    info::String
    _step_variables_block::Int
    _state_stage_eqs_block::Int
    _step_pathcons_block::Int

    # constructor
    function Midpoint(dim_NLP_steps, dim_NLP_x, dim_NLP_u, dim_NLP_v, dim_path_cons, dim_boundary_cons)

        # aux variables
        step_variables_block = dim_NLP_x + dim_NLP_u
        state_stage_eqs_block = dim_NLP_x
        step_pathcons_block = dim_path_cons

        # NLP variables size ([state, control, stage]_1..N, final state, variable)
        dim_NLP_variables = dim_NLP_steps * step_variables_block + dim_NLP_x + dim_NLP_v
        
        # NLP constraints size ([dynamics, path]_1..N, final path, boundary, variable)
        dim_NLP_constraints = dim_NLP_steps * (state_stage_eqs_block + step_pathcons_block) + step_pathcons_block + dim_boundary_cons

        disc = new("Implicit Midpoint aka Gauss-Legendre collocation for s=1, 2nd order, symplectic", step_variables_block, state_stage_eqs_block, step_pathcons_block)

        return disc, dim_NLP_variables, dim_NLP_constraints
    end
end


"""
$(TYPEDSIGNATURES)

Set work array for all dynamics and lagrange cost evaluations
"""
function setWorkArray(docp::DOCP{Midpoint}, xu, time_grid, v)
    
    work = similar(xu, docp.dim_NLP_x * docp.dim_NLP_steps)

    # loop over time steps
    for i = 1:docp.dim_NLP_steps
        offset = (i-1) * docp.dim_NLP_x
        ts = 0.5 * (time_grid[i] + time_grid[i+1])
        xs = 0.5 * (get_OCP_state_at_time_step(xu, docp, i) + get_OCP_state_at_time_step(xu, docp, i+1))
        ui = get_OCP_control_at_time_step(xu, docp, i)
        # OCP dynamics
        CTModels.dynamics(docp.ocp)((@view work[offset+1:offset+docp.dim_OCP_x]), ts, xs, ui, v)
        # lagrange cost
        if docp.has_lagrange
            work[offset+docp.dim_NLP_x] = CTModels.lagrange(docp.ocp)(ts, xs, ui, v)
        end   
    end

    return work
end


"""
$(TYPEDSIGNATURES)

Set the constraints corresponding to the state equation
Convention: 1 <= i <= dim_NLP_steps+1
"""
function setStepConstraints!(docp::DOCP{Midpoint}, c, xu, v, time_grid, i, work)
    
    # offset for previous steps
    offset = (i-1)*(docp.discretization._state_stage_eqs_block + docp.discretization._step_pathcons_block)

    # 0. variables
    ti = time_grid[i]
    xi = get_OCP_state_at_time_step(xu, docp, i)

    # 1. state equation
    if i <= docp.dim_NLP_steps
        # more variables
        tip1 = time_grid[i+1]
        xip1 = get_OCP_state_at_time_step(xu, docp, i+1)
        hi = tip1 - ti
        offset_dyn_i = (i-1)*docp.dim_NLP_x
       
        # state equation: midpoint rule
        @views @. c[offset+1:offset+docp.dim_OCP_x] = xip1 - (xi + hi * work[offset_dyn_i+1:offset_dyn_i+docp.dim_OCP_x])
        if docp.has_lagrange
            c[offset+docp.dim_NLP_x] = get_lagrange_state_at_time_step(xu, docp, i+1) - (get_lagrange_state_at_time_step(xu, docp, i) + hi * work[offset_dyn_i+docp.dim_NLP_x])
        end
        offset += docp.dim_NLP_x

    end
   
    # 2. path constraints
    if docp.dim_path_cons > 0
        ui = get_OCP_control_at_time_step(xu, docp, i)
        CTModels.path_constraints_nl(docp.ocp)[2]((@view c[offset+1:offset+docp.dim_path_cons]), ti, xi, ui, v)
        offset += docp.dim_path_cons
    end

end


"""
$(TYPEDSIGNATURES)

Build sparsity pattern for Jacobian of constraints
"""
function DOCP_Jacobian_pattern(docp::DOCP{Midpoint})

    # vector format for sparse matrix
    Is = Vector{Int}(undef, 0)
    Js = Vector{Int}(undef, 0)

    # index alias for v
    v_start = docp.dim_NLP_variables - docp.dim_NLP_v + 1
    v_end = docp.dim_NLP_variables

    # 1. main loop over steps
    for i = 1:docp.dim_NLP_steps

        # constraints block and offset: state equation, stage equation, path constraints
        c_block = docp.discretization._state_stage_eqs_block + docp.discretization._step_pathcons_block
        c_offset = (i-1)*c_block

        # contiguous variables blocks will be used when possible
        # x_i (l_i) u_i x_i+1 (l_i+1)
        var_offset = (i-1)*docp.discretization._step_variables_block
        xi_start = var_offset + 1
        xi_end = var_offset + docp.dim_OCP_x
        ui_start = var_offset + docp.dim_NLP_x + 1
        ui_end = var_offset + docp.dim_NLP_x + docp.dim_NLP_u
        xip1_end = var_offset + docp.discretization._step_variables_block + docp.dim_OCP_x
        var_end = var_offset + docp.discretization._step_variables_block + docp.dim_NLP_x

        # 1.1 state eq 0 = x_i+1 - (x_i + h_i * f(t_s, x_s, u_i, v))
        # with t_s = (t_i + t_i+1) / 2 and x_s = (x_i + x_i+1) / 2
        # depends on x_i, x_i+1, u_i, and v for h_i, t_s in variable times case
        # skip l_i
        add_nonzero_block!(Is, Js, c_offset+1, c_offset+docp.dim_OCP_x, xi_start, xi_end)
        add_nonzero_block!(Is, Js, c_offset+1, c_offset+docp.dim_OCP_x, ui_start, xip1_end)
        add_nonzero_block!(Is, Js, c_offset+1, c_offset+docp.dim_OCP_x, v_start, v_end)
        # 1.2 lagrange part 0 = l_i+1 - (l_i + h_i * l(t_s, x_s, u_i, v))
        # depends on l_i, l_i+1, x_i, x_i+1, u_i, and v for h_i, t_s in variable times case
        if docp.has_lagrange
            add_nonzero_block!(Is, Js, c_offset+docp.dim_NLP_x, c_offset+docp.dim_NLP_x, xi_start, var_end)
            add_nonzero_block!(Is, Js, c_offset+docp.dim_NLP_x, c_offset+docp.dim_NLP_x, v_start, v_end)
        end

        # 1.3 path constraint g(t_i, x_i, u_i, v)
        # depends on x_i, u_i, v; skip l_i
        add_nonzero_block!(Is, Js, c_offset+docp.dim_NLP_x+1, c_offset+c_block, xi_start, xi_end)
        add_nonzero_block!(Is, Js, c_offset+docp.dim_NLP_x+1, c_offset+c_block, ui_start, ui_end)
        add_nonzero_block!(Is, Js, c_offset+docp.dim_NLP_x+1, c_offset+c_block, v_start, v_end)
    end

    # 2. final path constraints (xf, uf, v)
    c_offset = docp.dim_NLP_steps * (docp.discretization._state_stage_eqs_block + docp.discretization._step_pathcons_block)
    c_block = docp.discretization._step_pathcons_block
    var_offset = docp.dim_NLP_steps*docp.discretization._step_variables_block
    xf_start = var_offset + 1
    xf_end = var_offset + docp.dim_OCP_x
    uf_start = var_offset-docp.discretization._step_variables_block + docp.dim_NLP_x + 1
    uf_end = var_offset-docp.discretization._step_variables_block + docp.dim_NLP_x + docp.dim_NLP_u
    add_nonzero_block!(Is, Js, c_offset+1, c_offset+c_block, xf_start, xf_end)
    add_nonzero_block!(Is, Js, c_offset+1,c_offset+c_block, uf_start, uf_end)
    add_nonzero_block!(Is, Js, c_offset+1, c_offset+c_block, v_start, v_end)

    # 3. boundary constraints (x0, xf, v)
    c_offset = docp.dim_NLP_steps * (docp.discretization._state_stage_eqs_block + docp.discretization._step_pathcons_block) + docp.discretization._step_pathcons_block
    c_block = docp.dim_boundary_cons
    x0_start = 1
    x0_end = docp.dim_OCP_x
    add_nonzero_block!(Is, Js, c_offset+1, c_offset+c_block, x0_start, x0_end)
    add_nonzero_block!(Is, Js, c_offset+1, c_offset+c_block, xf_start, xf_end)
    add_nonzero_block!(Is, Js, c_offset+1, c_offset+c_block, v_start, v_end)
    # 3.4 null initial condition for lagrangian cost state l0
    if docp.has_lagrange
        add_nonzero_block!(Is, Js, docp.dim_NLP_constraints, docp.dim_NLP_x)
    end

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

    # NB. need to provide full pattern for coloring, not just upper/lower part
    Is = Vector{Int}(undef, 0)
    Js = Vector{Int}(undef, 0)

    # index alias for v
    v_start = docp.dim_NLP_variables - docp.dim_NLP_v + 1
    v_end = docp.dim_NLP_variables

    # 0. objective
    # 0.1 mayer cost (x0, xf, v) 
    # -> grouped with term 3. for boundary conditions
    # 0.2 lagrange case (lf)
    # -> 2nd order term is zero
   
    # 1. main loop over steps
    # 1.0 v / v term
    add_nonzero_block!(Is, Js, v_start, v_end, v_start, v_end)

    for i = 1:docp.dim_NLP_steps

        # contiguous variables blocks will be used when possible
        # x_i (l_i) u_i x_i+1 (l_i+1)
        var_offset = (i-1)*docp.discretization._step_variables_block
        xi_start = var_offset + 1
        xi_end = var_offset + docp.dim_OCP_x
        xip1_end = var_offset + docp.discretization._step_variables_block + docp.dim_NLP_x
        ui_start = var_offset + docp.dim_NLP_x + 1

        # 1.1 state eq 0 = x_i+1 - (x_i + h_i * f(t_s, x_s, u_i, v))
        # with t_s = (t_i + t_i+1)/2    x_s = (x_i + x_i+1)/2
        # depends on x_i, u_i, x_i+1, and v; skip l_i
        add_nonzero_block!(Is, Js, xi_start, xi_end, xi_start, xi_end)
        add_nonzero_block!(Is, Js, ui_start, xip1_end, ui_start, xip1_end)
        add_nonzero_block!(Is, Js, xi_start, xi_end, ui_start, xip1_end; sym=true)
        add_nonzero_block!(Is, Js, xi_start, xi_end, v_start, v_end; sym=true)
        add_nonzero_block!(Is, Js, ui_start, xip1_end, v_start, v_end; sym=true)

        # 1.2 lagrange part 0 = l_i+1 - (l_i + h_i * l(t_s, x_s, u_i, v))
        # -> included in 1.1 since l_i and l_i+1 have no second order term

        # 1.3 path constraint g(t_i, x_i, u_i, v)
        # -> included in 1.1
    end

    # 2. final path constraints (xf, uf, v)
    # -> included in last loop iteration (with x_i+1 as x_j and u_i as u_f)

    # 3. boundary constraints (x0, xf, v) or mayer cost g0(x0, xf, v) (assume present)
    # -> x0 / x0, x0 / v terms included in first loop iteration
    # -> xf / xf, xf / v terms included in last loop iteration (with x_i+1 as x_f)
    x0_start = 1
    x0_end = docp.dim_OCP_x
    var_offset = docp.dim_NLP_steps*docp.discretization._step_variables_block
    xf_start = var_offset + 1
    xf_end = var_offset + docp.dim_OCP_x
    add_nonzero_block!(Is, Js, x0_start, x0_end, xf_start, xf_end; sym=true)

    # 3.1 null initial condition for lagrangian cost state l0
    # -> 2nd order term is zero
   
    # build and return sparse matrix
    nnzj = length(Is)
    Vs = ones(Bool, nnzj)
    return sparse(Is, Js, Vs, docp.dim_NLP_variables, docp.dim_NLP_variables)

end

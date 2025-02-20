#= Functions for generic implicit Runge Kutta discretization
Internal layout for NLP variables: 
[X_0, U_0, K_0^1..K_0^s, 
 X_1, U_1, K_1^1..K_1^s,
 .., 
 X_N-1, U_N-1, K_N-1^1..K_N-1^s,
 X_N, U_N, V]
with s the stage number and U piecewise constant equal to U_i in [t_i, t_i+1]
or, for methods with s>1, piecewise linear if option control_type set to :linear
NB. U_N may be removed at some point if we disable piecewise linear control
Path constraints are all evaluated at time steps, including final time.
=#


abstract type GenericIRK <: Discretization end

"""
$(TYPEDSIGNATURES)

Implicit Midpoint discretization, formulated as a generic IRK (ie Gauss Legendre 1)
NB. does not use the simplification xs = 0.5 * (xi + xip1) as in midpoint.jl
"""
struct Gauss_Legendre_1 <: GenericIRK

    info::String
    stage::Int
    butcher_a::Matrix{Float64}
    butcher_b::Vector{Float64}
    butcher_c::Vector{Float64}
    _step_variables_block::Int
    _state_stage_eqs_block::Int
    _step_pathcons_block::Int

    function Gauss_Legendre_1(dim_NLP_steps, dim_NLP_x, dim_NLP_u, dim_NLP_v, dim_path_cons, dim_boundary_cons)
        
        stage = 1

        step_variables_block, state_stage_eqs_block, step_pathcons_block, dim_NLP_variables, dim_NLP_constraints = IRK_dims(dim_NLP_steps, dim_NLP_x, dim_NLP_u, dim_NLP_v, dim_path_cons, dim_boundary_cons, stage)

        disc = new("Implicit Midpoint aka Gauss-Legendre collocation for s=1, 2nd order, symplectic, A-stable", stage, hcat(0.5), [1], [0.5], step_variables_block, state_stage_eqs_block, step_pathcons_block)

        return disc, dim_NLP_variables, dim_NLP_constraints
    end
end

"""
$(TYPEDSIGNATURES)

Gauss Legendre 2 discretization, formulated as a generic IRK
"""
struct Gauss_Legendre_2 <: GenericIRK

    info::String
    stage::Int
    butcher_a::Matrix{Float64}
    butcher_b::Vector{Float64}
    butcher_c::Vector{Float64}
    _step_variables_block::Int
    _state_stage_eqs_block::Int
    _step_pathcons_block::Int
    _control_type::Symbol

    function Gauss_Legendre_2(dim_NLP_steps, dim_NLP_x, dim_NLP_u, dim_NLP_v, dim_path_cons, dim_boundary_cons, control_type)
        
        stage = 2

        step_variables_block, state_stage_eqs_block, step_pathcons_block, dim_NLP_variables, dim_NLP_constraints =  IRK_dims(dim_NLP_steps, dim_NLP_x, dim_NLP_u, dim_NLP_v, dim_path_cons, dim_boundary_cons, stage)

        disc = new("Implicit Gauss-Legendre collocation for s=2, 4th order, symplectic, A-stable",stage,
        [0.25 (0.25-sqrt(3) / 6); (0.25+sqrt(3) / 6) 0.25],
        [0.5, 0.5],
        [(0.5 - sqrt(3) / 6), (0.5 + sqrt(3) / 6)],
        step_variables_block, state_stage_eqs_block, step_pathcons_block,
        control_type
        )

        return disc, dim_NLP_variables, dim_NLP_constraints
    end
end


"""
$(TYPEDSIGNATURES)

Gauss Legendre 3 discretization, formulated as a generic IRK
"""
struct Gauss_Legendre_3 <: GenericIRK

    info::String
    stage::Int
    butcher_a::Matrix{Float64}
    butcher_b::Vector{Float64}
    butcher_c::Vector{Float64}
    _step_variables_block::Int
    _state_stage_eqs_block::Int
    _step_pathcons_block::Int
    _control_type::Symbol

    function Gauss_Legendre_3(dim_NLP_steps, dim_NLP_x, dim_NLP_u, dim_NLP_v, dim_path_cons,  dim_boundary_cons, control_type)
        
        stage = 3

        step_variables_block, state_stage_eqs_block, step_pathcons_block, dim_NLP_variables, dim_NLP_constraints =  IRK_dims(dim_NLP_steps, dim_NLP_x, dim_NLP_u, dim_NLP_v, dim_path_cons, dim_boundary_cons, stage)

        disc = new("Implicit Gauss-Legendre collocation for s=3, 6th order, symplectic, A-stable",stage,
        [(5.0/36.0) (2/9 - sqrt(15) / 15) (5/36 - sqrt(15) / 30); 
        (5.0/36.0 + sqrt(15) / 24) (2.0/9.0) (5.0/36.0 - sqrt(15) / 24); 
        (5/36 + sqrt(15) / 30) (2/9 + sqrt(15) / 15) (5.0/36.0)],
        [5.0/18.0, 4.0/9.0, 5.0/18.0],
        [0.5 - 0.1*sqrt(15), 0.5, 0.5 + 0.1*sqrt(15)],
        step_variables_block, state_stage_eqs_block, step_pathcons_block, control_type
        )

        return disc, dim_NLP_variables, dim_NLP_constraints
    end
end


"""
$(TYPEDSIGNATURES)

Return the dimension of the NLP variables and constraints for a generic IRK discretizion, with the control taken constant per step (ie not distinct controls at time stages)
"""
function IRK_dims(dim_NLP_steps, dim_NLP_x, dim_NLP_u, dim_NLP_v, dim_path_cons, dim_boundary_cons, stage)

    # size of variables block for one step: x, u, k
    step_variables_block = dim_NLP_x + dim_NLP_u + dim_NLP_x * stage

    # size of state + stage equations for one step
    state_stage_eqs_block = dim_NLP_x * (1 + stage)

    # size of path constraints block for one step: u, x, xu 
    step_pathcons_block = dim_path_cons

    # NLP variables size ([state, control, stage]_1..N, final state and control, variable)
    dim_NLP_variables = dim_NLP_steps * step_variables_block + dim_NLP_x + dim_NLP_u + dim_NLP_v

    # NLP constraints size ([dynamics, stage, path]_1..N, final path, boundary, variable)
    dim_NLP_constraints = dim_NLP_steps * (state_stage_eqs_block + step_pathcons_block) + step_pathcons_block + dim_boundary_cons

    return step_variables_block, state_stage_eqs_block, step_pathcons_block, dim_NLP_variables, dim_NLP_constraints
end


"""
$(TYPEDSIGNATURES)

Retrieve state variables at given time step from the NLP variables.
Convention: 1 <= i <= dim_NLP_steps+1
Vector output
"""
function get_OCP_state_at_time_step(xu, docp::DOCP{ <: GenericIRK}, i)
    offset = (i-1) * docp.discretization._step_variables_block
    return @view xu[(offset + 1):(offset + docp.dim_OCP_x)]
end
"""
$(TYPEDSIGNATURES)

Retrieve state variable for lagrange cost at given time step from the NLP variables.
Convention: 1 <= i <= dim_NLP_steps+1   (no check for actual lagrange cost presence !)
"""
function get_lagrange_state_at_time_step(xu, docp::DOCP{ <: GenericIRK}, i)
    offset = (i-1) * docp.discretization._step_variables_block
    return xu[offset + docp.dim_NLP_x]
end
"""
$(TYPEDSIGNATURES)

Retrieve control variables at given time step (/stage) from the NLP variables.
Convention: 1 <= i <= dim_NLP_steps
Vector output
Step / Stage versions
"""
function get_OCP_control_at_time_step(xu, docp::DOCP{ <: GenericIRK}, i)
    offset = (i-1) * docp.discretization._step_variables_block + docp.dim_NLP_x
    return @view xu[(offset + 1):(offset + docp.dim_NLP_u)]
end
function get_OCP_control_at_time_stage(xu, docp::DOCP{ <: GenericIRK}, i, cj)
    if (docp.discretization.stage == 1) || (docp.discretization._constant_control)
        # constant interpolation on step
        return get_OCP_control_at_time_step(xu, docp, i)
    else
        # linear interpolation on step
        ui = get_OCP_control_at_time_step(xu, docp, i)
        uip1 = get_OCP_control_at_time_step(xu, docp, i+1)
        return @views @. (1 - cj) * ui + cj * uip1
    end
end


"""
$(TYPEDSIGNATURES)

Retrieve stage variables at given time step/stage from the NLP variables.
Convention: 1 <= i <= dim_NLP_steps,	1 <= j <= s
Vector output
"""
function get_stagevars_at_time_step(xu, docp::DOCP{ <: GenericIRK}, i, j)
    offset = (i-1) * docp.discretization._step_variables_block + docp.dim_NLP_x + docp.dim_NLP_u + (j-1)*docp.dim_NLP_x
    return @view xu[(offset + 1):(offset + docp.dim_NLP_x)]
end

"""
$(TYPEDSIGNATURES)

Set initial guess for state variables at given time step
Convention: 1 <= i <= dim_NLP_steps+1
"""
function set_state_at_time_step!(xu, x_init, docp::DOCP{ <: GenericIRK}, i)
    # initialize only actual state variables from OCP (not lagrange state)
    if !isnothing(x_init)
        offset = (i-1) * docp.discretization._step_variables_block
        xu[(offset + 1):(offset + docp.dim_OCP_x)] .= x_init
    end
end
"""
$(TYPEDSIGNATURES)

Set initial guess for control variables at given time step (/stage)
Convention: 1 <= i <= dim_NLP_steps+1
Step / stage versions
"""
function set_control_at_time_step!(xu, u_init, docp::DOCP{ <: GenericIRK}, i)
    if !isnothing(u_init)
        offset = (i-1) * docp.discretization._step_variables_block + docp.dim_NLP_x
        xu[(offset + 1):(offset + docp.dim_NLP_u)] .= u_init
    end
end

"""
$(TYPEDSIGNATURES)

Set work array for all dynamics and lagrange cost evaluations
"""
function setWorkArray(docp::DOCP{ <: GenericIRK}, xu, time_grid, v)
    # work array layout: [x_ij ; sum_bk]   
    work = similar(xu, docp.dim_OCP_x + docp.dim_NLP_x)
    return work
end

"""
$(TYPEDSIGNATURES)

Set the constraints corresponding to the state equation
Convention: 1 <= i <= dim_NLP_steps (+1)
"""
function setStepConstraints!(docp::DOCP{ <: GenericIRK}, c, xu, v, time_grid, i, work)

    # work array layout: [x_ij ; sum_bk]
    work_xij = @view work[1:docp.dim_OCP_x]
    work_sumbk = @view work[docp.dim_OCP_x+1:docp.dim_OCP_x+docp.dim_NLP_x]

    # offset for previous steps
    offset = (i-1)*(docp.discretization._state_stage_eqs_block + docp.discretization._step_pathcons_block)

    # 0. variables
    ti = time_grid[i]
    xi = get_OCP_state_at_time_step(xu, docp, i)

    # 1. state and stage equations
    if i <= docp.dim_NLP_steps
        # more variables
        tip1 = time_grid[i+1]
        xip1 = get_OCP_state_at_time_step(xu, docp, i+1)        
        hi = tip1 - ti
        offset_stage_eqs = docp.dim_NLP_x

        # loop over stages
        for j=1:docp.discretization.stage

            # time at stage: t_i^j = t_i + c[j] h_i
            cj = docp.discretization.butcher_c[j]
            tij = ti + cj * hi
            
            # stage variables
            kij = get_stagevars_at_time_step(xu, docp, i, j)
            
            # update sum b_j k_i^j (w/ lagrange term) for state equation after loop
            if j == 1
                @views @. work_sumbk[1:docp.dim_NLP_x] = docp.discretization.butcher_b[j] * kij[1:docp.dim_NLP_x]
            else
                @views @. work_sumbk[1:docp.dim_NLP_x] = work_sumbk[1:docp.dim_NLP_x] + docp.discretization.butcher_b[j] * kij[1:docp.dim_NLP_x]
            end

            # state at stage: x_i^j = x_i + h_i sum a_jl k_i^l
            # +++ still some allocations here
            @views @. work_xij[1:docp.dim_OCP_x] = xi
            for l = 1:docp.discretization.stage
                kil = get_stagevars_at_time_step(xu, docp, i, l)
                @views @. work_xij[1:docp.dim_OCP_x] = work_xij[1:docp.dim_OCP_x] + hi * docp.discretization.butcher_a[j,l] * kil[1:docp.dim_OCP_x]
            end

            # control at stage
            uij = get_OCP_control_at_time_stage(xu, docp, i, cj)

            # stage equations k_i^j = f(t_i^j, x_i^j, u_i, v) as c[] = k - f
            # NB. we skip the state equation here, which will be set below
            CTModels.dynamics(docp.ocp)((@view c[offset+offset_stage_eqs+1:offset+offset_stage_eqs+docp.dim_OCP_x]), tij, work_xij, uij, v)
            if docp.has_lagrange
                CTModels.lagrange(docp.ocp)((@view c[offset+offset_stage_eqs+docp.dim_NLP_x:offset+offset_stage_eqs+docp.dim_NLP_x]), tij, work_xij, uij, v)
            end
            @views @. c[offset+offset_stage_eqs+1:offset+offset_stage_eqs+docp.dim_NLP_x] = kij - c[offset+offset_stage_eqs+1:offset+offset_stage_eqs+docp.dim_NLP_x]
            offset_stage_eqs += docp.dim_NLP_x

        end

        # state equation x_i+1 = x_i + h_i sum b_j k_i^j
        @views @. c[offset+1:offset+docp.dim_OCP_x] = xip1 - (xi + hi * work_sumbk[1:docp.dim_OCP_x])
        if docp.has_lagrange
            c[offset+docp.dim_NLP_x] = get_lagrange_state_at_time_step(xu, docp, i+1) - (get_lagrange_state_at_time_step(xu, docp, i) + hi * work_sumbk[docp.dim_NLP_x])
        end

        # update offset for stage and state equations
        offset += docp.dim_NLP_x * (1 + docp.discretization.stage)
    end

    #2. path constraints
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
function DOCP_Jacobian_pattern(docp::DOCP{ <: GenericIRK})

    if docp.discretization._control_type != :constant
        error("Manual Jacobian sparsity pattern not supported for IRK scheme with piecewise linear control")
    end

    # vector format for sparse matrix
    Is = Vector{Int}(undef, 0)
    Js = Vector{Int}(undef, 0)

    s = docp.discretization.stage

    # index alias for v
    v_start = docp.dim_NLP_variables - docp.dim_NLP_v + 1
    v_end = docp.dim_NLP_variables

    # 1. main loop over steps
    for i = 1:docp.dim_NLP_steps

        # constraints block and offset: state equation, stage equations, path constraints
        c_block = docp.discretization._state_stage_eqs_block + docp.discretization._step_pathcons_block
        c_offset = (i-1)*c_block

        # contiguous variables blocks will be used when possible
        # x_i (l_i) u_i k_i x_i+1 (l_i+1)
        var_offset = (i-1)*docp.discretization._step_variables_block
        xi_start = var_offset + 1
        xi_end = var_offset + docp.dim_OCP_x
        ui_start = var_offset + docp.dim_NLP_x + 1
        ui_end = var_offset + docp.dim_NLP_x + docp.dim_NLP_u
        ki_start = var_offset + docp.dim_NLP_x + docp.dim_NLP_u + 1
        ki_end = var_offset + docp.discretization._step_variables_block
        xip1_end = var_offset + docp.discretization._step_variables_block + docp.dim_OCP_x
        li = var_offset + docp.dim_NLP_x
        lip1 = var_offset + docp.discretization._step_variables_block + docp.dim_NLP_x

        # 1.1 state eq 0 = x_i+1 - (x_i + h_i sum_j b_j k_ij)
        # depends on x_i, k_ij, x_i+1, and v for h_i in variable times case !
        # skip l_i, u_i; should skip k_i[n+1] also but annoying...
        add_nonzero_block!(Is, Js, c_offset+1, c_offset+docp.dim_OCP_x, xi_start, xi_end)
        add_nonzero_block!(Is, Js, c_offset+1, c_offset+docp.dim_OCP_x, ki_start, xip1_end)
        add_nonzero_block!(Is, Js, c_offset+1, c_offset+docp.dim_OCP_x, v_start, v_end)
        # 1.2 lagrange part l_i+1 = l_i + h_i (sum_j b_j k_ij)[n+1]
        # depends on l_i, k_ij[n+1], l_i+1, and v for h_i in variable times case !
        if docp.has_lagrange
            add_nonzero_block!(Is, Js, c_offset+docp.dim_NLP_x, li)
            add_nonzero_block!(Is, Js, c_offset+docp.dim_NLP_x, lip1)
            for i=1:s
                kij_l = var_offset + docp.dim_NLP_x + docp.dim_NLP_u + i*docp.dim_NLP_x 
                add_nonzero_block!(Is, Js, c_offset+docp.dim_NLP_x, kij_l)
            end
            add_nonzero_block!(Is, Js, c_offset+docp.dim_NLP_x, c_offset+docp.dim_NLP_x, v_start, v_end)
        end

        # 1.3 stage equations k_ij = f(t_ij, x_ij, u_ij, v) (with lagrange part)
        # with x_ij = x_i + sum_l a_il k_jl and assuming u_ij = u_i
        # depends on x_i, u_i, k_i, and v; skip l_i (could skip k_ij[n+1] too...)
        add_nonzero_block!(Is, Js, c_offset+docp.dim_NLP_x+1, c_offset+(s+1)*docp.dim_NLP_x, xi_start, xi_end)
        add_nonzero_block!(Is, Js, c_offset+docp.dim_NLP_x+1, c_offset+(s+1)*docp.dim_NLP_x, ui_start, ki_end)      
        add_nonzero_block!(Is, Js, c_offset+docp.dim_NLP_x+1, c_offset+(s+1)*docp.dim_NLP_x, v_start, v_end)

        # 1.4 path constraint g(t_i, x_i, u_i, v)
        # depends on x_i, u_i, v; skip l_i
        add_nonzero_block!(Is, Js, c_offset+(s+1)*docp.dim_NLP_x+1, c_offset+c_block, xi_start, xi_end)
        add_nonzero_block!(Is, Js, c_offset+(s+1)*docp.dim_NLP_x+1, c_offset+c_block, ui_start, ui_end)
        add_nonzero_block!(Is, Js, c_offset+(s+1)*docp.dim_NLP_x+1, c_offset+c_block, v_start, v_end)
    end

    # 2. final path constraints (xf, uf, v)
    c_offset = docp.dim_NLP_steps * (docp.discretization._state_stage_eqs_block + docp.discretization._step_pathcons_block)
    c_block = docp.discretization._step_pathcons_block
    var_offset = docp.dim_NLP_steps*docp.discretization._step_variables_block
    xf_start = var_offset + 1
    xf_end = var_offset + docp.dim_OCP_x
    uf_start = var_offset + docp.dim_NLP_x + 1
    uf_end = var_offset + docp.dim_NLP_x + docp.dim_NLP_u
    add_nonzero_block!(Is, Js, c_offset+1, c_offset+c_block, xf_start, xf_end)
    add_nonzero_block!(Is, Js, c_offset+1,c_offset+c_block, uf_start, uf_end)
    add_nonzero_block!(Is, Js, c_offset+1, c_offset+c_block, v_start, v_end)

    # 3. boundary constraints (x0, xf, v)
    c_offset = docp.dim_NLP_steps * (docp.discretization._state_stage_eqs_block + docp.discretization._step_pathcons_block) + docp.discretization._step_pathcons_block
    c_block = docp.dim_boundary_cons + docp.dim_v_cons
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
function DOCP_Hessian_pattern(docp::DOCP{ <: GenericIRK})

    if docp.discretization._control_type != :constant
        error("Manual Hessian sparsity pattern not supported for IRK scheme with piecewise linear control")
    end

    # NB. need to provide full pattern for coloring, not just upper/lower part
    Is = Vector{Int}(undef, 0)
    Js = Vector{Int}(undef, 0)
   
    s = docp.discretization.stage

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
        # x_i (l_i) u_i k_i x_i+1 (l_i+1)
        var_offset = (i-1)*docp.discretization._step_variables_block
        xi_start = var_offset + 1
        xi_end = var_offset + docp.dim_OCP_x
        ui_start = var_offset + docp.dim_NLP_x + 1
        ui_end = var_offset + docp.dim_NLP_x + docp.dim_NLP_u
        ki_start = var_offset + docp.dim_NLP_x + docp.dim_NLP_u + 1
        ki_end = var_offset + (s+1)*docp.dim_NLP_x + docp.dim_NLP_u

        # 1.1 state eq 0 = x_i+1 - (x_i + h_i sum_j b_j k_ij)
        # -> 2nd order terms are zero
        # 1.2 lagrange part 0 = l_i+1 - (l_i + h_i (sum_j b_j k_ij[n+1]))
        # -> 2nd order terms are zero

        # 1.3 stage equations 0 = k_ij - f(t_ij, x_ij, u_ij, v) (with lagrange part)
        # with x_ij = x_i + sum_l a_il k_jl and assuming u_ij = u_i
        # depends on x_i, u_i, k_i, and v; skip l_i (could skip k_ij[n+1] too...)
        add_nonzero_block!(Is, Js, xi_start, xi_end, xi_start, xi_end)
        add_nonzero_block!(Is, Js, ui_start, ki_end, ui_start, ki_end)
        add_nonzero_block!(Is, Js, xi_start, xi_end, ui_start, ki_end; sym=true)
        add_nonzero_block!(Is, Js, xi_start, xi_end, v_start, v_end; sym=true)
        add_nonzero_block!(Is, Js, ui_start, ki_end, v_start, v_end; sym=true)

        # 1.4 path constraint g(t_i, x_i, u_i, v)
        # -> included in 1.3
    end

    # 2. final path constraints (xf, uf, v) (assume present)
    var_offset = docp.dim_NLP_steps*docp.discretization._step_variables_block
    xf_start = var_offset + 1
    xf_end = var_offset + docp.dim_OCP_x
    uf_start = var_offset + docp.dim_NLP_x + 1
    uf_end = var_offset + docp.dim_NLP_x + docp.dim_NLP_u
    add_nonzero_block!(Is, Js, xf_start, xf_end, xf_start, xf_end)
    add_nonzero_block!(Is, Js, uf_start, uf_end, uf_start, uf_end)
    add_nonzero_block!(Is, Js, xf_start, xf_end, uf_start, uf_end; sym=true) 
    add_nonzero_block!(Is, Js, xf_start, xf_end, v_start, v_end; sym=true)
    add_nonzero_block!(Is, Js, uf_start, uf_end, v_start, v_end; sym=true)

    # 3. boundary constraints (x0, xf, v) or mayer cost g0(x0, xf, v) (assume present)
    # -> x0 / x0, x0 / v terms included in first loop iteration
    # -> xf / xf, xf / v terms included in 2.
    x0_start = 1
    x0_end = docp.dim_OCP_x
    add_nonzero_block!(Is, Js, x0_start, x0_end, xf_start, xf_end; sym=true)

    # 3.1 null initial condition for lagrangian cost state l0
    # -> 2nd order term is zero
   
    # build and return sparse matrix
    nnzj = length(Is)
    Vs = ones(Bool, nnzj)
    return sparse(Is, Js, Vs, docp.dim_NLP_variables, docp.dim_NLP_variables)

end

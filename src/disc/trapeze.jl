#= Functions for trapeze discretization scheme
Internal layout for NLP variables: 
[X_1,U_1, .., X_N+1,U_N+1, V]
=#

# NB. could also be defined as a generic IRK for testing
struct Trapeze <: Discretization

    info::String
    _step_variables_block::Int
    _state_stage_eqs_block::Int
    _step_pathcons_block::Int

    # constructor
    function Trapeze(dim_NLP_steps, dim_NLP_x, dim_NLP_u, dim_NLP_v, dim_u_cons, dim_x_cons, dim_xu_cons, dim_boundary_cons, dim_v_cons)

        # aux variables
        step_variables_block = dim_NLP_x + dim_NLP_u
        state_stage_eqs_block = dim_NLP_x
        step_pathcons_block = dim_u_cons + dim_x_cons + dim_xu_cons

        # NLP variables size ([state, control]_1..N+1, variable)
        dim_NLP_variables = (dim_NLP_steps + 1) * step_variables_block + dim_NLP_v

        # NLP constraints size ([dynamics, stage, path]_1..N, final path, boundary, variable)
        dim_NLP_constraints = dim_NLP_steps * (state_stage_eqs_block + step_pathcons_block) + step_pathcons_block + dim_boundary_cons + dim_v_cons

        disc = new("Implicit Trapeze aka Crank-Nicolson, 2nd order, A-stable", step_variables_block, state_stage_eqs_block, step_pathcons_block)

        return disc, dim_NLP_variables, dim_NLP_constraints
    end
end


"""
$(TYPEDSIGNATURES)

Retrieve state variables at given time step from the NLP variables.
Convention: 1 <= i <= dim_NLP_steps+1
Scalar / Vector output
"""
function get_OCP_state_at_time_step(xu, docp::DOCP{Trapeze, ScalVariable, <: ScalVect, <: ScalVect}, i)
    offset = (i-1) * (docp.dim_NLP_x + docp.dim_NLP_u)
    return xu[offset+1]
end
function get_OCP_state_at_time_step(xu, docp::DOCP{Trapeze, VectVariable, <: ScalVect, <: ScalVect}, i)
    offset = (i-1) * (docp.dim_NLP_x + docp.dim_NLP_u)
    return @view xu[(offset + 1):(offset + docp.dim_OCP_x)]
end
"""
$(TYPEDSIGNATURES)

Retrieve state variable for lagrange cost at given time step from the NLP variables.
Convention: 1 <= i <= dim_NLP_steps+1
"""
function get_lagrange_state_at_time_step(xu, docp::DOCP{Trapeze}, i)
    offset = (i-1) * (docp.dim_NLP_x + docp.dim_NLP_u)
    return xu[offset + docp.dim_NLP_x]
end
"""
$(TYPEDSIGNATURES)

Retrieve control variables at given time step from the NLP variables.
Convention: 1 <= i <= dim_NLP_steps+1
Scalar / Vector output
"""
function get_OCP_control_at_time_step(xu, docp::DOCP{Trapeze, <: ScalVect, ScalVariable, <: ScalVect}, i)
    offset = (i-1) * (docp.dim_NLP_x + docp.dim_NLP_u) + docp.dim_NLP_x
    return xu[offset+1]
end
function get_OCP_control_at_time_step(xu, docp::DOCP{Trapeze, <: ScalVect, VectVariable, <: ScalVect}, i)
    offset = (i-1) * (docp.dim_NLP_x + docp.dim_NLP_u) + docp.dim_NLP_x
    return @view xu[(offset + 1):(offset + docp.dim_NLP_u)]
end

"""
$(TYPEDSIGNATURES)

Set initial guess for state variables at given time step
Convention: 1 <= i <= dim_NLP_steps+1
"""
function set_state_at_time_step!(xu, x_init, docp::DOCP{Trapeze}, i)
    # initialize only actual state variables from OCP (not lagrange state)
    if !isnothing(x_init)
        offset = (i-1) * (docp.dim_NLP_x + docp.dim_NLP_u)
        xu[(offset + 1):(offset + docp.dim_OCP_x)] .= x_init
    end
end
"""
$(TYPEDSIGNATURES)

Set initial guess for control variables at given time step
Convention: 1 <= i <= dim_NLP_steps+1
"""
function set_control_at_time_step!(xu, u_init, docp::DOCP{Trapeze}, i)
    if !isnothing(u_init)
        offset = (i-1) * (docp.dim_NLP_x + docp.dim_NLP_u) + docp.dim_NLP_x
        xu[(offset + 1):(offset + docp.dim_NLP_u)] .= u_init
    end
end

"""
$(TYPEDSIGNATURES)

Set work array for all dynamics and lagrange cost evaluations
"""
function setWorkArray(docp::DOCP{Trapeze}, xu, time_grid, v)

    #= use work array to store all dynamics + lagrange costs
    NB. using a smaller work array to store a single dynamics between steps
    appears slower, maybe due to the copy involved ? =# 
    work = similar(xu, docp.dim_NLP_x * (docp.dim_NLP_steps+1))

    # loop over time steps
    for i = 1:docp.dim_NLP_steps+1
        offset = (i-1) * docp.dim_NLP_x
        ti = time_grid[i]
        xi = get_OCP_state_at_time_step(xu, docp, i)
        ui = get_OCP_control_at_time_step(xu, docp, i)
        # OCP dynamics
        docp.ocp.dynamics((@view work[offset+1:offset+docp.dim_OCP_x]), ti, xi, ui, v)
        # lagrange cost
        if docp.is_lagrange
            docp.ocp.lagrange((@view work[offset+docp.dim_NLP_x:offset+docp.dim_NLP_x]), ti, xi, ui, v)
        end
    end
    return work
end


"""
$(TYPEDSIGNATURES)

Set the constraints corresponding to the state equation
Convention: 1 <= i <= dim_NLP_steps+1
"""
function setStepConstraints!(docp::DOCP{Trapeze}, c, xu, v, time_grid, i, work)

    # offset for previous steps
    offset = (i-1)*(docp.discretization._state_stage_eqs_block + docp.discretization._step_pathcons_block)
    # c block version 
    # offset = 0

    # 0. variables
    ti = time_grid[i]
    xi = get_OCP_state_at_time_step(xu, docp, i)

    # 1. state equation
    if i <= docp.dim_NLP_steps
        # more variables
        tip1 = time_grid[i+1]
        xip1 = get_OCP_state_at_time_step(xu, docp, i+1)
        half_hi = 0.5 * (tip1 - ti)
        offset_dyn_i = (i-1)*docp.dim_NLP_x
        offset_dyn_ip1 = i*docp.dim_NLP_x

        # trapeze rule (no allocations ^^)
        @views @. c[offset+1:offset+docp.dim_OCP_x] = xip1 - (xi + half_hi * (work[offset_dyn_i+1:offset_dyn_i+docp.dim_OCP_x] + work[offset_dyn_ip1+1:offset_dyn_ip1+docp.dim_OCP_x]))

        if docp.is_lagrange
            c[offset+docp.dim_NLP_x] = get_lagrange_state_at_time_step(xu, docp, i+1) - (get_lagrange_state_at_time_step(xu, docp, i) + half_hi * (work[offset_dyn_i+docp.dim_NLP_x] + work[offset_dyn_ip1+docp.dim_NLP_x]))
        end
        offset += docp.dim_NLP_x
    end

    # 2. path constraints
    ui = get_OCP_control_at_time_step(xu, docp, i)
    setPathConstraints!(docp, c, ti, xi, ui, v, offset)

end


"""
$(TYPEDSIGNATURES)

Build sparsity pattern for Jacobian of constraints
"""
function DOCP_Jacobian_pattern(docp::DOCP{Trapeze})

    # vector format for sparse matrix
    Is = Vector{Int}(undef, 0)
    Js = Vector{Int}(undef, 0)

    # index alias for v
    v_start = docp.dim_NLP_variables - docp.dim_NLP_v + 1
    v_end = docp.dim_NLP_variables

    # 1. main loop over steps
    for i = 1:docp.dim_NLP_steps

        # constraints block and offset: state equation, path constraints
        c_block = docp.discretization._state_stage_eqs_block + docp.discretization._step_pathcons_block
        c_offset = (i-1)*c_block

        # contiguous variables blocks will be used when possible 
        # x_i (l_i) u_i x_i+1 (l_i+1) u_i+1
        var_block = docp.discretization._step_variables_block * 2
        var_offset = (i-1)*docp.discretization._step_variables_block
        xi_start = var_offset + 1
        xi_end = var_offset + docp.dim_OCP_x
        ui_start = var_offset + docp.dim_NLP_x + 1
        ui_end = var_offset + docp.dim_NLP_x + docp.dim_NLP_u
        xip1_end = var_offset + docp.dim_NLP_x + docp.dim_NLP_u + docp.dim_OCP_x
        uip1_start = var_offset + docp.dim_NLP_x*2 + docp.dim_NLP_u + 1
        uip1_end = var_offset + docp.dim_NLP_x*2 + docp.dim_NLP_u*2

        # 1.1 state eq 0 = x_i+1 - (x_i + h_i/2 (f(t_i,x_i,u_i,v) + f(t_i+1,x_i+1,u_i+1,v)))
        # depends on x_i, u_i, x_i+1, u_i+1; skip l_i, l_i+1; v cf 1.4
        add_nonzero_block!(Is, Js, c_offset+1, c_offset+docp.dim_OCP_x, xi_start, xi_end)
        add_nonzero_block!(Is, Js, c_offset+1, c_offset+docp.dim_OCP_x, ui_start, xip1_end)
        add_nonzero_block!(Is, Js, c_offset+1, c_offset+docp.dim_OCP_x, uip1_start, uip1_end)
        # 1.2 lagrange part 0 = l_i+1 - (l_i + h_i/2 (l(t_i,x_i,u_i,v) + l(t_i+1,x_i+1,u_i+1,v)))
        # depends on x_i, l_i, u_i, x_i+1, l_i+1, u_i+1 ie whole variable block; v cf 1.4
        if docp.is_lagrange
            add_nonzero_block!(Is, Js, c_offset+docp.dim_NLP_x, c_offset+docp.dim_NLP_x, var_offset+1, var_offset+var_block)
        end

        # 1.3 path constraint g(t_i, x_i, u_i, v) 
        # depends on x_i, u_i; skip l_i; v cf 1.4
        add_nonzero_block!(Is, Js, c_offset+docp.dim_NLP_x+1, c_offset+c_block, xi_start, xi_end)
        add_nonzero_block!(Is, Js, c_offset+docp.dim_NLP_x+1, c_offset+c_block, ui_start, ui_end)

        # 1.4 whole constraint block depends on v
        add_nonzero_block!(Is, Js, c_offset+1, c_offset+c_block, v_start, v_end)
    end

    # 2. final path constraints (xf, uf, v)
    c_offset = docp.dim_NLP_steps*(docp.discretization._state_stage_eqs_block + docp.discretization._step_pathcons_block)
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
    if docp.is_lagrange
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
function DOCP_Hessian_pattern(docp::DOCP{Trapeze})

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
        # x_i (l_i) u_i x_i+1 (l_i+1) u_i+1
        var_block = docp.discretization._step_variables_block * 2
        var_offset = (i-1)*docp.discretization._step_variables_block

        # 1.1 state eq 0 = x_i+1 - (x_i + h_i/2 (f(t_i,x_i,u_i,v) + f(t_i+1,x_i+1,u_i+1,v)))
        # depends on x_i, u_i, x_i+1, u_i+1, and v -> included in 1.2
        # 1.2 lagrange part 0 = l_i+1 - (l_i + h_i/2 (l(t_i,x_i,u_i,v) + l(t_i+1,x_i+1,u_i+1,v)))
        # depends on x_i, l_i, u_i, x_i+1, l_i+1, u_i+1, and v
        # -> use single block for all step variables
        add_nonzero_block!(Is, Js, var_offset+1, var_offset+var_block, var_offset+1, var_offset+var_block)
        add_nonzero_block!(Is, Js, var_offset+1, var_offset+var_block, v_start, v_end; sym=true)

        # 1.3 path constraint g(t_i, x_i, u_i, v)
        # -> included in 1.2
    end

    # 2. final path constraints (xf, uf, v)
    # -> included in last loop iteration

    # 3. boundary constraints (x0, xf, v)
    # -> (x0, v) terms included in first loop iteration
    # -> (xf, v) terms included in last loop iteration
    if docp.is_mayer || docp.dim_boundary_cons > 0
        var_offset = docp.dim_NLP_steps*docp.discretization._step_variables_block
        x0_start = 1
        x0_end = docp.dim_OCP_x
        xf_start = var_offset + 1
        xf_end = var_offset + docp.dim_OCP_x
        add_nonzero_block!(Is, Js, x0_start, x0_end, xf_start, xf_end; sym=true)
    end
    # 3.1 null initial condition for lagrangian cost state l0
    # -> 2nd order term is zero

    # build and return sparse matrix
    nnzj = length(Is)
    Vs = ones(Bool, nnzj)
    return sparse(Is, Js, Vs, docp.dim_NLP_variables, docp.dim_NLP_variables)

end


"""
$(TYPEDSIGNATURES)

Add block of nonzeros elements to a sparsity pattern 
Format: boolean matrix (M) or index vectors (Is, Js) 
Includes a more compact method for single element case
Option to add the symmetric block also (eg for Hessian)
Note: independent from discretization scheme
"""
function add_nonzero_block!(M, i_start, i_end, j_start, j_end; sym=false)
    M[i_start:i_end, j_start:j_end] .= true
    sym && (M[j_start:j_end, i_start:i_end] .= true)
    return
end
function add_nonzero_block!(M, i, j; sym=false)
    M[i,j] = true
    sym && (M[j,i] = true)
    return
end
function add_nonzero_block!(Is, Js, i_start, i_end, j_start, j_end; sym=false)
    for i=i_start:i_end
        for j=j_start:j_end
            push!(Is, i)
            push!(Js, j)
            sym && push!(Is, j)
            sym && push!(Js, i)
        end
    end
    return
end
function add_nonzero_block!(Is, Js, i, j; sym=false)
    push!(Is, i)
    push!(Js, j)
    sym && push!(Is, j)
    sym && push!(Js, i)
    return
end

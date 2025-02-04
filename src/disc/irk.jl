#= Functions for generic implicit Runge Kutta discretization
Internal layout for NLP variables: 
[X_0, U_0, K_0^1..K_0^s, 
 X_1, U_1, K_1^1..K_1^s,
 .., 
 X_N-1, U_N-1, K_N-1^1..K_N-1^s,
 X_N, U_N, V]
with s the stage number and U given by either linear interpolation in [t_i, t_i+1]
or constant interpolation for 1-stage methods or if specfied (U_N might end up unused)
Path constraints are all evaluated at time steps
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

    function Gauss_Legendre_1(dim_NLP_steps, dim_NLP_x, dim_NLP_u, dim_NLP_v, dim_u_cons, dim_x_cons, dim_xu_cons, dim_boundary_cons, dim_v_cons)
        
        stage = 1

        step_variables_block, state_stage_eqs_block, step_pathcons_block, dim_NLP_variables, dim_NLP_constraints = IRK_dims(dim_NLP_steps, dim_NLP_x, dim_NLP_u, dim_NLP_v, dim_u_cons, dim_x_cons, dim_xu_cons, dim_boundary_cons, dim_v_cons, stage)

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
    _constant_control::Bool

    function Gauss_Legendre_2(dim_NLP_steps, dim_NLP_x, dim_NLP_u, dim_NLP_v, dim_u_cons, dim_x_cons, dim_xu_cons, dim_boundary_cons, dim_v_cons, constant_control)
        
        stage = 2

        step_variables_block, state_stage_eqs_block, step_pathcons_block, dim_NLP_variables, dim_NLP_constraints =  IRK_dims(dim_NLP_steps, dim_NLP_x, dim_NLP_u, dim_NLP_v, dim_u_cons, dim_x_cons, dim_xu_cons, dim_boundary_cons, dim_v_cons, stage)

        disc = new("Implicit Gauss-Legendre collocation for s=2, 4th order, symplectic, A-stable",stage,
        [0.25 (0.25-sqrt(3) / 6); (0.25+sqrt(3) / 6) 0.25],
        [0.5, 0.5],
        [(0.5 - sqrt(3) / 6), (0.5 + sqrt(3) / 6)],
        step_variables_block, state_stage_eqs_block, step_pathcons_block,
        constant_control
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
    _constant_control::Bool

    function Gauss_Legendre_3(dim_NLP_steps, dim_NLP_x, dim_NLP_u, dim_NLP_v, dim_u_cons, dim_x_cons, dim_xu_cons, dim_boundary_cons, dim_v_cons, constant_control)
        
        stage = 3

        step_variables_block, state_stage_eqs_block, step_pathcons_block, dim_NLP_variables, dim_NLP_constraints =  IRK_dims(dim_NLP_steps, dim_NLP_x, dim_NLP_u, dim_NLP_v, dim_u_cons, dim_x_cons, dim_xu_cons, dim_boundary_cons, dim_v_cons, stage)

        disc = new("Implicit Gauss-Legendre collocation for s=3, 6th order, symplectic, A-stable",stage,
        [(5.0/36.0) (2/9 - sqrt(15) / 15) (5/36 - sqrt(15) / 30); 
        (5.0/36.0 + sqrt(15) / 24) (2.0/9.0) (5.0/36.0 - sqrt(15) / 24); 
        (5/36 + sqrt(15) / 30) (2/9 + sqrt(15) / 15) (5.0/36.0)],
        [5.0/18.0, 4.0/9.0, 5.0/18.0],
        [0.5 - 0.1*sqrt(15), 0.5, 0.5 + 0.1*sqrt(15)],
        step_variables_block, state_stage_eqs_block, step_pathcons_block, constant_control
        )

        return disc, dim_NLP_variables, dim_NLP_constraints
    end
end


"""
$(TYPEDSIGNATURES)

Return the dimension of the NLP variables and constraints for a generic IRK discretizion, with the control taken constant per step (ie not distinct controls at time stages)
"""
function IRK_dims(dim_NLP_steps, dim_NLP_x, dim_NLP_u, dim_NLP_v, dim_u_cons, dim_x_cons, dim_xu_cons, dim_boundary_cons, dim_v_cons, stage)

    # size of variables block for one step: x, u, k
    step_variables_block = dim_NLP_x + dim_NLP_u + dim_NLP_x * stage

    # size of state + stage equations for one step
    state_stage_eqs_block = dim_NLP_x * (1 + stage)

    # size of path constraints block for one step: u, x, xu 
    step_pathcons_block = dim_u_cons + dim_x_cons + dim_xu_cons

    # NLP variables size ([state, control, stage]_1..N, final state and control, variable)
    dim_NLP_variables = dim_NLP_steps * step_variables_block + dim_NLP_x + dim_NLP_u + dim_NLP_v

    # NLP constraints size ([dynamics, stage, path]_1..N, final path, boundary, variable)
    dim_NLP_constraints = dim_NLP_steps * (state_stage_eqs_block + step_pathcons_block) + step_pathcons_block + dim_boundary_cons + dim_v_cons

    return step_variables_block, state_stage_eqs_block, step_pathcons_block, dim_NLP_variables, dim_NLP_constraints
end


"""
$(TYPEDSIGNATURES)

Retrieve state variables at given time step from the NLP variables.
Convention: 1 <= i <= dim_NLP_steps+1
Scalar / Vector output
"""
function get_OCP_state_at_time_step(xu, docp::DOCP{ <: GenericIRK, ScalVariable, <: ScalVect, <: ScalVect}, i)
    offset = (i-1) * docp.discretization._step_variables_block
    return xu[offset+1]
end
function get_OCP_state_at_time_step(xu, docp::DOCP{ <: GenericIRK, VectVariable, <: ScalVect, <: ScalVect}, i)
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
Scalar / Vector output
Step / Stage versions
"""
function get_OCP_control_at_time_step(xu, docp::DOCP{ <: GenericIRK, <: ScalVect, ScalVariable, <: ScalVect}, i)
    offset = (i-1) * docp.discretization._step_variables_block + docp.dim_NLP_x
    return xu[offset+1]
end
function get_OCP_control_at_time_step(xu, docp::DOCP{ <: GenericIRK, <: ScalVect, VectVariable, <: ScalVect}, i)
    offset = (i-1) * docp.discretization._step_variables_block + docp.dim_NLP_x
    return @view xu[(offset + 1):(offset + docp.dim_NLP_u)]
end
function get_OCP_control_at_time_stage(xu, docp::DOCP{ <: GenericIRK, <: ScalVect, ScalVariable, <: ScalVect}, i, cj)
    if (docp.discretization.stage == 1) || (docp.discretization._constant_control)
        # constant interpolation on step
        return get_OCP_control_at_time_step(xu, docp, i)
    else
        # linear interpolation on step
        ui = get_OCP_control_at_time_step(xu, docp, i)
        uip1 = get_OCP_control_at_time_step(xu, docp, i+1)
        return (1 - cj) * ui + cj * uip1
    end
end
function get_OCP_control_at_time_stage(xu, docp::DOCP{ <: GenericIRK, <: ScalVect, VectVariable, <: ScalVect}, i, cj)
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
    # work array layout: [x_ij ; sum_bk ; u_ij] ?    
    work = similar(xu, docp.dim_OCP_x + docp.dim_NLP_x + docp.dim_NLP_u)
    return work
end

"""
$(TYPEDSIGNATURES)

Set the constraints corresponding to the state equation
Convention: 1 <= i <= dim_NLP_steps (+1)
"""
function setStepConstraints!(docp::DOCP{ <: GenericIRK}, c, xu, v, time_grid, i, work)

    # work array layout: [x_ij ; sum_bk ; u_ij] ?
    work_xij = @view work[1:docp.dim_OCP_x]
    work_sumbk = @view work[docp.dim_OCP_x+1:docp.dim_OCP_x+docp.dim_NLP_x]
    #work_sumbk .= zero(eltype(xu)) AD bug when affecting constant values...
    @views @. work_sumbk[1:docp.dim_NLP_x] = xu[1:docp.dim_NLP_x] * 0.
    #work_uij ?

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
            @views @. work_sumbk[1:docp.dim_NLP_x] = work_sumbk[1:docp.dim_NLP_x] + docp.discretization.butcher_b[j] * kij[1:docp.dim_NLP_x]

            # state at stage: x_i^j = x_i + h_i sum a_jl k_i^l
            # +++ still some allocations here
            @views @. work_xij[1:docp.dim_OCP_x] = xi
            for l = 1:docp.discretization.stage
                kil = get_stagevars_at_time_step(xu, docp, i, l)
                @views @. work_xij[1:docp.dim_OCP_x] = work_xij[1:docp.dim_OCP_x] + hi * docp.discretization.butcher_a[j,l] * kil[1:docp.dim_OCP_x]
            end
            if docp.dim_OCP_x == 1
                xij = work_xij[1]
            else
                xij = work_xij
            end

            # control at stage: interpolation between u_i and u_i+1 
            # +++ use work aray to reduce allocs ?
            uij = get_OCP_control_at_time_stage(xu, docp, i, cj)

            # stage equations k_i^j = f(t_i^j, x_i^j, u_i, v) as c[] = k - f
            # NB. we skip the state equation here, which will be set below
            docp.ocp.dynamics((@view c[offset+offset_stage_eqs+1:offset+offset_stage_eqs+docp.dim_OCP_x]), tij, xij, uij, v)
            if docp.is_lagrange
                docp.ocp.lagrange((@view c[offset+offset_stage_eqs+docp.dim_NLP_x:offset+offset_stage_eqs+docp.dim_NLP_x]), tij, xij, uij, v)
            end
            @views @. c[offset+offset_stage_eqs+1:offset+offset_stage_eqs+docp.dim_NLP_x] = kij - c[offset+offset_stage_eqs+1:offset+offset_stage_eqs+docp.dim_NLP_x]
            offset_stage_eqs += docp.dim_NLP_x

        end

        # state equation x_i+1 = x_i + h_i sum b_j k_i^j
        @views @. c[offset+1:offset+docp.dim_OCP_x] = xip1 - (xi + hi * work_sumbk[1:docp.dim_OCP_x])
        if docp.is_lagrange
            c[offset+docp.dim_NLP_x] = get_lagrange_state_at_time_step(xu, docp, i+1) - (get_lagrange_state_at_time_step(xu, docp, i) + hi * work_sumbk[docp.dim_NLP_x])
        end

        # update offset for stage and state equations
        offset += docp.dim_NLP_x * (1 + docp.discretization.stage)
    end

    #2. path constraints
    ui = get_OCP_control_at_time_step(xu, docp, i)
    setPathConstraints!(docp, c, ti, xi, ui, v, offset)

end


"""
$(TYPEDSIGNATURES)

Build sparsity pattern for Jacobian of constraints
"""
function DOCP_Jacobian_pattern(docp::DOCP{ <: GenericIRK})

    J = zeros(Bool, docp.dim_NLP_constraints, docp.dim_NLP_variables)

    s = docp.discretization.stage

    # 1. main loop over steps
    for i = 1:docp.dim_NLP_steps

        # constraints block and offset: state equation, path constraints
        c_block = docp.discretization._state_stage_eqs_block + docp.discretization._step_pathcons_block
        c_offset = (i-1)*c_block

        # variables block and offset: x_i (l_i) u_i k_i x_i+1 (l_i+1)
        var_block = docp.discretization._step_variables_block + docp.dim_NLP_x
        var_offset = (i-1)*docp.discretization._step_variables_block

        # state eq x_i+1 = x_i + h sum bj k_ij
        # 1.1 state eq wrt x_i
        J[c_offset+1:c_offset+docp.dim_OCP_x, var_offset+1:var_offset+docp.dim_OCP_x] .= true
        # 1.2 state eq wrt k_i, x_i+1 (skip l_i, u_i)
        J[c_offset+1:c_offset+docp.dim_OCP_x, var_offset+docp.dim_NLP_x+docp.dim_NLP_u+ 1:var_offset+docp.dim_NLP_x+docp.dim_NLP_u+s*docp.dim_NLP_x+docp.dim_OCP_x] .= true
        # 1.3 lagrange part l_i+1 = l_i + h (sum bj k_ij)[n+1]
        if docp.is_lagrange
            J[c_offset+docp.dim_NLP_x, var_offset+docp.dim_NLP_x] = true # l_i
            J[c_offset+docp.dim_NLP_x, var_offset+docp.discretization._step_variables_block+docp.dim_NLP_x] = true # l_i+1
            for i=1:s
                J[c_offset+docp.dim_NLP_x, var_offset+(s+1)*docp.dim_NLP_x] = true # k_ij[n+1]
            end
        end
        
        # 1.4 stage equations k_ij = f(t_ij, x_ij, u_ij, v)
        # with
        # x_ij depending on x_i and all k_ij
        # u_ij depending on u_i for piecewise constant or (u_i, u_i+1) for piecewise linear
        # ie whole block depends on x_i, u_i, k_i, and u_i+1 for piecewise linear control
        # NB we could skip l_i here...
        J[c_offset+docp.dim_NLP_x+1:c_offset+(s+1)*docp.dim_NLP_x , var_offset+1:var_offset+docp.discretization._step_variables_block] .= true
        if !docp.discretization._constant_control
            J[c_offset+docp.dim_NLP_x+1:c_offset+(s+1)*docp.dim_NLP_x , var_offset+docp.discretization._step_variables_block+docp.dim_NLP_x:var_offset+docp.discretization._step_variables_block+docp.dim_NLP_x+docp.dim_NLP_u] .= true
        end

        # 1.5 path constraint wrt x_i, u_i
        J[c_offset+(s+1)*docp.dim_NLP_x+1:c_offset+c_block, var_offset+1:var_offset+docp.dim_OCP_x] .= true
        J[c_offset+(s+1)*docp.dim_NLP_x+1:c_offset+c_block, var_offset+docp.dim_NLP_x+1:var_offset+docp.dim_NLP_x+docp.dim_NLP_u] .= true
        
        # 1.6 whole block wrt v
        J[c_offset+1:c_offset+c_block, docp.dim_NLP_variables-docp.dim_NLP_v+1:docp.dim_NLP_variables] .= true
    end

    # 2. final path constraints (xf, uf, v)
    c_offset = docp.dim_NLP_steps*(docp.discretization._state_stage_eqs_block + docp.discretization._step_pathcons_block)
    c_block = docp.discretization._step_pathcons_block
    var_offset = docp.dim_NLP_steps*docp.discretization._step_variables_block
    var_block = docp.discretization._step_variables_block
    # 2.1 wrt xf
    J[c_offset+1:c_offset+c_block, var_offset+1:var_offset+docp.dim_OCP_x] .= true
    # 2.2 wrt uf
    J[c_offset+1:c_offset+c_block, var_offset+docp.dim_NLP_x+1:var_offset+docp.dim_NLP_x+docp.dim_NLP_u] .= true
    # 2.3 wrt v
    J[c_offset+1:c_offset+c_block, docp.dim_NLP_variables-docp.dim_NLP_v+1:docp.dim_NLP_variables] .= true

    # 3. boundary constraints (x0, xf, v)
    c_offset = docp.dim_NLP_steps * (docp.discretization._state_stage_eqs_block + docp.discretization._step_pathcons_block) + docp.discretization._step_pathcons_block
    c_block = docp.dim_boundary_cons + docp.dim_v_cons
    J[c_offset+1:c_offset+c_block, 1:docp.dim_OCP_x] .= true # x0
    J[c_offset+1:c_offset+c_block, var_offset+1:var_offset+docp.dim_OCP_x] .= true # xf
    J[c_offset+1:c_offset+c_block, docp.dim_NLP_variables-docp.dim_NLP_v+1:docp.dim_NLP_variables] .= true # v
    # 3.4 null initial condition for lagrangian cost state l0
    if docp.is_lagrange
        J[docp.dim_NLP_constraints, docp.dim_NLP_x] = true
    end

    # return sparse matrix
    return sparse(J)
end


"""
$(TYPEDSIGNATURES)

Build sparsity pattern for Hessian of Lagrangian
"""
function DOCP_Hessian_pattern(docp::DOCP{ <: GenericIRK})

    # NB. need to provide full pattern for coloring, not just upper/lower part
    H = zeros(Bool, docp.dim_NLP_variables, docp.dim_NLP_variables)
   
    s = docp.discretization.stage

    # 0. objective
    # 0.1 mayer cost (x0, xf, v) 
    # -> grouped with term 3. for boundary conditions
    # 0.2 lagrange case (lf)
    if docp.is_lagrange
        lf_index = docp.dim_NLP_steps * docp.discretization._step_variables_block + docp.dim_NLP_x
        H[lf_index, lf_index] = true
    end
   
    # 1. main loop over steps
    # 1.0 v / v term
    H[docp.dim_NLP_variables-docp.dim_NLP_v+1:docp.dim_NLP_variables, docp.dim_NLP_variables-docp.dim_NLP_v+1:docp.dim_NLP_variables] .= true

    for i = 1:docp.dim_NLP_steps

        # variables block and offset: x_i (l_i) u_i k_i x_i+1 (l_i+1)
        var_block = docp.discretization._step_variables_block + docp.dim_NLP_x
        var_offset = (i-1)*docp.discretization._step_variables_block

        # 1.1 state eq x_i+1 = x_i + h sum bj k_ij
        # wrt x_i, k_i, x_i+1 (skip l_i, u_i)
        # -> included in 1.3 except x_i+1 terms
        H[var_offset+1:var_offset+docp.dim_OCP_x, var_offset+docp.discretization._step_variables_block+1:var_offset+docp.discretization._step_variables_block+docp.dim_OCP_x] .= true # x_i / x_i+1
        H[var_offset+docp.discretization._step_variables_block+1:var_offset+docp.discretization._step_variables_block+docp.dim_OCP_x , var_offset+1:var_offset+docp.dim_OCP_x] .= true # x_i+1 / x_i
        H[var_offset+docp.dim_NLP_x+1:var_offset+(s+1)*docp.dim_NLP_x, var_offset+docp.discretization._step_variables_block+1:var_offset+docp.discretization._step_variables_block+docp.dim_OCP_x] .= true # k_i / x_i+1
        H[var_offset+docp.discretization._step_variables_block+1:var_offset+docp.discretization._step_variables_block+docp.dim_OCP_x , var_offset+docp.dim_NLP_x+1:var_offset+(s+1)*docp.dim_NLP_x] .= true # x_i+1 / k_i    

        # 1.2 lagrange part l_i+1 = l_i + h (sum bj k_ij)[n+1]
        # -> included in 1.3 except l_i+1 terms
        # +++ could be done fully here and l_i skipped in 1.3
        H[var_offset+docp.dim_NLP_x , var_offset+docp.discretization._step_variables_block+docp.dim_NLP_x] = true # l_i / l_i+1
        H[ var_offset+docp.discretization._step_variables_block+docp.dim_NLP_x, var_offset+docp.dim_NLP_x] = true # l_i+1 / l_i
        for i=1:s
            H[var_offset+docp.dim_NLP_u+i*docp.dim_NLP_x , var_offset+docp.discretization._step_variables_block+docp.dim_NLP_x] = true # k_i[n+1] / l_i+1
            H[var_offset+docp.discretization._step_variables_block+docp.dim_NLP_x, var_offset+docp.dim_NLP_u+i*docp.dim_NLP_x] = true # l_i+1 / k_i[n+1]
        end

        # 1.3 stage equations k_ij = f(t_ij, x_ij, u_ij, v)
        # wrt x_i, u_i, k_i  (and u_i+1 for piecewise linear control)
        # NB. l_i terms for 1.2 are included but we have excess nnz eg l_i / x_i,u_i,k_i[1:n]
        if docp.discretization._constant_control
            H[var_offset+1:var_offset+docp.discretization._step_variables_block, var_offset+1:var_offset+docp.discretization._step_variables_block] .= true
        else 
            error("Manual Hessian sparsity pattern not supported for IRK scheme with piecewise linear control")
        end

        # 1.4 path constraint wrt x_i, u_i
        # -> included in 1.3

        # 1.5 whole block wrt v (NB. term v / v added before the loop)
        H[docp.dim_NLP_variables-docp.dim_NLP_v+1:docp.dim_NLP_variables, var_offset+1:var_offset+var_block] .= true # v / var block
        H[var_offset+1:var_offset+var_block, docp.dim_NLP_variables-docp.dim_NLP_v+1:docp.dim_NLP_variables] .= true # var block / v
    end

    # 2. final path constraints (xf, uf, v)
    # -> included in last iteration from loop

    # 3. boundary constraints (x0, xf, v)
    # -> x0 / x0, x0 / v, xf / xf, xf / v terms included in first/last iterations from loop
    if docp.is_mayer || docp.dim_boundary_cons > 0
        var_offset = docp.dim_NLP_steps*docp.discretization._step_variables_block
        H[1:docp.dim_OCP_x, var_offset+1:var_offset+docp.dim_OCP_x] .= true # x0 / xf
        H[var_offset+1:var_offset+docp.dim_OCP_x, 1:docp.dim_OCP_x] .= true # xf / x0
        H[docp.dim_NLP_variables-docp.dim_NLP_v+1:docp.dim_NLP_variables, 1:docp.dim_OCP_x] .= true # v / x0
        H[docp.dim_NLP_variables-docp.dim_NLP_v+1:docp.dim_NLP_variables, var_offset+1:var_offset+docp.dim_OCP_x] .= true # v / xf
    end
    # 3.1 null initial condition for lagrangian cost state l0
    if docp.is_lagrange
        H[docp.dim_NLP_x, docp.dim_NLP_x] = true
    end

    # return sparse matrix
    return sparse(H)

end
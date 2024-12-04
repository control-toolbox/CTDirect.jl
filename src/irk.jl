#= Functions for generic implicit Runge Kutta discretization
Internal layout for NLP variables: 
[X_0, U_0, K_0^1..K_0^s, 
 X_1, U_1, K_1^1..K_1^s,
 .., 
 X_N-1, U_N-1, K_N-1^1..K_N-1^s,
 X_N, U_N, V]
with s the stage number and U given by linear interpolation in [t_i, t_i+1]
Path constraints are all evaluated at time steps
=#


abstract type GenericIRK <: Discretization end

"""
$(TYPEDSIGNATURES)

Implicit Midpoint discretization, formulated as a generic IRK
NB. does not use the simplification xs = 0.5 * (xi + xip1) as in midpoint.jl
"""
struct Midpoint_IRK <: GenericIRK

    info::String
    stage::Int
    butcher_a::Matrix{Float64}
    butcher_b::Vector{Float64}
    butcher_c::Vector{Float64}
    _step_variables_block::Int
    _state_stage_eqs_block::Int
    _step_pathcons_block::Int

    function Midpoint_IRK(dim_NLP_steps, dim_NLP_x, dim_NLP_u, dim_NLP_v, dim_u_cons, dim_x_cons, dim_xu_cons, dim_boundary_cons, dim_v_cons)
        
        stage = 1

        step_variables_block, state_stage_eqs_block, step_pathcons_block, dim_NLP_variables, dim_NLP_constraints = IRK_dims(dim_NLP_steps, dim_NLP_x, dim_NLP_u, dim_NLP_v, dim_u_cons, dim_x_cons, dim_xu_cons, dim_boundary_cons, dim_v_cons, stage)

        disc = new("Implicit Midpoint aka Gauss-Legendre collocation for s=1, 2nd order, symplectic", stage, hcat(0.5), [1], [0.5], step_variables_block, state_stage_eqs_block, step_pathcons_block)

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

    function Gauss_Legendre_2(dim_NLP_steps, dim_NLP_x, dim_NLP_u, dim_NLP_v, dim_u_cons, dim_x_cons, dim_xu_cons, dim_boundary_cons, dim_v_cons)
        
        stage = 2

        step_variables_block, state_stage_eqs_block, step_pathcons_block, dim_NLP_variables, dim_NLP_constraints =  IRK_dims(dim_NLP_steps, dim_NLP_x, dim_NLP_u, dim_NLP_v, dim_u_cons, dim_x_cons, dim_xu_cons, dim_boundary_cons, dim_v_cons, stage)

        disc = new("Implicit Gauss-Legendre collocation for s=2, 4th order, symplectic",stage,
        [0.25 (0.25-sqrt(3) / 6); (0.25+sqrt(3) / 6) 0.25],
        [0.5, 0.5],
        [(0.5 - sqrt(3) / 6), (0.5 + sqrt(3) / 6)],
        step_variables_block, state_stage_eqs_block, step_pathcons_block
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
function get_OCP_control_at_time_stage(xu, docp::DOCP{ <: GenericIRK}, i, cj)
    # linear interpolation on step +++ check allocs
    ui = get_OCP_control_at_time_step(xu, docp, i)
    uip1 = get_OCP_control_at_time_step(xu, docp, i+1)
    return (1 - cj) * ui + cj * uip1
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

    # +++ PUT BACK COMPUTATIONS in setconstraintblock to reuse kij and xij better ??
    # use a new branch and compare performance, just ditch work array for irk

    # use work array to store all dynamics + lagrange costs
    # + all the xij variables
    # + one state/stage variable (including lagrange part for setConstraintsBlock)
    work = similar(xu, docp.dim_NLP_x * docp.discretization.stage * docp.dim_NLP_steps 
    + docp.dim_OCP_x * docp.discretization.stage * docp.dim_NLP_steps + docp.dim_NLP_x)
    offset_xij = docp.dim_NLP_x * docp.discretization.stage * docp.dim_NLP_steps

    # loop over time steps ans stages
    for i = 1:docp.dim_NLP_steps
        ti = time_grid[i]
        xi = get_OCP_state_at_time_step(xu, docp, i)
        tip1 = time_grid[i+1]
        hi = tip1 - ti
        for j = 1:docp.discretization.stage
            offset = (i-1) * docp.dim_NLP_x * docp.discretization.stage + (j-1) * docp.dim_NLP_x
            # time at stage: t_i^j = t_i + c[j] h_i
            cj = docp.discretization.butcher_c[j]
            tij = ti + cj * hi
            # control at stage: interpolation between u_i and u_i+1
            uij = get_OCP_control_at_time_stage(xu, docp, i, cj)
            # state at stage: x_i^j = x_i + h_i sum a_jl k_i^l
            @. work[offset_xij+1:offset_xij+docp.dim_OCP_x] = xi
            for l = 1:docp.discretization.stage
                kil = get_stagevars_at_time_step(xu, docp, i, l)
                @views @. work[offset_xij+1:offset_xij+docp.dim_OCP_x] = work[offset_xij+1:offset_xij+docp.dim_OCP_x] + hi * docp.discretization.butcher_a[j,l] * kil[1:docp.dim_OCP_x]
            end
            if docp.dim_OCP_x == 1
                xij = work[offset_xij+1]
            else
                @views xij = work[offset_xij+1:offset_xij+docp.dim_OCP_x] #seems to displace the allocs somehow...
            end
            offset_xij += docp.dim_OCP_x
            # dynamics for stage equation k_i^j = f(t_i^j, x_i^j, u_i, v)
            # +++ put -f in c directly ? just ditch work array for irk !
            docp.ocp.dynamics((@view work[offset+1:offset+docp.dim_OCP_x]), tij, xij, uij, v)
            # lagrange cost
            if docp.is_lagrange
                docp.ocp.lagrange((@view work[offset+docp.dim_NLP_x:offset+docp.dim_NLP_x]), tij, xij, uij, v)
            end
        end
    end
    return work
end

"""
$(TYPEDSIGNATURES)

Set the constraints corresponding to the state equation
Convention: 1 <= i <= dim_NLP_steps (+1)
"""
function setStepConstraints!(docp::DOCP{ <: GenericIRK}, c, xu, v, time_grid, i, work)

    # offset for previous steps
    offset = (i-1)*(docp.dim_NLP_x * (1+docp.discretization.stage) + docp.discretization._step_pathcons_block)

    # 0. variables
    ti = time_grid[i]
    xi = get_OCP_state_at_time_step(xu, docp, i)
    ui = get_OCP_control_at_time_step(xu, docp, i) # for constraints

    # 1. state and stage equations
        if i <= docp.dim_NLP_steps
        # more variables
        tip1 = time_grid[i+1]
        xip1 = get_OCP_state_at_time_step(xu, docp, i+1)        
        hi = tip1 - ti
        offset_dyn_i = (i-1) * docp.dim_NLP_x * docp.discretization.stage
        offset_sumbk = length(work) - docp.dim_NLP_x
        offset_stage_eqs = docp.dim_NLP_x

        # work array for sum b_j k_i^j (w/ lagrange term)
        #@. work[offset_x+1:offset_x+docp.dim_NLP_x] = 0 known AD bug with :optimized backend: cannot affect constants       
        @views @. work[offset_sumbk+1:offset_sumbk+docp.dim_NLP_x] = 0 * work[1:docp.dim_NLP_x]

        # loop over stages
        for j=1:docp.discretization.stage
            kij = get_stagevars_at_time_step(xu, docp, i, j)
            
            # update sum b_j k_i^j (w/ lagrange term) for state equation below
            @views @. work[offset_sumbk+1:offset_sumbk+docp.dim_NLP_x] = work[offset_sumbk+1:offset_sumbk+docp.dim_NLP_x] + docp.discretization.butcher_b[j] * kij[1:docp.dim_NLP_x]

            # stage equations k_i^j = f(t_i^j, x_i^j, u_i, v) cf setWorkArray()
            # NB. we skip the state equation here, which will be set below
            @views @. c[offset+offset_stage_eqs+1:offset+offset_stage_eqs+docp.dim_OCP_x] = kij[1:docp.dim_OCP_x] - work[offset_dyn_i+1:offset_dyn_i+docp.dim_OCP_x]
            if docp.is_lagrange
                c[offset+offset_stage_eqs+docp.dim_NLP_x] = kij[docp.dim_NLP_x] - work[offset_dyn_i+docp.dim_NLP_x]
            end
            offset_stage_eqs += docp.dim_NLP_x

        end

        # state equation x_i+1 = x_i + h_i sum b_j k_i^j
        @views @. c[offset+1:offset+docp.dim_OCP_x] = xip1 - (xi + hi * work[offset_sumbk+1:offset_sumbk+docp.dim_OCP_x])
        if docp.is_lagrange
            c[offset+docp.dim_NLP_x] = get_lagrange_state_at_time_step(xu, docp, i+1) - (get_lagrange_state_at_time_step(xu, docp, i) + hi * work[offset_sumbk+docp.dim_NLP_x])
        end

        # update offset for stage and state equations
        offset += docp.dim_NLP_x * (1 + docp.discretization.stage)
    end

    # 2. path constraints (always at step, move to generic function in problem.jl with arguments ti, xi, ui, v ?)
    # 2.1 control constraints
    if docp.dim_u_cons > 0
        docp.control_constraints[2]((@view c[offset+1:offset+docp.dim_u_cons]),ti, ui, v)
        offset += docp.dim_u_cons
    end

    # 2.2 state constraints
    if docp.dim_x_cons > 0 
        docp.state_constraints[2]((@view c[offset+1:offset+docp.dim_x_cons]),ti, xi, v)
        offset += docp.dim_x_cons
    end

    # 2.3 mixed constraints
    if docp.dim_xu_cons > 0
        docp.mixed_constraints[2]((@view c[offset+1:offset+docp.dim_xu_cons]), ti, xi, ui, v)
        offset += docp.dim_xu_cons
    end

end

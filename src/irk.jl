#= Functions for generic implicit Runge Kutta discretization
Internal layout for NLP variables: 
[X_0, U_0, K_0^1..K_0^s, 
 X_1, U_1, K_1^1..K_1^s,
 .., 
 X_N-1, U_N-1, K_N-1^1..K_N-1^s,
 XN, V]
with s the stage number and U_i taken constant per step
=#

abstract type GenericIRK <: Discretization end


"""
$(TYPEDSIGNATURES)

Return the dimension of the NLP variables and constraints for a generic IRK discretizion, with the control taken constant per step (ie not distinct controls at time stages)
"""
function IRK_dims(dim_NLP_steps, dim_NLP_x, dim_NLP_u, dim_NLP_v, dim_path_cons, dim_boundary_cons, dim_v_cons)
    
    # NLP variables size (state, control, variable, stage)
    dim_NLP_variables = (dim_NLP_steps + 1) * dim_NLP_x + dim_NLP_steps * dim_NLP_u + dim_NLP_v + dim_NLP_steps * dim_NLP_x * stage

    # NLP constraints size (dynamics, stage, path, boundary, variable)
    dim_NLP_constraints = dim_NLP_steps * (dim_NLP_x + (dim_NLP_x * stage) + dim_path_cons) + dim_path_cons + dim_boundary_cons + dim_v_cons

    # size of variables block for one step
    step_block = dim_NLP_x + dim_NLP_u + dim_NLP_x * stage

    return dim_NLP_variables, dim_NLP_constraints, step_block
 end


"""
$(TYPEDSIGNATURES)

Implicit Midpoint discretization, formulated as a generic IRK
"""
struct Midpoint_IRK <: GenericIRK

    stage::Int
    butcher_a::Matrix{Float64}
    butcher_b::Vector{Float64}
    butcher_c::Vector{Float64}

    _step_block::Int

    function Midpoint_IRK(dim_NLP_steps, dim_NLP_x, dim_NLP_u, dim_NLP_v, dim_path_cons, dim_boundary_cons, dim_v_cons)
        
        stage = 1

        dim_NLP_variables, dim_NLP_constraints, step_block = IRK_dims(dim_NLP_steps, dim_NLP_x, dim_NLP_u, dim_NLP_v, dim_path_cons, dim_boundary_cons, dim_v_cons)

        disc = new(stage, hcat(0.5), [1], [0.5], step_block,
        "Implicit Midpoint aka Gauss-Legendre collocation for s=1, 2nd order, symplectic")

        return disc, dim_NLP_variables, dim_NLP_constraints
    end
end

"""
$(TYPEDSIGNATURES)

Gauss Legendre 2 discretization, formulated as a generic IRK
"""
struct GaussLegendre2 <: GenericIRK

    stage::Int
    butcher_a::Matrix{Float64}
    butcher_b::Vector{Float64}
    butcher_c::Vector{Float64}

    function GaussLegendre2(dim_NLP_steps, dim_NLP_x, dim_NLP_u, dim_NLP_v, dim_path_cons, dim_boundary_cons, dim_v_cons)
        
        stage = 2

        disc = new(stage,
        [0.25 (0.25-sqrt(3) / 6); (0.25+sqrt(3) / 6) 0.25],
        [0.5, 0.5],
        [(0.5 - sqrt(3) / 6), (0.5 + sqrt(3) / 6)],
        "Implicit Gauss-Legendre collocation for s=2, 4th order, symplectic")

        dim_NLP_variables, dim_NLP_constraints = IRK_dims(dim_NLP_steps, dim_NLP_x, dim_NLP_u, dim_NLP_v, dim_path_cons, dim_boundary_cons, dim_v_cons)

        return disc, dim_NLP_variables, dim_NLP_constraints
    end
end


"""
$(TYPEDSIGNATURES)

Retrieve state variables at given time step from the NLP variables.
Convention: 1 <= i <= dim_NLP_steps+1
Scalar / Vector output
"""
function get_OCP_state_at_time_step(xu, docp::DOCP{GenericIRK, ScalVariable, <: ScalVect, <: ScalVect}, i)
    offset = (i-1) * docp.discretization._step_block
    return xu[offset+1]
end
function get_OCP_state_at_time_step(xu, docp::DOCP{GenericIRK, VectVariable, <: ScalVect, <: ScalVect}, i)
    offset = (i-1) * docp.discretization._step_block
    return @view xu[(offset + 1):(offset + docp.dim_OCP_x)]
end
"""
$(TYPEDSIGNATURES)

Retrieve state variable for lagrange cost at given time step from the NLP variables.
Convention: 1 <= i <= dim_NLP_steps+1
"""
function get_lagrange_state_at_time_step(xu, docp::DOCP{GenericIRK}, i)
    offset = (i-1) * docp.discretization._step_block
    return xu[offset + docp.dim_NLP_x]
end
"""
$(TYPEDSIGNATURES)

Retrieve control variables at given time step from the NLP variables.
Convention: 1 <= i <= dim_NLP_steps+1
Scalar / Vector output
"""
function get_OCP_control_at_time_step(xu, docp::DOCP{GenericIRK, <: ScalVect, ScalVariable, <: ScalVect}, i)
    offset = (i-1) * docp.discretization._step_block + docp.dim_NLP_x
    # use U_N-1 at tf
    (i == docp.dim_NLP_steps+1) && offset -= docp.discretization._step_block
    return xu[offset+1]
end
function get_OCP_control_at_time_step(xu, docp::DOCP{GenericIRK, <: ScalVect, VectVariable, <: ScalVect}, i)
    offset = (i-1) * docp.discretization._step_block + docp.dim_NLP_x
    # use U_N-1 at tf
    (i == docp.dim_NLP_steps+1) && offset -= docp.discretization._step_block
    return @view xu[(offset + 1):(offset + docp.dim_NLP_u)]
end

"""
$(TYPEDSIGNATURES)

Retrieve stage variables at given time step/stage from the NLP variables.
Convention: 1 <= i <= dim_NLP_steps+1,	1 <= j <= s
"""
function get_stagevars_at_time_step(xu, docp::DOCP{GenericIRK}, i, j)
    if i == docp.dim_NLP_steps+1
        return @view xu[1:docp.dim_NLP_x] # unused but keep same type !
    else
        offset = (i-1) * docp.discretization._step_block + docp.dim_NLP_x + docp.dim_NLP_u + (j-1)*docp.dim_NLP_x
        return @view xu[(offset + 1):(offset + docp.dim_NLP_x)]
    end
end

"""
$(TYPEDSIGNATURES)

Set initial guess for state variables at given time step
Convention: 1 <= i <= dim_NLP_steps+1
"""
function set_state_at_time_step!(xu, x_init, docp::DOCP{GenericIRK}, i)
    offset = (i-1) * docp.discretization._step_block
    # initialize only actual state variables from OCP (not lagrange state)
    if !isnothing(x_init)
        xu[(offset + 1):(offset + docp.dim_OCP_x)] .= x_init
    end
end
"""
$(TYPEDSIGNATURES)

Set initial guess for control variables at given time step
Convention: 1 <= i <= dim_NLP_steps+1
"""
function set_control_at_time_step!(xu, u_init, docp::DOCP{GenericIRK}, i)
    offset = (i-1) * docp.discretization._step_block + docp.dim_NLP_x
    # use U_N-1 at tf
    (i == docp.dim_NLP_steps+1) && offset -= docp.discretization._step_block
    if !isnothing(u_init)
        xu[(offset + 1):(offset + docp.dim_NLP_u)] .= u_init
    end
end

"""
$(TYPEDSIGNATURES)

Set work array for all dynamics and lagrange cost evaluations
"""
function setWorkArray(docp::DOCP{Midpoint}, xu, time_grid, v)

    # use work array to store all dynamics + lagrange costs
    work = similar(xu, docp.dim_NLP_x * (docp.dim_NLP_steps))
    if docp.dim_OCP_x > 1
        xs = similar(xu, docp.dim_OCP_x)
    end

    # loop over time steps ans stages
    for i = 1:docp.dim_NLP_steps
        ti = time_grid[i]
        xi = get_OCP_state_at_time_step(xu, docp, i)
        ui = get_OCP_control_at_time_step(xu, docp, i)
        tip1 = time_grid[i+1]
        hi = tip1 - ti
        for j = 1:docp.discretization.stage
            offset = (i-1) * docp.dim_NLP_x * docp.discretization.stage + (j-1) * docp.dim_NLP_x
            # time at stage: t_i^j = t_i + c[j] h_i
            tij = ti + docp.discretization.butcher_c[j] * hi 
            # state at stage: x_i^j = x_i + h_i sum a_jl k_i^l
            xij = xi
            # +++ check for allocs here !
            for l = 1:docp.discretization.stage
                xij += (hi * docp.discretization.butcher_a[j][l] * get_stagevars_at_time_step(xu, docp, i, l)
            end
            # dynamics for stage equation k_i^j = f(t_i^j, x_i^j, u_i, v) 
            docp.ocp.dynamics((@view work[offset+1:offset+docp.dim_OCP_x]), tij, xij, ui, v)
            # lagrange cost
            if docp.is_lagrange
                docp.ocp.lagrange((@view work[offset+docp.dim_NLP_x:offset+docp.dim_NLP_x]), tij, xij, ui, v)
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
function setConstraintBlock!(docp::DOCP{Midpoint}, c, xu, v, time_grid, i, work)

    # offset for previous steps
    offset = (i-1)*(docp.dim_NLP_x * (1+docp.discretization.stage) + docp.dim_path_cons)

    # 0. variables
    ti = time_grid[i]
    xi = get_OCP_state_at_time_step(xu, docp, i)
    ui = get_OCP_control_at_time_step(xu, docp, i)
    
    # 1. state and stage equations
    if i <= docp.dim_NLP_steps
        # more variables
        tip1 = time_grid[i+1]
        xip1 = get_OCP_state_at_time_step(xu, docp, i+1)        
        hi = tip1 - ti
        offset_dyn_i = (i-1) * docp.dim_NLP_x * docp.discretization.stage

        # state equation x_i+1 = x_i + h_i sum b_j k_i^j
        # probably one alloc here. may be removed by rewriting the state equation
        # ie starting with storing sumbk in c ?
        # maybe fuse the 2 loops over j as well ?
        sum_bki = docp.discretization.butcher_b[1] * get_stagevars_at_time_step(xu, docp, i, 1)
        for j=2:docp.discretization.stage
            sum_bki += docp.discretization.butcher_b[j] * get_stagevars_at_time_step(xu, docp, i, j)
        end

        @views @. c[offset+1:offset+docp.dim_OCP_x] = xip1 - (xi + hi * sum_bki[1:docp.dim_OCP_x])
        if docp.is_lagrange
            c[offset+docp.dim_NLP_x] = get_lagrange_state_at_time_step(xu, docp, i+1) - (get_lagrange_state_at_time_step(xu, docp, i) + hi * sum_bki[docp.dim_NLP_x])
        end
        offset += docp.dim_NLP_x

        # stage equations k_i^j = f(t_i^j, x_i^j, u_i, v) cf setWorkArray()
        for j=1:docp.discretization.stage
            kij = get_stagevars_at_time_step(xu, docp, i, j)
            @views @. c[offset+1:offset+docp.dim_OCP_x] = kij[1:docp.dim_OCP_x] - work[offset_dyn_i+1:offset_dyn_i+docp.dim_OCP_x]
            if docp.is_lagrange
                c[offset+docp.dim_NLP_x] = kij[docp.dim_NLP_x] - work[offset_work+docp.dim_NLP_x]
            end
            offset += docp.dim_NLP_x
        end

    end

    # 2. path constraints +++ use a function in problem.jl ?
    if docp.dim_u_cons > 0
        docp.control_constraints[2]((@view c[offset+1:offset+docp.dim_u_cons]),ti, ui, v)
    end
    if docp.dim_x_cons > 0 
        docp.state_constraints[2]((@view c[offset+docp.dim_u_cons+1:offset+docp.dim_u_cons+docp.dim_x_cons]),ti, xi, v)
    end
    if docp.dim_mixed_cons > 0 
        docp.mixed_constraints[2]((@view c[offset+docp.dim_u_cons+docp.dim_x_cons+1:offset+docp.dim_u_cons+docp.dim_x_cons+docp.dim_mixed_cons]), ti, xi, ui, v)
    end

end

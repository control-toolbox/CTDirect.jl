#= Functions for implicit midpoint discretization scheme
Internal layout for NLP variables: 
[X_1,U_1,K_1 .., X_N,U_N,K_N, X_N+1, V]
with the convention u([t_i,t_i+1[) = U_i and u(tf) is NOT defined
=#

struct Midpoint <: Discretization

    info::String
    _step_pathcons_block::Int
    _state_stage_eqs_block::Int

    # constructor
    function Midpoint(dim_NLP_steps, dim_NLP_x, dim_NLP_u, dim_NLP_v, dim_u_cons, dim_x_cons, dim_xu_cons, dim_boundary_cons, dim_v_cons)

        # NLP variables size ([state, control, stage]_1..N, final state and control, variable)
        dim_NLP_variables = dim_NLP_steps * (dim_NLP_x + dim_NLP_u + dim_NLP_x) + dim_NLP_x + dim_NLP_u + dim_NLP_v
        
        # Path constraints (control, state, mixed)
        state_stage_eqs_block = dim_NLP_x * 2
        step_pathcons_block = dim_u_cons + dim_x_cons + dim_xu_cons

        # NLP constraints size ([dynamics, stage, path]_1..N, final path, boundary, variable)
        dim_NLP_constraints = dim_NLP_steps * (state_stage_eqs_block + step_pathcons_block) + step_pathcons_block + dim_boundary_cons + dim_v_cons

        disc = new("Implicit Midpoint aka Gauss-Legendre collocation for s=1, 2nd order, symplectic", step_pathcons_block, state_stage_eqs_block)

        return disc, dim_NLP_variables, dim_NLP_constraints
    end
end

# NB. some of these are actually the same as for IRK
# later compare performance of this one vs the irk midpoint
# and remove this one if no noticeable gain

"""
$(TYPEDSIGNATURES)

Retrieve state variables at given time step from the NLP variables.
Convention: 1 <= i <= dim_NLP_steps+1
Scalar / Vector output
"""
function get_OCP_state_at_time_step(xu, docp::DOCP{Midpoint, ScalVariable, <: ScalVect, <: ScalVect}, i)
    offset = (i-1) * (docp.dim_NLP_x*2 + docp.dim_NLP_u)
    return xu[offset+1]
end
function get_OCP_state_at_time_step(xu, docp::DOCP{Midpoint, VectVariable, <: ScalVect, <: ScalVect}, i)
    offset = (i-1) * (docp.dim_NLP_x*2+ docp.dim_NLP_u)
    return @view xu[(offset + 1):(offset + docp.dim_OCP_x)]
end
"""
$(TYPEDSIGNATURES)

Retrieve state variable for lagrange cost at given time step from the NLP variables.
Convention: 1 <= i <= dim_NLP_steps+1
"""
function get_lagrange_state_at_time_step(xu, docp::DOCP{Midpoint}, i)
    offset = (i-1) * (docp.dim_NLP_x*2 + docp.dim_NLP_u)
    return xu[offset + docp.dim_NLP_x]
end
"""
$(TYPEDSIGNATURES)

Retrieve control variables at given time step from the NLP variables.
Convention: 1 <= i <= dim_NLP_steps
Scalar / Vector output
"""
function get_OCP_control_at_time_step(xu, docp::DOCP{Midpoint, <: ScalVect, ScalVariable, <: ScalVect}, i)
    offset = (i-1) * (docp.dim_NLP_x*2 + docp.dim_NLP_u) + docp.dim_NLP_x
    return xu[offset+1]
end
function get_OCP_control_at_time_step(xu, docp::DOCP{Midpoint, <: ScalVect, VectVariable, <: ScalVect}, i)
    offset = (i-1) * (docp.dim_NLP_x*2 + docp.dim_NLP_u) + docp.dim_NLP_x
    return @view xu[(offset + 1):(offset + docp.dim_NLP_u)]
end

"""
$(TYPEDSIGNATURES)

Retrieve stage variables at given time step from the NLP variables.
Convention: 1 <= i <= dim_NLP_steps
"""
function get_stagevars_at_time_step(xu, docp::DOCP{Midpoint}, i)
    offset = (i-1) * (docp.dim_NLP_x*2 + docp.dim_NLP_u) + docp.dim_NLP_x + docp.dim_NLP_u
    return @view xu[(offset + 1):(offset + docp.dim_NLP_x)]
end

"""
$(TYPEDSIGNATURES)

Set initial guess for state variables at given time step
Convention: 1 <= i <= dim_NLP_steps+1
"""
function set_state_at_time_step!(xu, x_init, docp::DOCP{Midpoint}, i)
    # initialize only actual state variables from OCP (not lagrange state)
    if !isnothing(x_init)
        offset = (i-1) * (docp.dim_NLP_x*2 + docp.dim_NLP_u)
        xu[(offset + 1):(offset + docp.dim_OCP_x)] .= x_init
    end
end
"""
$(TYPEDSIGNATURES)

Set initial guess for control variables at given time step
Convention: 1 <= i <= dim_NLP_steps+1
"""
function set_control_at_time_step!(xu, u_init, docp::DOCP{Midpoint}, i)
    if !isnothing(u_init)
        offset = (i-1) * (docp.dim_NLP_x*2 + docp.dim_NLP_u) + docp.dim_NLP_x
        xu[(offset + 1):(offset + docp.dim_NLP_u)] .= u_init
    end
end

"""
$(TYPEDSIGNATURES)

Set work array for all dynamics and lagrange cost evaluations
"""
function setWorkArray(docp::DOCP{Midpoint}, xu, time_grid, v)
    
    # NB. compare with version with everything in setStepcontraints

    # use work array to store all dynamics + lagrange costs
    work = similar(xu, docp.dim_NLP_x * (docp.dim_NLP_steps))
    if docp.dim_OCP_x > 1
        xs = similar(xu, docp.dim_OCP_x)
    end

    # loop over time steps
    for i = 1:docp.dim_NLP_steps
        offset = (i-1) * docp.dim_NLP_x
        ti = time_grid[i]
        xi = get_OCP_state_at_time_step(xu, docp, i)
        ui = get_OCP_control_at_time_step(xu, docp, i)
        tip1 = time_grid[i+1]
        xip1 = get_OCP_state_at_time_step(xu, docp, i+1)
        ts = 0.5 * (ti + tip1)
        # add a new getter ?
        if docp.dim_OCP_x == 1
            xs = 0.5 * (xi + xip1)
        else
            @. xs = 0.5 * (xi + xip1) # not ok for dim 1...
        end
        # OCP dynamics
        docp.ocp.dynamics((@view work[offset+1:offset+docp.dim_OCP_x]), ts, xs, ui, v)
        # lagrange cost
        if docp.is_lagrange
            docp.ocp.lagrange((@view work[offset+docp.dim_NLP_x:offset+docp.dim_NLP_x]), ts, xs, ui, v)
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
    offset = (i-1)*(docp.dim_NLP_x * 2 + docp.discretization._step_pathcons_block)

    # 0. variables
    ti = time_grid[i]
    xi = get_OCP_state_at_time_step(xu, docp, i)
    ui = get_OCP_control_at_time_step(xu, docp, i)

    # 1. state equation
    if i <= docp.dim_NLP_steps
        # more variables
        tip1 = time_grid[i+1]
        xip1 = get_OCP_state_at_time_step(xu, docp, i+1)
        hi = tip1 - ti
        ki = get_stagevars_at_time_step(xu, docp, i)
        offset_dyn_i = (i-1)*docp.dim_NLP_x

        # midpoint rule
        @views @. c[offset+1:offset+docp.dim_OCP_x] = xip1 - (xi + hi * ki[1:docp.dim_OCP_x])
        if docp.is_lagrange
            c[offset+docp.dim_NLP_x] = get_lagrange_state_at_time_step(xu, docp, i+1) - (get_lagrange_state_at_time_step(xu, docp, i) + hi * ki[docp.dim_NLP_x])
        end
        offset += docp.dim_NLP_x

        # stage equation at mid-step
        @views @. c[offset+1:offset+docp.dim_OCP_x] = ki[1:docp.dim_OCP_x] - work[offset_dyn_i+1:offset_dyn_i+docp.dim_OCP_x]
        if docp.is_lagrange
            c[offset+docp.dim_NLP_x] = ki[docp.dim_NLP_x] - work[offset_dyn_i+docp.dim_NLP_x]
        end
        offset += docp.dim_NLP_x
    end
    
    # 2. path constraints (control, state, mixed)
    if docp.dim_u_cons > 0
        docp.control_constraints[2]((@view c[offset+1:offset+docp.dim_u_cons]),ti, ui, v)
    end
    if docp.dim_x_cons > 0 
        docp.state_constraints[2]((@view c[offset+docp.dim_u_cons+1:offset+docp.dim_u_cons+docp.dim_x_cons]),ti, xi, v)
    end
    if docp.dim_xu_cons > 0 
        docp.mixed_constraints[2]((@view c[offset+docp.dim_u_cons+docp.dim_x_cons+1:offset+docp.dim_u_cons+docp.dim_x_cons+docp.dim_xu_cons]), ti, xi, ui, v)
    end

end






#######
# +++ ditch work array and put everything in setconstraintblock.
# compare bench vs previous midpoint

struct Midpoint2 <: Discretization

    stage::Int
    _step_pathcons_block::Int
    control_disc::Symbol
    info::String

    # constructor
    function Midpoint2(dim_NLP_steps, dim_NLP_x, dim_NLP_u, dim_NLP_v, dim_u_cons, dim_x_cons, dim_xu_cons, dim_boundary_cons, dim_v_cons)

        stage = 1

        # NLP variables size (state, control, variable, stage)
        dim_NLP_variables = (dim_NLP_steps + 1) * dim_NLP_x + dim_NLP_steps * dim_NLP_u + dim_NLP_v + dim_NLP_steps * dim_NLP_x * stage
        
        # Path constraints (control, state, mixed) 
        step_pathcons_block = dim_u_cons + dim_x_cons + dim_xu_cons

        # NLP constraints size (dynamics, stage, path, boundary, variable)
        dim_NLP_constraints = dim_NLP_steps * (dim_NLP_x + (dim_NLP_x * stage) + step_pathcons_block) + step_pathcons_block + dim_boundary_cons + dim_v_cons

        disc = new(stage, step_pathcons_block, :step, "Implicit Midpoint aka Gauss-Legendre collocation for s=1, 2nd order, symplectic")

        return disc, dim_NLP_variables, dim_NLP_constraints
    end
end


function setWorkArray(docp::DOCP{Midpoint2}, xu, time_grid, v)
  
end

function setStepConstraints!(docp::DOCP{Midpoint2}, c, xu, v, time_grid, i, work)

end
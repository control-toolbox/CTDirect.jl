#= Functions for trapeze discretization scheme
Internal layout for NLP variables: 
[X_1,U_1, .., X_N+1,U_N+1, V]
=#

# NB. could be defined as a generic IRK
struct Trapeze <: Discretization

    info::String
    _step_pathcons_block::Int
    _state_stage_eqs_block::Int

    # constructor
    function Trapeze(dim_NLP_steps, dim_NLP_x, dim_NLP_u, dim_NLP_v, dim_u_cons, dim_x_cons, dim_xu_cons, dim_boundary_cons, dim_v_cons)

        # NLP variables size ([state, control]_1..N+1, variable)
        dim_NLP_variables = (dim_NLP_steps + 1) * (dim_NLP_x + dim_NLP_u) + dim_NLP_v

        # Path constraints (control, state, mixed)
        state_stage_eqs_block = dim_NLP_x
        step_pathcons_block = dim_u_cons + dim_x_cons + dim_xu_cons

        # NLP constraints size ([dynamics, stage, path]_1..N, final path, boundary, variable)
        dim_NLP_constraints = dim_NLP_steps * (state_stage_eqs_block + step_pathcons_block) + step_pathcons_block + dim_boundary_cons + dim_v_cons

        disc = new("Implicit Trapeze aka Crank-Nicolson, 2nd order, A-stable", step_pathcons_block, state_stage_eqs_block)

        return disc, dim_NLP_variables, dim_NLP_constraints
    end
end

# not only for Trapeze, but may be redefined if needed
function get_OCP_variable(xu, docp::DOCP{<: Discretization, <: ScalVect, <: ScalVect, ScalVariable})
    return xu[docp.dim_NLP_variables]
end
function get_OCP_variable(xu, docp::DOCP{<: Discretization, <: ScalVect, <: ScalVect, VectVariable})
    return @view xu[(docp.dim_NLP_variables - docp.dim_NLP_v + 1):docp.dim_NLP_variables]
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
    # use work array to store all dynamics + lagrange costs
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
    offset = (i-1)*(docp.dim_NLP_x + docp.discretization._step_pathcons_block)

    # 0. variables
    ti = time_grid[i]
    xi = get_OCP_state_at_time_step(xu, docp, i)
    ui = get_OCP_control_at_time_step(xu, docp, i)

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

    # 2. path constraints (control, state, mixed)
    setPathConstraints!(docp, c, ti, xi, ui, v, offset)

end

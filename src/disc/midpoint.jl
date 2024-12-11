#= Functions for implicit midpoint discretization scheme
Internal layout for NLP variables: 
[X_1,U_1, .., X_N,U_N, X_N+1, V]
with the convention u([t_i,t_i+1[) = U_i and u(tf) = U_N
NB. stage variables are removed via the simplification x_s = (x_i + x_i+1) / 2
=#

struct Midpoint <: Discretization

    info::String
    _step_variables_block::Int
    _state_stage_eqs_block::Int
    _step_pathcons_block::Int
    _kvars::Bool

    # constructor
    function Midpoint(dim_NLP_steps, dim_NLP_x, dim_NLP_u, dim_NLP_v, dim_u_cons, dim_x_cons, dim_xu_cons, dim_boundary_cons, dim_v_cons; kvars=false)

        # aux variables
        if kvars
            step_variables_block = dim_NLP_x * 2 + dim_NLP_u
            state_stage_eqs_block = dim_NLP_x * 2
        else
            step_variables_block = dim_NLP_x + dim_NLP_u
            state_stage_eqs_block = dim_NLP_x
        end

        # NLP variables size ([state, control]_1..N, final state, variable)
        dim_NLP_variables = dim_NLP_steps * step_variables_block + dim_NLP_x + dim_NLP_v
        
        # Path constraints (control, state, mixed)
        step_pathcons_block = dim_u_cons + dim_x_cons + dim_xu_cons

        # NLP constraints size ([dynamics, path]_1..N, final path, boundary, variable)
        dim_NLP_constraints = dim_NLP_steps * (state_stage_eqs_block + step_pathcons_block) + step_pathcons_block + dim_boundary_cons + dim_v_cons

        disc = new("Implicit Midpoint aka Gauss-Legendre collocation for s=1, 2nd order, symplectic", step_variables_block, state_stage_eqs_block, step_pathcons_block, kvars)

        return disc, dim_NLP_variables, dim_NLP_constraints
    end
end


"""
$(TYPEDSIGNATURES)

Retrieve state variables at given time step from the NLP variables.
Convention: 1 <= i <= dim_NLP_steps+1
Scalar / Vector output
"""
function get_OCP_state_at_time_step(xu, docp::DOCP{Midpoint, ScalVariable, <: ScalVect, <: ScalVect}, i)
    offset = (i-1) * docp.discretization._step_variables_block
    return xu[offset+1]
end
function get_OCP_state_at_time_step(xu, docp::DOCP{Midpoint, VectVariable, <: ScalVect, <: ScalVect}, i)
    offset = (i-1) * docp.discretization._step_variables_block
    return @view xu[(offset + 1):(offset + docp.dim_OCP_x)]
end
"""
$(TYPEDSIGNATURES)

Retrieve state variable for lagrange cost at given time step from the NLP variables.
Convention: 1 <= i <= dim_NLP_steps+1
"""
function get_lagrange_state_at_time_step(xu, docp::DOCP{Midpoint}, i)
    offset = (i-1) * docp.discretization._step_variables_block
    return xu[offset + docp.dim_NLP_x]
end
"""
$(TYPEDSIGNATURES)

Retrieve control variables at given time step from the NLP variables.
Convention: 1 <= i <= dim_NLP_steps(+1), with convention u(tf) = U_N
Scalar / Vector output
"""
function get_OCP_control_at_time_step(xu, docp::DOCP{Midpoint, <: ScalVect, ScalVariable, <: ScalVect}, i)
    # final time case
    (i == docp.dim_NLP_steps + 1) && (i = docp.dim_NLP_steps)
    offset = (i-1) * docp.discretization._step_variables_block + docp.dim_NLP_x
    return xu[offset+1]
end
function get_OCP_control_at_time_step(xu, docp::DOCP{Midpoint, <: ScalVect, VectVariable, <: ScalVect}, i)
    # final time case
    (i == docp.dim_NLP_steps + 1) && (i = docp.dim_NLP_steps)
    offset = (i-1) * docp.discretization._step_variables_block + docp.dim_NLP_x
    return @view xu[(offset + 1):(offset + docp.dim_NLP_u)]
end

"""
$(TYPEDSIGNATURES)

Retrieve stage variables at given time step/stage from the NLP variables.
Convention: 1 <= i <= dim_NLP_steps,	1 <= j <= s
Vector output
"""
function get_stagevars_at_time_step(xu, docp::DOCP{Midpoint}, i, j)
    offset = (i-1) * docp.discretization._step_variables_block + docp.dim_NLP_x + docp.dim_NLP_u + (j-1)*docp.dim_NLP_x
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
        offset = (i-1) * docp.discretization._step_variables_block
        xu[(offset + 1):(offset + docp.dim_OCP_x)] .= x_init
    end
end
"""
$(TYPEDSIGNATURES)

Set initial guess for control variables at given time step
Convention: 1 <= i <= dim_NLP_steps
"""
function set_control_at_time_step!(xu, u_init, docp::DOCP{Midpoint}, i)
    if i <= docp.dim_NLP_steps && !isnothing(u_init)
        offset = (i-1) * docp.discretization._step_variables_block + docp.dim_NLP_x
        xu[(offset + 1):(offset + docp.dim_NLP_u)] .= u_init
    end
end

"""
$(TYPEDSIGNATURES)

Set work array for all dynamics and lagrange cost evaluations
"""
function setWorkArray(docp::DOCP{Midpoint}, xu, time_grid, v)
    
    # NB. recheck performance vs using stage variables again
    work = similar(xu, docp.dim_OCP_x)
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
    ui = get_OCP_control_at_time_step(xu, docp, i)

    # 1. state equation
    if i <= docp.dim_NLP_steps
        # more variables
        tip1 = time_grid[i+1]
        xip1 = get_OCP_state_at_time_step(xu, docp, i+1)
        hi = tip1 - ti
        ts = 0.5 * (ti + tip1)
        @views @. work[1:docp.dim_OCP_x] = 0.5 * (xi + xip1)
        if docp.dim_OCP_x == 1
            xs = work[1]
        else
            xs = work
        end
       
        if docp.discretization._kvars
            # state equation: midpoint rule
            ki = get_stagevars_at_time_step(xu, docp, i, 1)
            @views @. c[offset+1:offset+docp.dim_OCP_x] = xip1 - (xi + hi * ki[1:docp.dim_OCP_x])
            if docp.is_lagrange
            c[offset+docp.dim_NLP_x] = get_lagrange_state_at_time_step(xu, docp, i+1) - (get_lagrange_state_at_time_step(xu, docp, i) + hi * ki[docp.dim_NLP_x])
            end
            offset += docp.dim_NLP_x

            # stage equation at mid step k_i = f(x_s)
            docp.ocp.dynamics((@view c[offset+1:offset+docp.dim_OCP_x]), ts, xs, ui, v)
            if docp.is_lagrange
                docp.ocp.lagrange((@view c[offset+docp.dim_NLP_x:offset+docp.dim_NLP_x]), ts, xs, ui, v)
            end            
            @views @. c[offset+1:offset+docp.dim_NLP_x] = ki - c[offset+1:offset+docp.dim_NLP_x]
            offset += docp.dim_NLP_x

        else
            # No stage variables, compute dynamics directly for state equation
            # OCP dynamics
            docp.ocp.dynamics((@view c[offset+1:offset+docp.dim_OCP_x]), ts, xs, ui, v)
            # lagrange cost
            if docp.is_lagrange
                docp.ocp.lagrange((@view c[offset+docp.dim_NLP_x:offset+docp.dim_NLP_x]), ts, xs, ui, v)
            end

            # midpoint rule
            @views @. c[offset+1:offset+docp.dim_OCP_x] = xip1 - (xi + hi * c[offset+1:offset+docp.dim_OCP_x])
            if docp.is_lagrange
            c[offset+docp.dim_NLP_x] = get_lagrange_state_at_time_step(xu, docp, i+1) - (get_lagrange_state_at_time_step(xu, docp, i) + hi * c[offset+docp.dim_NLP_x])
            end
            offset += docp.dim_NLP_x
        end

    end
   
    # 2. path constraints
    setPathConstraints!(docp, c, ti, xi, ui, v, offset)
    
end



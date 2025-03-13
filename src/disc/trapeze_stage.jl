#= Functions for trapeze discretization scheme
Internal layout for NLP variables: 
[X_1,U_1,K_1, ... , X_N+1,U_N+1,K_N+1, V]
with 'stage' variables K_i = f(T_i, X_i, U_i, V)
No work array
=#

struct Trapeze_stage <: Discretization

    info::String
    _step_variables_block::Int
    _state_stage_eqs_block::Int
    _step_pathcons_block::Int

    # constructor
    function Trapeze_stage(dim_NLP_steps, dim_NLP_x, dim_NLP_u, dim_NLP_v, dim_u_cons, dim_x_cons, dim_xu_cons, dim_boundary_cons, dim_v_cons)

        # aux variables
        step_variables_block = dim_NLP_x * 2 + dim_NLP_u
        state_stage_eqs_block = dim_NLP_x * 2
        step_pathcons_block = dim_u_cons + dim_x_cons + dim_xu_cons

        # NLP variables size ([state, control]_1..N+1, variable)
        dim_NLP_variables = (dim_NLP_steps + 1) * step_variables_block + dim_NLP_v

        # NLP constraints size ([dynamics, stage, path]_1..N, final stage and path, boundary, variable)
        dim_NLP_constraints = dim_NLP_steps * (state_stage_eqs_block + step_pathcons_block) + dim_NLP_x + step_pathcons_block + dim_boundary_cons + dim_v_cons

        disc = new("Implicit Trapeze aka Crank-Nicolson, 2nd order, A-stable", step_variables_block, state_stage_eqs_block, step_pathcons_block)

        return disc, dim_NLP_variables, dim_NLP_constraints
    end
end


"""
$(TYPEDSIGNATURES)

Retrieve control variables at given time step from the NLP variables.
Convention: 1 <= i <= dim_NLP_steps+1
Scalar / Vector output
"""
function get_OCP_control_at_time_step(xu, docp::DOCP{Trapeze_stage, <: ScalVect, ScalVariable, <: ScalVect}, i)
    offset = (i-1) * docp.discretization._step_variables_block + docp.dim_NLP_x
    return xu[offset+1]
end
function get_OCP_control_at_time_step(xu, docp::DOCP{Trapeze_stage, <: ScalVect, VectVariable, <: ScalVect}, i)
    offset = (i-1) * docp.discretization._step_variables_block + docp.dim_NLP_x
    return @view xu[(offset + 1):(offset + docp.dim_NLP_u)]
end

function get_stagevars_at_time_step(xu, docp::DOCP{Trapeze_stage}, i)
    offset = (i-1) * docp.discretization._step_variables_block + docp.dim_NLP_x + docp.dim_NLP_u
    return @view xu[(offset + 1):(offset + docp.dim_NLP_x)]
end


"""
$(TYPEDSIGNATURES)

Set initial guess for control variables at given time step
Convention: 1 <= i <= dim_NLP_steps+1
"""
function set_control_at_time_step!(xu, u_init, docp::DOCP{Trapeze_stage}, i)
    if !isnothing(u_init)
        offset = (i-1) * docp.discretization._step_variables_block + docp.dim_NLP_x
        xu[(offset + 1):(offset + docp.dim_NLP_u)] .= u_init
    end
end


"""
$(TYPEDSIGNATURES)

Set the constraints corresponding to the state equation
Convention: 1 <= i <= dim_NLP_steps+1
"""
function setStepConstraints!(docp::DOCP{Trapeze_stage}, c, xu, v, time_grid, i, work)

    # offset for previous steps
    offset = (i-1)*(docp.discretization._state_stage_eqs_block + docp.discretization._step_pathcons_block)

    # 0. variables
    ti = time_grid[i]
    xi = get_OCP_state_at_time_step(xu, docp, i)
    fi = get_stagevars_at_time_step(xu, docp, i)

    # 1.1 state equation
    if i <= docp.dim_NLP_steps
        # more variables
        tip1 = time_grid[i+1]
        xip1 = get_OCP_state_at_time_step(xu, docp, i+1)
        fip1 = get_stagevars_at_time_step(xu, docp, i+1)
        half_hi = 0.5 * (tip1 - ti)

        # trapeze rule
        @views @. c[offset+1:offset+docp.dim_OCP_x] = xip1 - (xi + half_hi * (fi[1:docp.dim_OCP_x] + fip1[1:docp.dim_OCP_x]))
        if docp.is_lagrange
            c[offset+docp.dim_NLP_x] = get_lagrange_state_at_time_step(xu, docp, i+1) - (get_lagrange_state_at_time_step(xu, docp, i) + half_hi * (fi[docp.dim_NLP_x] + fip1[docp.dim_NLP_x]))
        end
        offset += docp.dim_NLP_x
    end

    # 1.2 stage equations k_i = f(t_i, x_i, u_i, v) including at tf
    ui = get_OCP_control_at_time_step(xu, docp, i)
    docp.ocp.dynamics((@view c[offset+1:offset+docp.dim_OCP_x]), ti, xi, ui, v)
    if docp.is_lagrange
        docp.ocp.lagrange((@view c[offset+docp.dim_NLP_x:offset+docp.dim_NLP_x]), ti, xi, ui, v)
    end            
    @views @. c[offset+1:offset+docp.dim_NLP_x] = fi - c[offset+1:offset+docp.dim_NLP_x]
    offset += docp.dim_NLP_x

    # 2. path constraints
    setPathConstraints!(docp, c, ti, xi, ui, v, offset)

end

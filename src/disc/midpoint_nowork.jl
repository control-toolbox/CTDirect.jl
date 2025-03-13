#= Functions for implicit midpoint discretization scheme
Internal layout for NLP variables: 
[X_1,U_1, .., X_N,U_N, X_N+1, V]
with the convention u([t_i,t_i+1[) = U_i and u(tf) = U_N
NB. This version give identical performance to the one with work array
=#

struct Midpoint_nowork <: Discretization

    info::String
    _step_variables_block::Int
    _state_stage_eqs_block::Int
    _step_pathcons_block::Int

    # constructor
    function Midpoint_nowork(dim_NLP_steps, dim_NLP_x, dim_NLP_u, dim_NLP_v, dim_u_cons, dim_x_cons, dim_xu_cons, dim_boundary_cons, dim_v_cons)

        # aux variables
        step_variables_block = dim_NLP_x + dim_NLP_u
        state_stage_eqs_block = dim_NLP_x
        step_pathcons_block = dim_u_cons + dim_x_cons + dim_xu_cons

        # NLP variables size ([state, control, stage]_1..N, final state, variable)
        dim_NLP_variables = dim_NLP_steps * step_variables_block + dim_NLP_x + dim_NLP_v
        
        # NLP constraints size ([dynamics, path]_1..N, final path, boundary, variable)
        dim_NLP_constraints = dim_NLP_steps * (state_stage_eqs_block + step_pathcons_block) + step_pathcons_block + dim_boundary_cons + dim_v_cons

        disc = new("Implicit Midpoint aka Gauss-Legendre collocation for s=1, 2nd order, symplectic", step_variables_block, state_stage_eqs_block, step_pathcons_block)

        return disc, dim_NLP_variables, dim_NLP_constraints
    end
end


"""
$(TYPEDSIGNATURES)

Set the constraints corresponding to the state equation
Convention: 1 <= i <= dim_NLP_steps+1
"""
function setStepConstraints!(docp::DOCP{Midpoint_nowork}, c, xu, v, time_grid, i, work)
    
    # offset for previous steps
    offset = (i-1)*(docp.discretization._state_stage_eqs_block + docp.discretization._step_pathcons_block)

    # 0. variables
    ti = time_grid[i]
    ui = get_OCP_control_at_time_step(xu, docp, i)
    xi = get_OCP_state_at_time_step(xu, docp, i)

    # 1. state equation
    if i <= docp.dim_NLP_steps
        # more variables
        tip1 = time_grid[i+1]
        xip1 = get_OCP_state_at_time_step(xu, docp, i+1)
        hi = tip1 - ti
        offset_dyn_i = (i-1)*docp.dim_NLP_x
       
        # state equation: midpoint rule
        ts = 0.5 * (ti + tip1)
        xs = 0.5 * (xi + xip1)
        docp.ocp.dynamics((@view c[offset+1:offset+docp.dim_OCP_x]), ts, xs, ui, v)
        @views @. c[offset+1:offset+docp.dim_OCP_x] = xip1 - (xi + hi * c[offset+1:offset+docp.dim_OCP_x])

        if docp.is_lagrange
            docp.ocp.lagrange((@view c[offset+docp.dim_NLP_x:offset+docp.dim_NLP_x]), ts, xs, ui, v)
            c[offset+docp.dim_NLP_x] = get_lagrange_state_at_time_step(xu, docp, i+1) - (get_lagrange_state_at_time_step(xu, docp, i) + hi * c[offset+docp.dim_NLP_x])
        end
        offset += docp.dim_NLP_x

    end
   
    # 2. path constraints
    if docp.discretization._step_pathcons_block > 0
        setPathConstraints!(docp, c, ti, xi, ui, v, offset)
    end
    
end



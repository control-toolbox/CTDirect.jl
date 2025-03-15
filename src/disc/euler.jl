#= Functions for explicit and implicit euler discretization scheme
Internal layout for NLP variables: 
[X_1,U_1,K_1 .., X_N,U_N,K_N, X_N+1, V]
with the convention 
- Explicit Euler: u([t_i,t_i+1[) = U_i and u(tf) = U_N
- Implicit Euler: u(]t_i,t_i+1]) = U_i and u(t0) = U_1
Note that both the explicit and implicit versions therefore use the same variables layout.
=#

struct Euler <: Discretization

    info::String
    _step_variables_block::Int
    _state_stage_eqs_block::Int
    _step_pathcons_block::Int
    _final_control::Bool
    _explicit::Bool

    # constructor
    function Euler(dim_NLP_steps, dim_NLP_x, dim_NLP_u, dim_NLP_v, dim_path_cons, dim_boundary_cons; explicit=true)

        # aux variables
        step_variables_block = dim_NLP_x + dim_NLP_u
        state_stage_eqs_block = dim_NLP_x
        step_pathcons_block = dim_path_cons

        # NLP variables size ([state, control]_1..N, final state, variable)
        dim_NLP_variables = dim_NLP_steps * step_variables_block + dim_NLP_x + dim_NLP_v
        
        # NLP constraints size ([dynamics, path]_1..N, final path, boundary, variable)
        dim_NLP_constraints = dim_NLP_steps * (state_stage_eqs_block + step_pathcons_block) + step_pathcons_block + dim_boundary_cons

        if explicit 
            info = "Euler (explicit), 1st order"
        else
            info = "Euler (implicit), 1st order"
        end
        disc = new(info, step_variables_block, state_stage_eqs_block, step_pathcons_block, false, explicit)

        return disc, dim_NLP_variables, dim_NLP_constraints
    end
end


"""
$(TYPEDSIGNATURES)

Retrieve control variables at given time step from the NLP variables.
Convention: see above for acplicit / implicit versions
Vector output
"""
function get_OCP_control_at_time_step(xu, docp::DOCP{Euler}, i)
    if docp.discretization._explicit
        # final time case
        (i == docp.dim_NLP_steps + 1) && (i = docp.dim_NLP_steps)
        offset = (i-1) * docp.discretization._step_variables_block + docp.dim_NLP_x
    else 
        # initial time case
        (i == 1) && (i = 2)
        offset = (i-2) * docp.discretization._step_variables_block + docp.dim_NLP_x
    end
    return @view xu[(offset + 1):(offset + docp.dim_NLP_u)]
end


"""
$(TYPEDSIGNATURES)

Set work array for all dynamics and lagrange cost evaluations
"""
function setWorkArray(docp::DOCP{Euler}, xu, time_grid, v)

    work = similar(xu, docp.dim_NLP_x * docp.dim_NLP_steps)

    # loop over time steps
    for i = 1:docp.dim_NLP_steps
        offset = (i-1) * docp.dim_NLP_x

        # get variables at t_i or t_i+1
        if docp.discretization._explicit
            index = i
        else
            index = i+1
        end
        t = time_grid[index]
        x = get_OCP_state_at_time_step(xu, docp, index)
        u = get_OCP_control_at_time_step(xu, docp, index)

        # OCP dynamics
        CTModels.dynamics(docp.ocp)((@view work[offset+1:offset+docp.dim_OCP_x]), t, x, u, v)
        # lagrange cost
        if docp.has_lagrange
            work[offset+docp.dim_NLP_x] = CTModels.lagrange(docp.ocp)(t, x, u, v)
        end   
    end
    return work
end


"""
$(TYPEDSIGNATURES)

Set the constraints corresponding to the state equation
Convention: 1 <= i <= dim_NLP_steps+1
"""
function setStepConstraints!(docp::DOCP{Euler}, c, xu, v, time_grid, i, work)
    
    # offset for previous steps
    offset = (i-1)*(docp.discretization._state_stage_eqs_block + docp.discretization._step_pathcons_block)

    # 0. variables
    ti = time_grid[i]
    xi = get_OCP_state_at_time_step(xu, docp, i)

    # 1. state equation
    if i <= docp.dim_NLP_steps
        # more variables
        tip1 = time_grid[i+1]
        xip1 = get_OCP_state_at_time_step(xu, docp, i+1)
        hi = tip1 - ti
        offset_dyn_i = (i-1)*docp.dim_NLP_x

        # state equation: euler rule
        @views @. c[offset+1:offset+docp.dim_OCP_x] = xip1 - (xi + hi * work[offset_dyn_i+1:offset_dyn_i+docp.dim_OCP_x])
        if docp.has_lagrange
            c[offset+docp.dim_NLP_x] = get_lagrange_state_at_time_step(xu, docp, i+1) - (get_lagrange_state_at_time_step(xu, docp, i) + hi * work[offset_dyn_i+docp.dim_NLP_x])
        end
        offset += docp.dim_NLP_x

    end
   
    # 2. path constraints
    if docp.discretization._step_pathcons_block > 0
        ui = get_OCP_control_at_time_step(xu, docp, i)
        CTModels.path_constraints_nl(docp.ocp)[2]((@view c[offset+1:offset+docp.dim_path_cons]), ti, xi, ui, v)
        offset += docp.dim_path_cons
    end
    
end

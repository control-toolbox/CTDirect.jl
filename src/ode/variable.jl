# variable step ode solvers for direct shooting
import DifferentialEquations as DE

struct VariableStepODE <: Scheme
    info::String
    _step_variables_block::Int
    _state_stage_eqs_block::Int
    _step_pathcons_block::Int
    _final_control::Bool

    function VariableStepODE(dims::DOCPdims, time::DOCPtime)
        
        step_variables_block = dims.NLP_x + dims.NLP_u * time.control_steps
        state_stage_eqs_block = dims.NLP_x
        step_pathcons_block = dims.path_cons

        # NLP variables size ([state, controls]_1..N, final state, variable)
        dim_NLP_variables = time.steps * step_variables_block + dims.NLP_x + dims.NLP_v

        # NLP constraints size ([state eq, path]_1..N, final_path, boundary)
        dim_NLP_constraints = time.steps * (dims.NLP_x + step_pathcons_block) + 
        step_pathcons_block + dims.boundary_cons
        
        disc = new(
            "Variable step ODE solver",
            step_variables_block,
            state_stage_eqs_block,
            step_pathcons_block,
            false,
        )

        return disc, dim_NLP_variables, dim_NLP_constraints
    end
end

# +++ move later to state parametrization
# get state at any given time (use piecewise constant and later linear parametrization)
function get_OCP_state_at_time(xu, docp::DOCP{VariableStepODE}, t; mode=:constant)

    # retrieve actual time grid (needed for variable time case)
    time_grid = get_time_grid(xu, docp)

    # find index i such that t is in [t_i, t_i+1]
    i = searchsortedlast(t, time_grid)

    return get_OCP_state_at_time_step(xu, docp, i)

end

# +++ move later to control parametrization
# get control at any given time (use piecewise constant and later linear parametrization)
function get_OCP_control_at_time(xu, docp::DOCP{VariableStepODE}, t; mode=:constant)
    
    # retrieve actual time grid (needed for variable time case)
    time_grid = get_time_grid(xu, docp)

    # find index i such that t is in [t_i, t_i+1]
    i = searchsortedlast(t, time_grid)

    if (i == docp.time.steps + 1) || (docp.time.control_steps == 1)
        return get_OCP_control_at_time_step(xu, docp, i)
    else
        # find index j for piecewise constant control over [t_i, t_i+1]
        # nb control getter will later use single index instead of i,j pair
        control_step = (time_grid[i+1] - time_grid[i]) / docp.time.control_steps
        j = floor((t - t[i]) / control_step) + 1
        (j > docp.time.control_steps) && (j = docp.time.control_steps)
        return get_OCP_control_at_time_step(xu, docp, i, j)
    end

end

function setWorkArray(docp::DOCP{VariableStepODE}, xu, time_grid, v)
    return
end

# NB. here we will use approximate values for the state 
# this is not equivalent to reformulating lagrange as mayer if lagrange depends on x !
function runningCost(docp::DOCP{VariableStepODE}, xu, v, time_grid)

    function lagrange_closure(value, v, t)
        x = get_OCP_state_at_time(xu, docp, t)
        u = get_OCP_control_at_time(xu, docp, t)
        return CTModels.lagrange(ocp_model(docp))(t,x,u,v)
    end

    prob = DE.ODEProblem(lagrange_closure, 0.0, (time_grid[1], time_grid[end]))
    sol = DE.solve(prob)

    return sol
end

# +++ we could also solve once over [t0,tf] and use dense output for state constraints ?
# this would be a different version of __constraints though...
function stepStateConstraints!(docp::DOCP{VariableStepODE}, c, xu, v, time_grid, i, work)

    # variables
    disc = disc_model(docp)
    ti = time_grid[i]
    xi = get_OCP_state_at_time_step(xu, docp, i)
    tip1 = time_grid[i+1]
    xip1 = get_OCP_state_at_time_step(xu, docp, i+1)
    offset_c = (i-1)*(disc._state_stage_eqs_block + disc._step_pathcons_block)

    function dynamics_closure(x, v, t)
        u = get_OCP_control_at_time(xu, docp, t)
        return CTModels.dynamics(ocp_model(docp))(t,x,u,v)
    end

    # integrate dynamics from t_i to t_i+1     +++ later try to use c directly for x_next ?
    prob = DE.ODEProblem(dynamics_closure, xi, (ti, tip1))
    x_next = DE.solve(prob)

    # set constraint
    @views @. c[(offset_c + 1):(offset_c + docp.dims.NLP_x)] = xip1 - x_next
    return
end



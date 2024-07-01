"""
$(TYPEDSIGNATURES)

Retrieve optimization variables from the NLP variables
"""
function get_variable(xu, docp)
    if is_variable_dependent(docp.ocp)
        if docp.dim_NLP_v == 1
            return xu[end]
        else
            return xu[end-docp.dim_NLP_v+1:end]
        end
    else
        return Float64[]
    end
end


"""
$(TYPEDSIGNATURES)

Retrieve state variables at given time step from the NLP variables
"""
function get_state_at_time_step(xu, docp, i::Int64)
    nx = docp.dim_NLP_x
    n = docp.ocp.state_dimension
    N = docp.dim_NLP_steps
    @assert i <= N "trying to get x(t_i) for i > N"
    if n == 1
        return xu[i*nx + 1]
    else
        return xu[i*nx + 1 : i*nx + n]
    end
end


"""
$(TYPEDSIGNATURES)

Retrieve the additional state variable corresponding to the lagrange (running) cost at given time step from the NLP variables
"""
function get_lagrange_cost_at_time_step(xu, docp, i)
    nx = docp.dim_NLP_x
    N = docp.dim_NLP_steps
    @assert i <= N "trying to get lagrange cost at t_i for i > N"
    return xu[(i+1)*nx]
end

# internal vector version
function vget_state_at_time_step(xu, docp, i)
    nx = docp.dim_NLP_x
    N = docp.dim_NLP_steps
    @assert i <= N "trying to get x(t_i) for i > N"
    return xu[i*nx + 1 : (i+1)*nx]
end


"""
$(TYPEDSIGNATURES)

Retrieve control variables at given time step from the NLP variables
"""
function get_control_at_time_step(xu, docp, i)
    nx = docp.dim_NLP_x
    m = docp.dim_NLP_u
    N = docp.dim_NLP_steps
    @assert i <= N "trying to get u(t_i) for i > N"
    if m == 1
        return xu[(N+1)*nx + i*m + 1]
    else
        return xu[(N+1)*nx + i*m + 1 : (N+1)*nx + (i+1)*m]
    end
end

# internal vector version
function vget_control_at_time_step(xu, docp, i)
    nx = docp.dim_NLP_x
    m = docp.dim_NLP_u
    N = docp.dim_NLP_steps
    @assert i <= N "trying to get u(t_i) for i > N"
    return xu[(N+1)*nx + i*m + 1 : (N+1)*nx + (i+1)*m]
end


"""
$(TYPEDSIGNATURES)

Retrieve initial time for OCP (may be fixed or variable)
"""
function get_initial_time(xu, docp)
    if has_free_initial_time(docp.ocp)
        v = get_variable(xu, docp)
        return v[docp.ocp.initial_time]
    else
        return docp.ocp.initial_time
    end
end


"""
$(TYPEDSIGNATURES)

Retrieve final time for OCP (may be fixed or variable)
"""
function get_final_time(xu, docp)
    if has_free_final_time(docp.ocp)
        v = get_variable(xu, docp)
        return v[docp.ocp.final_time]
    else
        return docp.ocp.final_time
    end
end


"""
$(TYPEDSIGNATURES)

Get actual (un-normalized) time value
"""
function get_unnormalized_time(xu, docp, t_normalized)
    t0 = get_initial_time(xu, docp)
    tf = get_final_time(xu, docp)    
    return t0 + t_normalized * (tf - t0)
end


"""
$(TYPEDSIGNATURES)

Get actual (un-normalized) time at give time step
"""
function get_time_at_time_step(xu, docp, i)
    N = docp.dim_NLP_steps
    @assert i <= N "trying to get t_i for i > N"
    t_normalized = docp.NLP_normalized_time_grid[i+1]
    return get_unnormalized_time(xu, docp, t_normalized)
end


"""
$(TYPEDSIGNATURES)

Set state variables at given time step in the NLP variables (for initial guess)
"""
function set_state_at_time_step!(xu, x_init, docp, i)
    nx = docp.dim_NLP_x
    n = docp.ocp.state_dimension
    N = docp.dim_NLP_steps
    @assert i <= N "trying to set init for x(t_i) with i > N"
    # NB. only set first the actual state variables from the OCP (not the possible additional state for lagrange cost)
    if n == 1
        xu[i*n + 1] = x_init[]
    else
        xu[i*nx + 1 : i*nx + n] = x_init
    end
end


"""
$(TYPEDSIGNATURES)

Set control variables at given time step in the NLP variables (for initial guess)
"""
function set_control_at_time_step!(xu, u_init, docp, i)
    nx = docp.dim_NLP_x
    m = docp.dim_NLP_u
    N = docp.dim_NLP_steps
    @assert i <= N "trying to set init for u(t_i) with i > N"
    offset = (N+1)*nx
    if m == 1
        xu[offset + i*m + 1] = u_init[]
    else        
        xu[offset + i*m + 1 : offset + i*m + m] = u_init
    end
end


"""
$(TYPEDSIGNATURES)

Set optimization variables in the NLP variables (for initial guess)
"""
function set_variable!(xu, v_init, docp)
    if docp.dim_NLP_v == 1
        xu[end] = v_init[]
    else
        xu[end-docp.dim_NLP_v+1 : end] = v_init
    end
end


# local init
function checkDim(actual_dim, target_dim)
    if !isnothing(target_dim) && actual_dim != target_dim
        error("Init dimension mismatch: got ",actual_dim," instead of ",target_dim )
    end
end

function isaVectVect(data)
    return (data isa Vector) && (data[1] isa ctVector)
end

function formatData(data)
    if data isa Matrix
        return matrix2vec(data,1)
    else
        return data
    end
end

function formatTimeGrid(time)
    if isnothing(time)
        return nothing
    elseif time isa ctVector
        return time
    else
        return vec(time)
    end
end

function buildFunctionalInit(data::Nothing, time, dim)
    # fallback to method-dependent default initialization
    return t-> nothing
end

function buildFunctionalInit(data::Function, time, dim)
    # functional initialization
    checkDim(length(data(0)),dim)
    return t -> data(t)
end

function buildFunctionalInit(data::ctVector, time, dim)
    if !isnothing(time) && (length(data) == length(time))
        # interpolation vs time, dim 1 case
        itp = ctinterpolate(time, data)
        return t -> itp(t)
    else
        # constant initialization
        checkDim(length(data), dim)
        return t -> data
    end
end

function buildFunctionalInit(data, time, dim)
    if isaVectVect(data)
        # interpolation vs time, general case
        itp = ctinterpolate(time, data)
        checkDim(length(itp(0)), dim)
        return t -> itp(t)
    else
        error("Unrecognized initialization argument: ",typeof(data))
    end
end

function buildVectorInit(data, dim)
    if isnothing(data)
        return data
    else
        checkDim(length(data),dim)
        return data
    end
end

mutable struct _OptimalControlInit
   
    state_init::Function
    control_init::Function
    variable_init::Union{Nothing, ctVector}
    #costate_init::Function
    #multipliers_init::Union{Nothing, ctVector}

    # base constructor with explicit arguments
    function _OptimalControlInit(; state=nothing, control=nothing, variable=nothing, time=nothing, state_dim=nothing, control_dim=nothing, variable_dim=nothing)
        
        init = new()
    
        # some matrix / vector conversions
        time = formatTimeGrid(time)
        state = formatData(state)
        control = formatData(control)
        
        # set initialization for x, u, v
        init.state_init = buildFunctionalInit(state, time, state_dim)
        init.control_init = buildFunctionalInit(control, time, control_dim)
        init.variable_init = buildVectorInit(variable, variable_dim)

        return init

    end

    # version with arguments as named tuple or dict
    function _OptimalControlInit(init_data; state_dim=nothing, control_dim=nothing, variable_dim=nothing)

        # trivial case: default init
        x_init = nothing
        u_init = nothing
        v_init = nothing
        t_init = nothing
        x_dim = nothing
        u_dim = nothing
        v_dim = nothing

        # parse arguments and call base constructor
        if !isnothing(init_data)
            for key in keys(init_data)
                if key == :state
                    x_init = init_data[:state]
                elseif key == :control
                    u_init = init_data[:control]
                elseif key == :variable
                    v_init = init_data[:variable]
                elseif key == :time
                    t_init = init_data[:time]
                else
                    error("Unknown key in initialization data (allowed: state, control, variable, time, state_dim, control_dim, variable_dim): ", key)
                end
            end
        end

        return _OptimalControlInit(state=x_init, control=u_init, variable=v_init, time=t_init, state_dim=state_dim, control_dim=control_dim, variable_dim=variable_dim)
    
    end

    # warm start from solution
    function _OptimalControlInit(sol::OptimalControlSolution; unused_kwargs...)
        return _OptimalControlInit(state=sol.state, control=sol.control, variable=sol.variable, state_dim=sol.state_dimension, control_dim=sol.control_dimension, variable_dim=sol.variable_dimension)
    end

end


"""
$(TYPEDSIGNATURES)

Build initial guess for discretized problem
"""
function DOCP_initial_guess(docp,
    init::_OptimalControlInit=_OptimalControlInit())

    # default initialization
    # note: internal variables (lagrange cost, k_i for RK schemes) will keep these default values 
    xuv = 0.1 * ones(docp.dim_NLP_variables)

    # set variables if provided
    # (needed first in case of free times !)
    if !isnothing(init.variable_init)
        set_variable!(xuv, init.variable_init, docp)
    end

    # set state / control variables if provided
    for i in 0:docp.dim_NLP_steps
        ti = get_time_at_time_step(xuv, docp, i)
        if !isnothing(init.state_init(ti))
            set_state_at_time_step!(xuv, init.state_init(ti), docp, i)
        end
        if !isnothing(init.control_init(ti))
            set_control_at_time_step!(xuv, init.control_init(ti), docp, i)
        end
    end

    return xuv
end


#struct for interpolated ocp solution with only basic data types that can be exported as json
# +++todo:
# - add more fields from OptimalControlSolution
# - constructor to recreate OptimalControlSolution from this one
mutable struct OCP_Solution_discrete

    times
    #initial_time_name::Union{String, Nothing}=nothing
    #final_time_name::Union{String, Nothing}=nothing
    #time_name::Union{String, Nothing}=nothing
    control_dimension
    #control_components_names::Union{Vector{String}, Nothing}=nothing
    #control_name::Union{String, Nothing}=nothing
    control
    state_dimension
    #state_components_names::Union{Vector{String}, Nothing}=nothing
    #state_name::Union{String, Nothing}=nothing
    state
    variable_dimension
    #variable_components_names::Union{Vector{String}, Nothing}=nothing
    #variable_name::Union{String, Nothing}=nothing
    variable
    costate
    objective
    #iterations::Union{Nothing, Integer}=nothing
    #stopping::Union{Nothing, Symbol}=nothing # the stopping criterion
    #message::Union{Nothing, String}=nothing # the message corresponding to the stopping criterion
    #success::Union{Nothing, Bool}=nothing # whether or not the method has finished successfully: CN1, stagnation vs iterations max
    #infos::Dict{Symbol, Any}=Dict{Symbol, Any}()
    #OCP_Solution_discrete() = new() # for StructTypes / JSON
    
    function OCP_Solution_discrete(solution::OptimalControlSolution)
        solution_d = new()

        # raw copy +++ reuse the copy! in CTBase ?
        solution_d.objective = solution.objective
        solution_d.times = solution.times
        solution_d.state_dimension = solution.state_dimension
        solution_d.control_dimension = solution.control_dimension
        solution_d.variable_dimension = solution.variable_dimension
        solution_d.variable = solution.variable

        # interpolate functions into vectors
        solution_d.state = solution.state.(solution_d.times)
        solution_d.control = solution.control.(solution_d.times)
        solution_d.costate = solution.costate.(solution_d.times)
        return solution_d
    end
end

# placeholders (see CTDirectExt)
function save_OCP_solution end
function load_OCP_solution end
function export_OCP_solution end
function read_OCP_solution end

#=
"""
$(TYPEDSIGNATURES)
 
Parse interpolated OCP solution saved in JSON format (NOT IMPLEMENTED)
"""
function parse_JSON_solution(json_solution)
    println("parse_JSON_solution not implemented yet")
        #+++ parse JSON object containing interpolated solution
        #+++ then call raw constructor for OCP solution
end
=#

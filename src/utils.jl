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
    if docp.dim_OCP_x == 1
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

# NLP full state version
function get_NLP_state_at_time_step(xu, docp, i)
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
    if docp.dim_NLP_u == 1
        return xu[(N+1)*nx + i*m + 1]
    else
        return xu[(N+1)*nx + i*m + 1 : (N+1)*nx + (i+1)*m]
    end
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
    return get_unnormalized_time(xu, docp, docp.NLP_normalized_time_grid[i+1])
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
    xu[i*nx + 1 : i*nx + n] .= x_init
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
    xu[offset + i*m + 1 : offset + i*m + m] .= u_init
end


"""
$(TYPEDSIGNATURES)

Set optimization variables in the NLP variables (for initial guess)
"""
function set_variable!(xu, v_init, docp)
    xu[end-docp.dim_NLP_v+1 : end] .= v_init
end


"""
$(TYPEDSIGNATURES)

Build initial guess for discretized problem
"""
function DOCP_initial_guess(docp,
    init::OptimalControlInit=OptimalControlInit())

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
mutable struct OCPDiscreteSolution

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
    
    function OCPDiscreteSolution(solution::OptimalControlSolution)
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
function export_ocp_solution end
function import_ocp_solution end


"""
$(TYPEDSIGNATURES)

Retrieve optimization variables from the NLP variables
"""
function get_variable(xu, docp)
    if docp.has_variable
        if docp.variable_dimension == 1
            return xu[end]
        else
            return xu[end-docp.variable_dimension+1:end]
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
    nx = docp.dim_NLP_state
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
    nx = docp.dim_NLP_state
    N = docp.dim_NLP_steps
    @assert i <= N "trying to get lagrange cost at t_i for i > N"
    return xu[(i+1)*nx]
end

# internal vector version
function vget_state_at_time_step(xu, docp, i)
    nx = docp.dim_NLP_state
    N = docp.dim_NLP_steps
    @assert i <= N "trying to get x(t_i) for i > N"
    return xu[i*nx + 1 : (i+1)*nx]
end


"""
$(TYPEDSIGNATURES)

Retrieve control variables at given time step from the NLP variables
"""
function get_control_at_time_step(xu, docp, i)
    nx = docp.dim_NLP_state
    m = docp.ocp.control_dimension
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
    nx = docp.dim_NLP_state
    m = docp.ocp.control_dimension
    N = docp.dim_NLP_steps
    @assert i <= N "trying to get u(t_i) for i > N"
    return xu[(N+1)*nx + i*m + 1 : (N+1)*nx + (i+1)*m]
end


"""
$(TYPEDSIGNATURES)

Retrieve initial time for OCP (may be fixed or variable)
"""
function get_initial_time(xu, docp)
    if docp.has_free_initial_time
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
    if docp.has_free_final_time
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
    nx = docp.dim_NLP_state
    n = docp.ocp.state_dimension
    N = docp.dim_NLP_steps
    @assert i <= N "trying to set init for x(t_i) with i > N"
    # NB. only set first n components of state variable (ie the possible additional state for lagrange cost keeps the default init since it is not available from the OCP solution)
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
    nx = docp.dim_NLP_state
    m = docp.ocp.control_dimension
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
    if docp.variable_dimension == 1
        xu[end] = v_init[]
    else
        xu[end-docp.variable_dimension+1 : end] = v_init
    end
end


"""
$(TYPEDSIGNATURES)

Build initial guess for discretized problem
"""
function DOCP_initial_guess(docp, init::OCPInit=OCPInit())

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



#+++ to be moved to CTBase !
"""
$(TYPEDSIGNATURES)

Implement implicit call to OCPInit constructor if needed
"""
# +++ rename OCPInit as OptimalControlInit
# +++ replace this function by a trivial constructor method that returns its argument ?
function implicitInit(init_passed)
    # if argument passed for init is
    # - nothing or initialization data (vectors, functions, solution): call constructor
    # - already an OCPInit: just transmit
    if !(init_passed isa OCPInit)
        return OCPInit(init_passed)
    else
        return init_passed
    end
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

"""
$(TYPEDSIGNATURES)
  
Save OCP solution in JLD2/JSON format
"""
# +++ split as export/import for JSON
function save_OCP_solution(sol::OptimalControlSolution; filename_prefix="solution", format="JLD2")
    if format == "JLD2"
        save_object(filename_prefix * ".jld2", sol)
    elseif format == "JSON"
        open(filename_prefix * ".json", "w") do io
            JSON3.pretty(io, OCP_Solution_discrete(sol))
        end
    else
        println("ERROR: save_OCP_solution: format should be JLD2 or JSON, received ", format)
    end
    return nothing
end
    
"""
$(TYPEDSIGNATURES)
 
Load OCP solution in JLD2/JSON format
"""
function load_OCP_solution(filename_prefix="solution"; format="JLD2")
    if format == "JLD2"
        return load_object(filename_prefix * ".jld2")
    elseif format == "JSON"
        json_string = read(filename_prefix * ".json", String)
        return JSON3.read(json_string)
        #+++ parse JSON object containing interpolated solution
        #+++ then call raw constructor for OCP solution
    else
        println("ERROR: save_OCP_solution: format should be JLD2 or JSON, received ", format)
        return nothing
    end
end

#=
"""
$(TYPEDSIGNATURES)
 
Parse interpolated OCP solution saved in JSON format (NOT IMPLEMENTED)
"""
function parse_JSON_solution(json_solution)
    println("parse_JSON_solution not implemented yet")
end
=#

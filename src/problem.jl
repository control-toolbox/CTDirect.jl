# Internal layout for NLP variables: 
# [X0,U0, X1,U1, .., XN,UN,V]

"""
$(TYPEDSIGNATURES)

Struct for discretized optimal control problem DOCP

Contains:
- a copy of the original OCP
- a NLP formulation of the DOCP
- data required to link the two problems
"""
struct DOCP

    ## OCP
    ocp::OptimalControlModel
    control_constraints
    state_constraints
    mixed_constraints
    boundary_constraints
    variable_constraints
    control_box
    state_box
    variable_box

    # faster than calling ocp functions each time ?
    has_free_t0::Bool
    has_free_tf::Bool
    has_lagrange::Bool
    has_mayer::Bool
    has_variable::Bool
    has_maximization::Bool

    dim_x_box::Int
    dim_u_box::Int
    dim_v_box::Int
    dim_path_cons::Int
    dim_x_cons::Int
    dim_u_cons::Int
    dim_v_cons::Int
    dim_mixed_cons::Int
    dim_boundary_cons::Int

    ## NLP    
    dim_NLP_x::Int  # possible lagrange cost
    dim_NLP_u::Int
    dim_NLP_v::Int
    dim_OCP_x::Int  # original OCP state
    dim_NLP_steps::Int
    NLP_normalized_time_grid::Vector{Float64}
    dim_NLP_variables::Int
    dim_NLP_constraints::Int

    # lower and upper bounds for variables and constraints
    var_l::Vector{Float64}
    var_u::Vector{Float64}
    con_l::Vector{Float64}
    con_u::Vector{Float64}

    # constructor
    function DOCP(ocp::OptimalControlModel, grid_size::Integer, time_grid)       

        # time grid
        if time_grid == nothing
            NLP_normalized_time_grid = collect(LinRange(0, 1, grid_size+1))
            dim_NLP_steps = grid_size
        else
            # check strictly increasing
            if !issorted(time_grid,lt=<=)
                throw(ArgumentError("given time grid is not strictly increasing. Aborting..."))
                return nothing
            end
            # normalize input grid if needed
            if (time_grid[1] != 0) || (time_grid[end] != 1)
                #println("INFO: normalizing given time grid...")
                t0 = time_grid[1]
                tf = time_grid[end]
                time_grid = (time_grid .- t0) ./ (tf - t0) 
            end
            NLP_normalized_time_grid = time_grid
            dim_NLP_steps = length(time_grid) - 1
        end

        # additional indicators (faster ?)
        has_free_t0 = has_free_initial_time(ocp)
        has_free_tf = has_free_final_time(ocp)
        has_lagrange = has_lagrange_cost(ocp)
        has_mayer = has_mayer_cost(ocp)
        has_variable = is_variable_dependent(ocp)
        has_maximization = is_max(ocp)

        # dimensions
        if has_lagrange
            dim_NLP_x = ocp.state_dimension + 1
        else
            dim_NLP_x = ocp.state_dimension
        end
        dim_NLP_u = ocp.control_dimension
        if has_variable
            dim_NLP_v = ocp.variable_dimension
        else
            dim_NLP_v = 0 # dim in ocp would be Nothing
        end
        dim_OCP_x = ocp.state_dimension

        N = dim_NLP_steps

        # NLP unknown (state + control + variable)
        dim_NLP_variables = (N + 1) * (dim_NLP_x + dim_NLP_u) + dim_NLP_v

        # NLP constraints 
        # parse NLP constraints (and initialize dimensions)
        control_constraints, state_constraints, mixed_constraints, boundary_constraints, variable_constraints, control_box, state_box, variable_box = CTBase.nlp_constraints!(ocp)

        dim_x_box = dim_state_range(ocp)
        dim_u_box = dim_control_range(ocp)
        dim_v_box = dim_variable_range(ocp)
        dim_path_cons = dim_path_constraints(ocp)
        dim_x_cons = dim_state_constraints(ocp)
        dim_u_cons = dim_control_constraints(ocp)
        dim_v_cons = dim_variable_constraints(ocp)
        dim_mixed_cons = dim_mixed_constraints(ocp)
        dim_boundary_cons = dim_boundary_constraints(ocp)

        # lagrange to mayer transformation
        dim_NLP_constraints = N * (dim_NLP_x + dim_path_cons) + dim_path_cons + dim_boundary_cons + dim_v_cons
        if has_lagrange
            dim_NLP_constraints += 1          
        end

        # call constructor with const fields
        docp = new(ocp, control_constraints, state_constraints, mixed_constraints, boundary_constraints, variable_constraints, control_box, state_box, variable_box, has_free_t0, has_free_tf, has_lagrange, has_mayer, has_variable,
        has_maximization, dim_x_box, dim_u_box, dim_v_box, dim_path_cons, dim_x_cons, dim_u_cons, dim_v_cons,dim_mixed_cons, dim_boundary_cons, dim_NLP_x, dim_NLP_u, dim_NLP_v, dim_OCP_x, dim_NLP_steps, NLP_normalized_time_grid, dim_NLP_variables, dim_NLP_constraints,-Inf * ones(dim_NLP_variables), Inf * ones(dim_NLP_variables), zeros(dim_NLP_constraints), zeros(dim_NLP_constraints))

        return docp

    end

end


"""
$(TYPEDSIGNATURES)

Check if an OCP is solvable by the method [`solve`](@ref).
"""
function is_solvable(ocp)
    solvable = true
    return solvable
end


"""
$(TYPEDSIGNATURES)

Build upper and lower bounds vectors for the DOCP nonlinear constraints.
"""
function constraints_bounds!(docp::DOCP)

    lb = docp.con_l
    ub = docp.con_u

    index = 1 # counter for the constraints
    for i in 0:docp.dim_NLP_steps-1
        # skip (ie leave 0) for equality dynamics constraint
        index = index + docp.dim_NLP_x
        # path constraints
        index = setPathBounds!(docp, index, lb, ub)
    end
    
    # path constraints at final time
    index = setPathBounds!(docp, index, lb, ub) 

    # boundary and variable constraints
     index = setPointBounds!(docp, index, lb, ub)

    return lb, ub
end


"""
$(TYPEDSIGNATURES)

Build upper and lower bounds vectors for the DOCP variable box constraints.
"""
function variables_bounds!(docp::DOCP)

    N = docp.dim_NLP_steps
    var_l = docp.var_l
    var_u = docp.var_u
    ocp = docp.ocp

    # NB. keep offset for each block since they are optional !

    # build ordered bounds vectors for state and control
    x_lb = -Inf * ones(docp.dim_OCP_x)
    x_ub = Inf * ones(docp.dim_OCP_x)
    for j in 1:docp.dim_x_box
        indice = docp.state_box[2][j]
        x_lb[indice] = docp.state_box[1][j]
        x_ub[indice] = docp.state_box[3][j]
    end
    u_lb = -Inf * ones(docp.dim_NLP_u)
    u_ub = Inf * ones(docp.dim_NLP_u)
    for j in 1:docp.dim_u_box
        indice = docp.control_box[2][j]
        u_lb[indice] = docp.control_box[1][j]
        u_ub[indice] = docp.control_box[3][j]
    end

    # apply bounds for NLP variables
    for i in 0:N
        set_variables_at_time_step!(var_l, x_lb, u_lb, docp, i)
        set_variables_at_time_step!(var_u, x_ub, u_ub, docp, i)
    end

    # variable box
    offset = (N+1) * (docp.dim_NLP_x + docp.dim_NLP_u)
    if docp.dim_v_box > 0
        for j in 1:docp.dim_v_box
            indice = docp.variable_box[2][j]
            var_l[offset+indice] = docp.variable_box[1][j]
            var_u[offset+indice] = docp.variable_box[3][j]
        end
    end

    return var_l, var_u
end


"""
$(TYPEDSIGNATURES)

Compute the objective for the DOCP problem.
"""
# DOCP objective
function DOCP_objective(xu, docp::DOCP)

    obj = 0.
    N = docp.dim_NLP_steps
    ocp = docp.ocp

    # mayer cost
    if docp.has_mayer
        v = get_variable(xu, docp)
        t0 = get_initial_time(xu, docp)
        tf = get_final_time(xu, docp)
        x0,u0,xl0 = get_variables_at_time_step(xu, docp, 0)
        xf,uf,xlf = get_variables_at_time_step(xu, docp, N)
        #obj = obj + ocp.mayer(t0, tf, x0, xf, v)
        obj = obj + ocp.mayer(x0, xf, v)
    end
    
    # lagrange cost
    if docp.has_lagrange
        xf,uf,xlf = get_variables_at_time_step(xu, docp, N)
        obj = obj + xlf
    end

    # maximization problem
    if docp.has_maximization
        obj = -obj
    end 
    
    return obj
end


"""
$(TYPEDSIGNATURES)

Compute the constraints C for the DOCP problem (modeled as LB <= C(X) <= UB).
"""
function DOCP_constraints!(c, xu, docp::DOCP)    
    """
    compute the constraints for the NLP : 
        - discretization of the dynamics via the trapeze method
        - boundary conditions
    inputs
    ocp :: ocp model
    xu :: 
    return
    c :: 
    """

    # could use a single v, however v is set only in Args constructor, not update, so not much is wasted
    v = get_variable(xu, docp)

    # main loop on time steps
    args_i = ArgsAtTimeStep(xu, docp, 0, v)
    args_ip1 = ArgsAtTimeStep(xu, docp, 1, v)

    index = 1 # counter for the constraints
    for i in 0:docp.dim_NLP_steps-1

        # state equation
        index = setStateEquation!(docp, c, index, (args_i, args_ip1))
        # path constraints 
        index = setPathConstraints!(docp, c, index, args_i, v)

        # smart update for next iteration
        if i < docp.dim_NLP_steps-1
            args_i = args_ip1
            args_ip1 = ArgsAtTimeStep(xu, docp, i+2, v)
        end
    end

    # path constraints at final time
    args_0 = ArgsAtTimeStep(xu, docp, 0, v)
    args_f = ArgsAtTimeStep(xu, docp, docp.dim_NLP_steps, v)
    index = setPathConstraints!(docp, c, index, args_f, v)

    # boundary conditions and variable constraints
    index = setPointConstraints!(docp, c, index, args_0, args_f, v)

    # needed even for inplace version, AD error otherwise
    # may be because actual return would be index above ?
    return c 
end


# +++ later use abstract interface, based on Args variants ?
# this struct and the related functions will change according to discretization scheme
# put this part in trapeze.jl file
"""
$(TYPEDSIGNATURES)

Useful values at a time step: time, state, control, dynamics...
"""
struct ArgsAtTimeStep

    time
    state
    control
    dynamics
    lagrange_state
    lagrange_cost

    function ArgsAtTimeStep(xu, docp::DOCP, i::Int, v)

        # variables
        ti = get_time_at_time_step(xu, docp, i)
        #xi = get_state_at_time_step(xu, docp, i)
        #ui = get_control_at_time_step(xu, docp, i)
        xi, ui, xli = get_variables_at_time_step(xu, docp, i)

        # dynamics and lagrange cost
        fi = docp.ocp.dynamics(ti, xi, ui, v)

        if docp.has_lagrange
            #xli = get_lagrange_cost_at_time_step(xu, docp, i)
            li = docp.ocp.lagrange(ti, xi, ui, v)
            args = new(ti, xi, ui, fi, xli, li)
        else
            args = new(ti, xi, ui, fi)
        end

        return args
    end

end


"""
$(TYPEDSIGNATURES)

Set the constraints corresponding to the state equation
"""
function setStateEquation!(docp::DOCP, c, index::Int, args_trapeze)
    
    ocp = docp.ocp
    args_i = args_trapeze[1]
    args_ip1 = args_trapeze[2]
    hi = args_ip1.time - args_i.time
    
    # trapeze rule
    c[index:index+docp.dim_OCP_x-1] .= args_ip1.state .- (args_i.state .+ 0.5*hi*(args_i.dynamics .+ args_ip1.dynamics))

    if docp.has_lagrange
        c[index+docp.dim_OCP_x] = args_ip1.lagrange_state - (args_i.lagrange_state + 0.5*hi*(args_i.lagrange_cost + args_ip1.lagrange_cost))
    end
    
    index = index + docp.dim_NLP_x
    return index
end


"""
$(TYPEDSIGNATURES)

Set the path constraints at given time step
"""
function setPathConstraints!(docp::DOCP, c, index::Int, args::ArgsAtTimeStep, v)

    ocp = docp.ocp

    ti = args.time
    xi = args.state
    ui = args.control

    # pure control constraints
    if docp.dim_u_cons > 0
        c[index:index+docp.dim_u_cons-1] = docp.control_constraints[2](ti, ui, v)
        index = index + docp.dim_u_cons
    end

    # pure state constraints
    if docp.dim_x_cons > 0
        c[index:index+docp.dim_x_cons-1] = docp.state_constraints[2](ti, xi ,v)
        index = index + docp.dim_x_cons
    end

    # mixed state / control constraints
    if docp.dim_mixed_cons > 0
        c[index:index+docp.dim_mixed_cons-1] = docp.mixed_constraints[2](ti, xi, ui, v)
        index = index + docp.dim_mixed_cons
    end
    
    return index

end

"""
$(TYPEDSIGNATURES)

Set bounds for the path constraints at given time step
"""
function setPathBounds!(docp::DOCP, index::Int, lb, ub)

    ocp = docp.ocp

    # pure control constraints
    if docp.dim_u_cons > 0
        lb[index:index+docp.dim_u_cons-1] = docp.control_constraints[1]
        ub[index:index+docp.dim_u_cons-1] = docp.control_constraints[3]
        index = index + docp.dim_u_cons
    end

    # pure state constraints
    if docp.dim_x_cons > 0
        lb[index:index+docp.dim_x_cons-1] = docp.state_constraints[1]
        ub[index:index+docp.dim_x_cons-1] = docp.state_constraints[3]
        index = index + docp.dim_x_cons
    end

    # mixed state / control constraints
    if docp.dim_mixed_cons > 0
        lb[index:index+docp.dim_mixed_cons-1] = docp.mixed_constraints[1]
        ub[index:index+docp.dim_mixed_cons-1] = docp.mixed_constraints[3]
        index = index + docp.dim_mixed_cons
    end
    
    return index

end


"""
$(TYPEDSIGNATURES)

Set the boundary and variable constraints
"""
function setPointConstraints!(docp::DOCP, c, index::Int, args_0::ArgsAtTimeStep, args_f::ArgsAtTimeStep, v)

    ocp = docp.ocp

    x0 = args_0.state
    xf = args_f.state

    # boundary constraints
    if docp.dim_boundary_cons > 0
        c[index:index+docp.dim_boundary_cons-1] = docp.boundary_constraints[2](x0, xf, v)
        index = index + docp.dim_boundary_cons
    end

    # variable constraints
    if docp.dim_v_cons > 0
        c[index:index+docp.dim_v_cons-1] = docp.variable_constraints[2](v)
        index = index + docp.dim_v_cons
    end

    # null initial condition for lagrangian cost state
    if docp.has_lagrange
        c[index] = args_0.lagrange_state
        index = index + 1
    end

    return index

end


"""
$(TYPEDSIGNATURES)

Set bounds for the boundary and variable constraints
"""
function setPointBounds!(docp::DOCP, index::Int, lb, ub)

    ocp = docp.ocp

    # boundary constraints
    if docp.dim_boundary_cons > 0
        lb[index:index+docp.dim_boundary_cons-1] = docp.boundary_constraints[1]
        ub[index:index+docp.dim_boundary_cons-1] = docp.boundary_constraints[3]
        index = index + docp.dim_boundary_cons
    end

    # variable constraints
    if docp.dim_v_cons > 0
        lb[index:index+docp.dim_v_cons-1] = docp.variable_constraints[1]
        ub[index:index+docp.dim_v_cons-1] = docp.variable_constraints[3]
        index = index + docp.dim_v_cons
    end

    # null initial condition for lagrangian cost state
    if docp.has_lagrange
        lb[index] = 0.
        ub[index] = 0.
        index = index + 1
    end

    return index

end


"""
$(TYPEDSIGNATURES)

Check the nonlinear constraints violation for the DOCP problem. 
"""
function DOCP_constraints_check!(cb, constraints, docp)

    # +++ todo add a single utils function check_bounds(v,lb,ub) that returns the error vector

    # check constraints vs bounds
    # by construction only one of the two can be active
    for i in 1:docp.dim_NLP_constraints
        if constraints[i] < docp.con_l[i]
            cb[i] = constraints[i] - docp.con_l[i]
        end
        if constraints[i] > docp.con_u[i]
            cb[i] = constraints[i] - docp.con_u[i]
        end
    end
    return nothing
end


"""
$(TYPEDSIGNATURES)

Check the variables box constraints violation for the DOCP problem. 
"""
function DOCP_variables_check!(vb, variables, docp)
    # check variables vs bounds
    # by construction only one of the two can be active
    for i in 1:docp.dim_NLP_variables
        if variables[i] < docp.var_l[i]
            vb[i] = variables[i] - docp.var_l[i] # < 0
        end
        if variables[i] > docp.var_u[i]
            vb[i] = variables[i] - docp.var_u[i] # > 0
        end
    end
    return nothing
end
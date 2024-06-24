#=
NLP variables layout: [X0, X1 .. XN, U0, U1 .. UN, V]

additional state variable x_{n+1}(t) for the objective (Lagrange to Mayer formulation)
    [x_1(t_0), ... , x_{n+1}(t_0),
    ... , 
    x_{1}(t_N), ... , x_{n+1}(t_N),
    u_1(t_0), ... , u_m(t_0), 
    ... , 
    u_m(t_N), ..., u_m(t_N),
    v]
=#

"""
$(TYPEDSIGNATURES)

Struct for discretized optimal control problem DOCP

Contains:
- a copy of the original OCP
- a NLP formulation of the DOCP
- data required to link the two problems
"""
mutable struct DOCP

    #+++ constructor could use sub-functions ?
    #+++ add types, use const when possible

    ## OCP
    const ocp::OptimalControlModel

    # OCP variables and functions
    #variable_dimension::Int64

    # functions
    control_constraints
    state_constraints
    mixed_constraints
    boundary_constraints
    variable_constraints
    control_box
    state_box
    variable_box

    ## NLP
    dim_NLP_x::Int64  # possible additional lagrange cost
    dim_NLP_u::Int64
    dim_NLP_v::Int64
    dim_NLP_constraints::Int64
    dim_NLP_variables::Int64
    dim_NLP_steps::Int64
    NLP_normalized_time_grid

    # lower and upper bounds for variables and constraints
    var_l
    var_u
    con_l
    con_u

    # NLP model for solver
    nlp

    # constructor
    function DOCP(ocp::OptimalControlModel, grid_size::Integer, time_grid)       

        # +++ try to put here more const members (indicators etc), also move some parts to CTBase
        docp = new(ocp)

        ## Optimal Control Problem OCP
        # time grid
        if time_grid == nothing
            docp.NLP_normalized_time_grid = collect(LinRange(0, 1, grid_size+1))
            docp.dim_NLP_steps = grid_size
        else
            # check strictly increasing
            if !issorted(time_grid,lt=<=)
                throw(ArgumentError("given time grid is not strictly increasing. Aborting..."))
                return docp
            end
            # normalize input grid if needed
            if (time_grid[1] != 0) || (time_grid[end] != 1)
                println("INFO: normalizing given time grid...")
                t0 = time_grid[1]
                tf = time_grid[end]
                time_grid = (time_grid .- t0) ./ (tf - t0) 
            end
            docp.NLP_normalized_time_grid = time_grid
            docp.dim_NLP_steps = length(time_grid) - 1
        end
        N = docp.dim_NLP_steps

        # parse NLP constraints
        docp.control_constraints, docp.state_constraints, docp.mixed_constraints, docp.boundary_constraints, docp.variable_constraints, docp.control_box, docp.state_box, docp.variable_box = nlp_constraints!(ocp)

        # set dimensions
        # Mayer to Lagrange: additional state with Lagrange cost as dynamics and null initial condition
        if has_lagrange_cost(ocp)
            docp.dim_NLP_x = docp.ocp.state_dimension + 1  
            docp.dim_NLP_constraints = N * (docp.dim_NLP_x + dim_path_constraints(ocp)) + dim_path_constraints(ocp) + dim_boundary_constraints(ocp) + dim_variable_constraints(ocp) + 1           
        else
            docp.dim_NLP_x = docp.ocp.state_dimension  
            docp.dim_NLP_constraints = N * (docp.dim_NLP_x + dim_path_constraints(ocp)) + dim_path_constraints(ocp) + dim_boundary_constraints(ocp) + dim_variable_constraints(ocp)
        end

        docp.dim_NLP_u = ocp.control_dimension
        
        if is_variable_dependent(ocp)
            docp.dim_NLP_v = ocp.variable_dimension
        else
            docp.dim_NLP_v = 0 # dim in ocp would be 'nothing'
        end

        # NLP unknown (state + control + variable)
        docp.dim_NLP_variables = (N + 1) * (docp.dim_NLP_x + docp.dim_NLP_u) + docp.dim_NLP_v

        return docp

    end

end


"""
$(TYPEDSIGNATURES)

Check if an OCP is solvable by the method [`solve`](@ref).
"""
function is_solvable(ocp)
    solvable = true
    # +++ note: non-autonomous mayer case is not supported
    return solvable
end


"""
$(TYPEDSIGNATURES)

Build upper and lower bounds vectors for the DOCP nonlinear constraints.
"""
function constraints_bounds(docp)

    lb = zeros(docp.dim_NLP_constraints)
    ub = zeros(docp.dim_NLP_constraints)

    index = 1 # counter for the constraints
    for i in 0:docp.dim_NLP_steps-1
        # skip (ie leave 0) for equality dynamics constraint
        index = index + docp.dim_NLP_x
        # path constraints
        index = setPathConstraintsAtTimeStep!(docp, index, :bounds; lb=lb, ub=ub)
    end
    
    # path constraints at final time
    index = setPathConstraintsAtTimeStep!(docp, index, :bounds; lb=lb, ub=ub) 

    # boundary and variable constraints
     index = setPunctualConditions!(docp, index, :bounds; lb=lb, ub=ub)

    return lb, ub
end


"""
$(TYPEDSIGNATURES)

Build upper and lower bounds vectors for the DOCP variable box constraints.
"""
function variables_bounds(docp)

    N = docp.dim_NLP_steps
    l_var = -Inf * ones(docp.dim_NLP_variables)
    u_var = Inf * ones(docp.dim_NLP_variables)
    ocp = docp.ocp

    # NB. keep offset for each block since they are optional !
    # Also, not practical to reuse the setters for x,u,v due to the non-ordered indices and possibly not full dimension

    # state box
    offset = 0
    if dim_state_box(ocp) > 0
        for i in 0:N
            for j in 1:dim_state_box(ocp)
                indice = docp.state_box[2][j]
                l_var[offset+indice] = docp.state_box[1][j]
                u_var[offset+indice] = docp.state_box[3][j]
            end
            offset = offset + docp.dim_NLP_x
        end
    end

    # control box
    offset = (N+1) * docp.dim_NLP_x
    if dim_control_box(ocp) > 0
        for i in 0:N
            for j in 1:dim_control_box(ocp)
                indice = docp.control_box[2][j]
                l_var[offset+indice] = docp.control_box[1][j]
                u_var[offset+indice] = docp.control_box[3][j]
            end
            offset = offset + docp.dim_NLP_u
        end
    end

    # variable box
    offset = (N+1) * (docp.dim_NLP_x + docp.dim_NLP_u)
    if dim_variable_box(ocp) > 0
        for j in 1:dim_variable_box(ocp)
            indice = docp.variable_box[2][j]
            l_var[offset+indice] = docp.variable_box[1][j]
            u_var[offset+indice] = docp.variable_box[3][j]
        end
    end

    return l_var, u_var
end


"""
$(TYPEDSIGNATURES)

Compute the objective for the DOCP problem.
"""
# DOCP objective
function DOCP_objective(xu, docp)

    obj = 0
    N = docp.dim_NLP_steps
    ocp = docp.ocp

    # note: non-autonomous mayer case is not supported
    if has_mayer_cost(ocp)
        v = get_variable(xu, docp)
        x0 = vec2scal(get_state_at_time_step(xu, docp, 0))
        xf = vec2scal(get_state_at_time_step(xu, docp, N))
        obj = obj + ocp.mayer(x0, xf, v)
    end
    
    if has_lagrange_cost(ocp)
        obj = obj + xu[(N+1)*docp.dim_NLP_x]
    end

    if is_min(ocp)
        return obj
    else
        return -obj
    end
end


"""
$(TYPEDSIGNATURES)

Compute the constraints C for the DOCP problem (modeled as LB <= C(X) <= UB).
"""
function DOCP_constraints!(c, xu, docp)    
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

    # t,x,u,f,... at t_0
    args_0 = ArgsAtTimeStep(xu, docp, 0)
    args_i = args_0

    # main loop on time steps
    index = 1 # counter for the constraints
    for i in 0:docp.dim_NLP_steps-1

        # t,x,u,f,... at t_{i+1}
        args_ip1 = ArgsAtTimeStep(xu, docp, i+1)

        # state equation
        index = setStateEquationAtTimeStep!(docp, c, index, args_i, args_ip1)

        # path constraints 
        index = setPathConstraintsAtTimeStep!(docp, index, :constraints; c=c, args=args_i)

        # updates for next iteration
        args_i = args_ip1
    end

    # path constraints at final time
    args_f = args_i
    index = setPathConstraintsAtTimeStep!(docp, index, :constraints; c=c, args=args_f)

    # boundary conditions and variable constraints
    index = setPunctualConditions!(docp, index, :constraints; c=c, args_0=args_0, args_f=args_f)

    # needed even for inplace version, AD error otherwise oO
    return c 
end


# +++ later use abstract interface
"""
$(TYPEDSIGNATURES)

Useful values at a time step: time, state, control, dynamics...
"""
mutable struct ArgsAtTimeStep
    time
    state
    control
    variable
    dynamics
    lagrange_state
    lagrange_cost

    function ArgsAtTimeStep(xu, docp, i)
        args = new()
        args.time = get_time_at_time_step(xu, docp, i)
        args.state = vec2scal(get_state_at_time_step(xu, docp, i))
        args.control = vec2scal(get_control_at_time_step(xu, docp, i))
        args.variable = get_variable(xu, docp)
        args.dynamics = docp.ocp.dynamics(args.time, args.state, args.control, args.variable)
        if has_lagrange_cost(docp.ocp)
            args.lagrange_state = get_lagrange_cost_at_time_step(xu, docp, i)
            args.lagrange_cost = docp.ocp.lagrange(args.time, args.state, args.control, args.variable)
        else
            args.lagrange_state = 0
            args.lagrange_cost = 0
        end
        return args
    end
end

"""
$(TYPEDSIGNATURES)

Set the constraints corresponding to the state equation
"""
function setStateEquationAtTimeStep!(docp, c, index, args_i, args_ip1)
    
    ocp = docp.ocp
    hi = args_ip1.time - args_i.time
    
    if ocp.state_dimension == 1
        c[index] = args_ip1.state - (args_i.state + 0.5*hi*(args_i.dynamics + args_ip1.dynamics))            
    else
        toto = args_i.state .+ 0.5*hi*(args_i.dynamics + args_ip1.dynamics)
        c[index:index+ocp.state_dimension-1] = args_ip1.state - toto
    end

    if has_lagrange_cost(ocp)
        c[index+ocp.state_dimension] = args_ip1.lagrange_state - (args_i.lagrange_state + 0.5*hi*(args_i.lagrange_cost + args_ip1.lagrange_cost))
    end
    
    index = index + docp.dim_NLP_x
    return index
end


"""
$(TYPEDSIGNATURES)

Set the path constraints / bounds for given time step
target = :constraints | :bounds
"""
function setPathConstraintsAtTimeStep!(docp, index, target; c=nothing, args=nothing, lb=nothing, ub=nothing)

    ocp = docp.ocp
    if target == :constraints
        ti = args.time
        xi = args.state
        ui = args.control
        v = args.variable
    end

    # pure control constraints
    if dim_control_constraints(ocp) > 0
        if target == :constraints
            c[index:index+dim_control_constraints(ocp)-1] = docp.control_constraints[2](ti, ui, v)
        else
            lb[index:index+dim_control_constraints(ocp)-1] = docp.control_constraints[1]
            ub[index:index+dim_control_constraints(ocp)-1] = docp.control_constraints[3]
        end
        index = index + dim_control_constraints(ocp)
    end

    # pure state constraints
    if dim_state_constraints(ocp) > 0
        if target == :constraints
            c[index:index+dim_state_constraints(ocp)-1] = docp.state_constraints[2](ti, xi ,v)
        else
            lb[index:index+dim_state_constraints(ocp)-1] = docp.state_constraints[1]
            ub[index:index+dim_state_constraints(ocp)-1] = docp.state_constraints[3]
        end
        index = index + dim_state_constraints(ocp)
    end

    # mixed state / control constraints
    if dim_mixed_constraints(ocp) > 0
        if target == :constraints
            c[index:index+dim_mixed_constraints(ocp)-1] = docp.mixed_constraints[2](ti, xi, ui, v)
        else
            lb[index:index+dim_mixed_constraints(ocp)-1] = docp.mixed_constraints[1]
            ub[index:index+dim_mixed_constraints(ocp)-1] = docp.mixed_constraints[3]
        end
        index = index + dim_mixed_constraints(ocp)
    end
    
    return index

end


"""
$(TYPEDSIGNATURES)

Set the boundary and variable constraints / their bounds
target = :constraints | :bounds
"""
function setPunctualConditions!(docp, index, target; c=nothing, args_0=nothing, args_f=nothing, lb=nothing, ub=nothing)

    ocp = docp.ocp

    # boundary constraints
    if dim_boundary_constraints(ocp) > 0
        if target == :constraints
            c[index:index+dim_boundary_constraints(ocp)-1] = docp.boundary_constraints[2](args_0.state, args_f.state, args_0.variable)
        else
            lb[index:index+dim_boundary_constraints(ocp)-1] = docp.boundary_constraints[1]
            ub[index:index+dim_boundary_constraints(ocp)-1] = docp.boundary_constraints[3]
        end
        index = index + dim_boundary_constraints(ocp)
    end

    # variable constraints
    if dim_variable_constraints(ocp) > 0
        if target == :constraints
            c[index:index+dim_variable_constraints(ocp)-1] = docp.variable_constraints[2](args_0.variable)
        else
            lb[index:index+dim_variable_constraints(ocp)-1] = docp.variable_constraints[1]
            ub[index:index+dim_variable_constraints(ocp)-1] = docp.variable_constraints[3]
        end
        index = index + dim_variable_constraints(ocp)
    end

    # null initial condition for lagrangian cost state
    if has_lagrange_cost(ocp)
        if target == :constraints
            c[index] = args_0.lagrange_state
        else
            lb[index] = 0.
            ub[index] = 0.
        end
        index = index + 1
    end

    return index

end


"""
$(TYPEDSIGNATURES)

Check the nonlinear constraints violation for the DOCP problem. 
"""
function DOCP_constraints_check!(cb, constraints, docp)

    # +++ todo unify in a single utils function check_bounds(v,lb,ub) that returns the error vector

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
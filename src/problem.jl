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
    variable_dimension::Int64
    #has_free_initial_time::Bool
    #has_free_final_time::Bool
    #has_variable::Bool
    #has_lagrange_cost::Bool
    #has_mayer_cost::Bool

    # OCP constraints
    # indicators
    has_control_constraints::Bool
    has_state_constraints::Bool
    has_mixed_constraints::Bool
    has_boundary_conditions::Bool
    has_variable_constraints::Bool
    has_control_box::Bool
    has_state_box::Bool
    has_variable_box::Bool

    # dimensions
    dim_control_constraints::Int64
    dim_state_constraints::Int64
    dim_mixed_constraints::Int64
    dim_path_constraints::Int64
    dim_boundary_conditions::Int64
    dim_variable_constraints::Int64
    dim_control_box::Int64
    dim_state_box::Int64
    dim_variable_box::Int64
    
    # functions
    control_constraints
    state_constraints
    mixed_constraints
    boundary_conditions
    variable_constraints
    control_box
    state_box
    variable_box

    ## NLP
    # NLP problem
    dim_NLP_state::Int64  # 'augmented state' with possible additional component for lagrange cost
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

        # +++ try to put here more const members (indicators etc)
        # +++ also move some parts to CTBase.OptimalControlProblem
        docp = new(ocp)

        ## Optimal Control Problem OCP
        # time grid
        #docp.has_free_initial_time = (typeof(ocp.initial_time)==Index)
        #docp.has_free_final_time = (typeof(ocp.final_time)==Index)
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

        # dimensions and functions
        #docp.has_variable = !isnothing(ocp.variable_dimension)
        #if docp.has_variable
        if is_variable_dependent(ocp)
            docp.variable_dimension = ocp.variable_dimension
        else
            docp.variable_dimension = 0
        end
        #docp.has_lagrange_cost = !isnothing(ocp.lagrange)
        #docp.has_mayer_cost = !isnothing(ocp.mayer)
        
        # constraints
        docp.control_constraints, docp.state_constraints, docp.mixed_constraints, docp.boundary_conditions, docp.variable_constraints, docp.control_box, docp.state_box, docp.variable_box = nlp_constraints(ocp)

        docp.dim_control_constraints = length(docp.control_constraints[1])
        docp.dim_state_constraints = length(docp.state_constraints[1])
        docp.dim_mixed_constraints = length(docp.mixed_constraints[1])
        docp.dim_path_constraints = docp.dim_control_constraints + docp.dim_state_constraints + docp.dim_mixed_constraints
        docp.dim_boundary_conditions = length(docp.boundary_conditions[1])
        docp.dim_variable_constraints = length(docp.variable_constraints[1])
        docp.dim_control_box = length(docp.control_box[1])
        docp.dim_state_box = length(docp.state_box[1])
        docp.dim_variable_box = length(docp.variable_box[1])
        docp.has_control_constraints = !isempty(docp.control_constraints[1])
        docp.has_state_constraints = !isempty(docp.state_constraints[1])
        docp.has_mixed_constraints = !isempty(docp.mixed_constraints[1])
        docp.has_boundary_conditions = !isempty(docp.boundary_conditions[1])
        docp.has_variable_constraints = !isempty(docp.variable_constraints[1])
        docp.has_control_box = !isempty(docp.control_box[1])
        docp.has_state_box = !isempty(docp.state_box[1])
        docp.has_variable_box = !isempty(docp.variable_box[1])
  
        ## Non Linear Programming NLP

        # Mayer to Lagrange reformulation: 
        # additional state with Lagrange cost as dynamics and null initial condition
        if has_lagrange_cost(ocp)
            docp.dim_NLP_state = docp.ocp.state_dimension + 1  
            docp.dim_NLP_constraints = N * (docp.dim_NLP_state + docp.dim_path_constraints) + docp.dim_path_constraints + docp.dim_boundary_conditions + docp.dim_variable_constraints + 1           
        else
            docp.dim_NLP_state = docp.ocp.state_dimension  
            docp.dim_NLP_constraints = N * (docp.dim_NLP_state + docp.dim_path_constraints) + docp.dim_path_constraints + docp.dim_boundary_conditions + docp.dim_variable_constraints
        end

        # set dimension for NLP unknown (state + control + variable)
        docp.dim_NLP_variables = (N + 1) * (docp.dim_NLP_state + docp.ocp.control_dimension) + docp.variable_dimension

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


# +++ for the aux functions manipulating X and C(X) and bounds
# move to a more abstract level to make scheme change easier
# eg use addDynamicsConstraint, addBoundsBlock etc
# with abstract interfaces, to be implemented for each scheme

"""
$(TYPEDSIGNATURES)

Build upper and lower bounds vectors for the DOCP nonlinear constraints.
"""
function constraints_bounds(docp)

    N = docp.dim_NLP_steps
    lb = zeros(docp.dim_NLP_constraints)
    ub = zeros(docp.dim_NLP_constraints)
    ocp = docp.ocp

    index = 1 # counter for the constraints
    for i in 0:N-1
        # skip (ie leave 0) bound for equality dynamics constraint
        index = index + docp.dim_NLP_state
        # path constraints 
        if docp.has_control_constraints
            lb[index:index+docp.dim_control_constraints-1] = docp.control_constraints[1]
            ub[index:index+docp.dim_control_constraints-1] = docp.control_constraints[3]
            index = index + docp.dim_control_constraints
        end
        if docp.has_state_constraints
            lb[index:index+docp.dim_state_constraints-1] = docp.state_constraints[1]
            ub[index:index+docp.dim_state_constraints-1] = docp.state_constraints[3]
            index = index + docp.dim_state_constraints
        end
        if docp.has_mixed_constraints
            lb[index:index+docp.dim_mixed_constraints-1] = docp.mixed_constraints[1]
            ub[index:index+docp.dim_mixed_constraints-1] = docp.mixed_constraints[3]
            index = index + docp.dim_mixed_constraints
        end
    end
    
    # path constraints at final time
    if docp.has_control_constraints
        lb[index:index+docp.dim_control_constraints-1] = docp.control_constraints[1]
        ub[index:index+docp.dim_control_constraints-1] = docp.control_constraints[3]
        index = index + docp.dim_control_constraints
    end
    if docp.has_state_constraints
        lb[index:index+docp.dim_state_constraints-1] = docp.state_constraints[1]
        ub[index:index+docp.dim_state_constraints-1] = docp.state_constraints[3]
        index = index + docp.dim_state_constraints
    end
    if docp.has_mixed_constraints
        lb[index:index+docp.dim_mixed_constraints-1] = docp.mixed_constraints[1]
        ub[index:index+docp.dim_mixed_constraints-1] = docp.mixed_constraints[3]
        index = index + docp.dim_mixed_constraints
    end
    
    # boundary conditions
    if docp.has_boundary_conditions
        lb[index:index+docp.dim_boundary_conditions-1] = docp.boundary_conditions[1]
        ub[index:index+docp.dim_boundary_conditions-1] = docp.boundary_conditions[3]
        index = index + docp.dim_boundary_conditions
    end

    # variable constraints
    if docp.has_variable_constraints
        lb[index:index+docp.dim_variable_constraints-1] = docp.variable_constraints[1]
        ub[index:index+docp.dim_variable_constraints-1] = docp.variable_constraints[3]
        index = index + docp.dim_variable_constraints
    end 

    # lagrange cost (set integral to 0 at t0)
    if has_lagrange_cost(ocp)
        lb[index] = 0.
        ub[index] = 0.
        index = index + 1
    end

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

    # NLP variables layout: [X0, X1 .. XN, U0, U1 .. UN, V]
    # NB. keep offset for each block since blocks are optional !

    # +++ we could use the setters here (build local vectors then call setter ?!)

    # state box
    offset = 0
    if docp.has_state_box
        for i in 0:N
            for j in 1:docp.dim_state_box
                indice = docp.state_box[2][j]
                l_var[offset+indice] = docp.state_box[1][j]
                u_var[offset+indice] = docp.state_box[3][j]
            end
            offset = offset + docp.dim_NLP_state
        end
    end

    # control box
    offset = (N+1) * docp.dim_NLP_state
    if docp.has_control_box
        for i in 0:N
            for j in 1:docp.dim_control_box
                indice = docp.control_box[2][j]
                l_var[offset+indice] = docp.control_box[1][j]
                u_var[offset+indice] = docp.control_box[3][j]
            end
            offset = offset + docp.ocp.control_dimension
        end
    end

    # variable box
    offset = (N+1) * (docp.dim_NLP_state + docp.ocp.control_dimension)
    if docp.has_variable_box
        for j in 1:docp.dim_variable_box
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
        x0 = get_state_at_time_step(xu, docp, 0)
        xf = get_state_at_time_step(xu, docp, N)
        obj = obj + ocp.mayer(x0, xf, v)
    end
    
    if has_lagrange_cost(ocp)
        obj = obj + xu[(N+1)*docp.dim_NLP_state]
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
        layout of the nlp unknown xu for trapeze discretization 
        additional state variable x_{n+1}(t) for the objective (Lagrange to Mayer formulation)
        [x_1(t_0), ... , x_{n+1}(t_0),
        ... , 
        x_{1}(t_N), ... , x_{n+1}(t_N),
        u_1(t_0), ... , u_m(t_0), 
        ... , 
        u_m(t_N), ..., u_m(t_N),
        v]
    return
    c :: 
    """

    # initialize main loop on time steps
    N = docp.dim_NLP_steps
    v = get_variable(xu, docp)
    ocp = docp.ocp

    # time, state and control at t_0
    ti = get_time_at_time_step(xu, docp, 0)
    xi = get_state_at_time_step(xu, docp, 0)
    ui = get_control_at_time_step(xu, docp, 0)
    fi = ocp.dynamics(ti, xi, ui, v)
    if has_lagrange_cost(ocp)
        xli = get_lagrange_cost_at_time_step(xu, docp, 0)
        li = ocp.lagrange(ti, xi, ui, v)
    end

    # main loop on time steps
    index = 1 # counter for the constraints
    for i in 0:N-1

        # time, state and control at t_{i+1}
        tip1 = get_time_at_time_step(xu, docp, i+1)
        xip1 = get_state_at_time_step(xu, docp, i+1)
        uip1 = get_control_at_time_step(xu, docp, i+1)
        fip1 = ocp.dynamics(tip1, xip1, uip1, v)
        hi = tip1 - ti

        # state equation
        if ocp.state_dimension == 1
            c[index] = xip1 - (xi + 0.5*hi*(fi + fip1))            
        else
            c[index:index+ocp.state_dimension-1] = xip1 - (xi + 0.5*hi*(fi + fip1))
        end
        if has_lagrange_cost(ocp)
            xlip1 = get_lagrange_cost_at_time_step(xu, docp, i+1)
            lip1 = ocp.lagrange(tip1, xip1, uip1, v)
            c[index+ocp.state_dimension] = xlip1 - (xli + 0.5*hi*(li + lip1))
            xli = xlip1
            li = lip1
        end
        index = index + docp.dim_NLP_state

        # path constraints 
        # +++use aux function for block, see solution also
        if docp.has_control_constraints
            c[index:index+docp.dim_control_constraints-1] = docp.control_constraints[2](ti, ui, v)
            index = index + docp.dim_control_constraints
        end
        if docp.has_state_constraints
            c[index:index+docp.dim_state_constraints-1] = docp.state_constraints[2](ti, xi ,v)
            index = index + docp.dim_state_constraints
        end
        if docp.has_mixed_constraints
            c[index:index+docp.dim_mixed_constraints-1] = docp.mixed_constraints[2](ti, xi, ui, v)
            index = index + docp.dim_mixed_constraints
        end

        # updates for next iteration
        ti = tip1
        xi = xip1
        ui = uip1
        fi = fip1
    end

    # path constraints at final time
    tf = get_time_at_time_step(xu, docp, N)
    xf = get_state_at_time_step(xu, docp, N)
    uf = get_control_at_time_step(xu, docp, N)
    if docp.has_control_constraints
        c[index:index+docp.dim_control_constraints-1] = docp.control_constraints[2](tf, uf, v)      
        index = index + docp.dim_control_constraints
    end  
    if docp.has_state_constraints
        c[index:index+docp.dim_state_constraints-1] = docp.state_constraints[2](tf, xf, v)      
        index = index + docp.dim_state_constraints
    end 
    if docp.has_mixed_constraints
        c[index:index+docp.dim_mixed_constraints-1] = docp.mixed_constraints[2](tf, xf, uf, v)
        index = index + docp.dim_mixed_constraints
    end

    # boundary conditions
    if docp.has_boundary_conditions
        x0 = get_state_at_time_step(xu, docp, 0)
        c[index:index+docp.dim_boundary_conditions-1] = docp.boundary_conditions[2](x0, xf, v)
        index = index + docp.dim_boundary_conditions
    end

    # variable constraints
    if docp.has_variable_constraints
        c[index:index+docp.dim_variable_constraints-1] = docp.variable_constraints[2](v)
        index = index + docp.dim_variable_constraints
    end

    # null initial condition for augmented state (reformulated lagrangian cost)
    if has_lagrange_cost(ocp)
        c[index] = get_lagrange_cost_at_time_step(xu, docp, 0)
        index = index + 1
    end
    return c # needed even for inplace version, AD error otherwise oO
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
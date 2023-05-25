# struct for ocop/nlp info
mutable struct CTDirect_data

    ## OCP
    # OCP variables and functions
    initial_time
    final_time
    state_dimension
    control_dimension
    variable_dimension
    dynamics
    mayer
    lagrange
    criterion_min_max
    has_free_initial_time
    has_free_final_time
    has_variable
    has_lagrange_cost
    has_mayer_cost

    # OCP constraints
    # indicators
    has_control_constraints
    has_state_constraints
    has_mixed_constraints
    has_boundary_conditions
    has_variable_constraints
    has_control_box
    has_state_box
    has_variable_box
    has_scalar_state   
    has_scalar_control

    # dimensions
    dim_control_constraints
    dim_state_constraints
    dim_mixed_constraints
    dim_path_constraints
    dim_boundary_conditions
    dim_variable_constraints
    dim_control_box
    dim_state_box
    dim_variable_box

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
    dim_NLP_state  # 'augmented state'
    dim_NLP_constraints
    dim_NLP_variables
    dim_NLP_steps
    NLP_init

    # NLP solution
    NLP_solution
    NLP_objective
    NLP_sol_constraints
    NLP_constraints_violation
    NLP_iterations
    NLP_stats       # remove later ? type is https://juliasmoothoptimizers.github.io/SolverCore.jl/stable/reference/#SolverCore.GenericExecutionStats

    function CTDirect_data(ocp::OptimalControlModel, N::Integer, init=nothing)

        ctd = new()

        ## Optimal Control Problem OCP
        # time
        ctd.initial_time = ocp.initial_time
        ctd.final_time = ocp.final_time
        ctd.has_free_initial_time = (typeof(ocp.initial_time)==Index)
        ctd.has_free_final_time = (typeof(ocp.final_time)==Index)
        
        # dimensions and functions
        ctd.state_dimension = ocp.state_dimension
        ctd.has_scalar_state = (ctd.state_dimension == 1)    
        ctd.control_dimension = ocp.control_dimension
        ctd.has_scalar_control = (ctd.control_dimension == 1)   
        ctd.has_variable = !isnothing(ocp.variable_dimension)
        if ctd.has_variable
            ctd.variable_dimension = ocp.variable_dimension
        else
            ctd.variable_dimension = 0
        end
        println("Free t0: ", ctd.has_free_initial_time, " Free tf: ", ctd.has_free_final_time, " Variable dim: ", ctd.variable_dimension, " has scalar state: ", ctd.has_scalar_state, " has scalar control: ", ctd.has_scalar_control)
        ctd.dynamics = ocp.dynamics
        ctd.has_lagrange_cost = !isnothing(ocp.lagrange)
        ctd.lagrange = ocp.lagrange
        ctd.has_mayer_cost = !isnothing(ocp.mayer)
        ctd.mayer = ocp.mayer
        
        # constraints
        ctd.control_constraints, ctd.state_constraints, ctd.mixed_constraints, ctd.boundary_conditions, ctd.variable_constraints, ctd.control_box, ctd.state_box, ctd.variable_box = nlp_constraints(ocp)
        ctd.dim_control_constraints = length(ctd.control_constraints[1])
        ctd.dim_state_constraints = length(ctd.state_constraints[1])
        ctd.dim_mixed_constraints = length(ctd.mixed_constraints[1])
        ctd.dim_path_constraints = ctd.dim_control_constraints + ctd.dim_state_constraints + ctd.dim_mixed_constraints
        ctd.dim_boundary_conditions = length(ctd.boundary_conditions[1])
        ctd.dim_variable_constraints = length(ctd.variable_constraints[1])
        ctd.dim_control_box = length(ctd.control_box[1])
        ctd.dim_state_box = length(ctd.state_box[1])
        ctd.dim_variable_box = length(ctd.variable_box[1])
        ctd.has_control_constraints = !isempty(ctd.control_constraints[1])
        ctd.has_state_constraints = !isempty(ctd.state_constraints[1])
        ctd.has_mixed_constraints = !isempty(ctd.mixed_constraints[1])
        ctd.has_boundary_conditions = !isempty(ctd.boundary_conditions[1])
        ctd.has_variable_constraints = !isempty(ctd.variable_constraints[1])
        ctd.has_control_box = !isempty(ctd.control_box[1])
        ctd.has_state_box = !isempty(ctd.state_box[1])
        ctd.has_variable_box = !isempty(ctd.variable_box[1])

        ## Non Linear Programming NLP
        ctd.dim_NLP_steps = N
        ctd.NLP_init = init

        # Mayer to Lagrange reformulation: 
        # additional state with Lagrange cost as dynamics and null initial condition
        if ctd.has_lagrange_cost
            ctd.dim_NLP_state = ctd.state_dimension + 1  
            ctd.dim_NLP_constraints = N * (ctd.dim_NLP_state + ctd.dim_path_constraints) + ctd.dim_path_constraints + ctd.dim_boundary_conditions + ctd.dim_variable_constraints + 1           
        else
            ctd.dim_NLP_state = ctd.state_dimension  
            ctd.dim_NLP_constraints = N * (ctd.dim_NLP_state + ctd.dim_path_constraints) + ctd.dim_path_constraints + ctd.dim_boundary_conditions + ctd.dim_variable_constraints
        end

        # min or max problem (unused ?)
        ctd.criterion_min_max = ocp.criterion

        # set dimension for NLP unknown (state + control + variable)
        ctd.dim_NLP_variables = (N + 1) * (ctd.dim_NLP_state + ctd.control_dimension) + ctd.variable_dimension

        return ctd

    end

end


function is_solvable(ocp)
    solvable = true
    return solvable
end


# bounds for the constraints
function constraints_bounds(ctd)

    N = ctd.dim_NLP_steps
    lb = zeros(ctd.dim_NLP_constraints)
    ub = zeros(ctd.dim_NLP_constraints)

    index = 1 # counter for the constraints
    for i in 0:N-1
        # skip (ie leave 0) bound for equality dynamics constraint
        index = index + ctd.dim_NLP_state
        # path constraints 
        if ctd.has_control_constraints
            lb[index:index+ctd.dim_control_constraints-1] = ctd.control_constraints[1]
            ub[index:index+ctd.dim_control_constraints-1] = ctd.control_constraints[3]
            index = index + ctd.dim_control_constraints
        end
        if ctd.has_state_constraints
            lb[index:index+ctd.dim_state_constraints-1] = ctd.state_constraints[1]
            ub[index:index+ctd.dim_state_constraints-1] = ctd.state_constraints[3]
            index = index + ctd.dim_state_constraints
        end
        if ctd.has_mixed_constraints
            lb[index:index+ctd.dim_mixed_constraints-1] = ctd.mixed_constraints[1]
            ub[index:index+ctd.dim_mixed_constraints-1] = ctd.mixed_constraints[3]
            index = index + ctd.dim_mixed_constraints
        end
    end
    
    # path constraints at final time
    if ctd.has_control_constraints
        lb[index:index+ctd.dim_control_constraints-1] = ctd.control_constraints[1]
        ub[index:index+ctd.dim_control_constraints-1] = ctd.control_constraints[3]
        index = index + ctd.dim_control_constraints
    end
    if ctd.has_state_constraints
        lb[index:index+ctd.dim_state_constraints-1] = ctd.state_constraints[1]
        ub[index:index+ctd.dim_state_constraints-1] = ctd.state_constraints[3]
        index = index + ctd.dim_state_constraints
    end
    if ctd.has_mixed_constraints
        lb[index:index+ctd.dim_mixed_constraints-1] = ctd.mixed_constraints[1]
        ub[index:index+ctd.dim_mixed_constraints-1] = ctd.mixed_constraints[3]
        index = index + ctd.dim_mixed_constraints
    end
    
    # boundary conditions
    if ctd.has_boundary_conditions
        lb[index:index+ctd.dim_boundary_conditions-1] = ctd.boundary_conditions[1]
        ub[index:index+ctd.dim_boundary_conditions-1] = ctd.boundary_conditions[3]
        index = index + ctd.dim_boundary_conditions
    end

    # variable constraints
    if ctd.has_variable_constraints
        lb[index:index+ctd.dim_variable_constraints-1] = ctd.variable_constraints[1]
        ub[index:index+ctd.dim_variable_constraints-1] = ctd.variable_constraints[3]
        index = index + ctd.dim_variable_constraints
    end 

    # lagrange cost (set integral to 0 at t0)
    if ctd.has_lagrange_cost
        lb[index] = 0.
        ub[index] = 0.
        index = index + 1
    end

    return lb, ub
end

# box constraints for variables
function variables_bounds(ctd)

    N = ctd.dim_NLP_steps
    l_var = -Inf*ones(ctd.dim_NLP_variables)
    u_var = Inf*ones(ctd.dim_NLP_variables)

    # NLP variables layout: [X0, X1 .. XN, U0, U1 .. UN, V]
    # NB. keep offset for each block since blocks are optional !

    # state box
    offset = 0
    if ctd.has_state_box
        for i in 0:N
            for j in 1:ctd.dim_state_box
                indice = ctd.state_box[2][j]
                l_var[offset+indice] = ctd.state_box[1][j]
                u_var[offset+indice] = ctd.state_box[3][j]
            end
            offset = offset + ctd.dim_NLP_state
        end
    end

    # control box
    offset = (N+1) * ctd.dim_NLP_state
    if ctd.has_control_box
        for i in 0:N
            for j in 1:ctd.dim_control_box
                indice = ctd.control_box[2][j]
                l_var[offset+indice] = ctd.control_box[1][j]
                u_var[offset+indice] = ctd.control_box[3][j]
            end
            offset = offset + ctd.control_dimension
        end
    end

    # variable box
    offset = (N+1) * (ctd.dim_NLP_state + ctd.control_dimension)
    if ctd.has_variable_box
        for j in 1:ctd.dim_variable_box
            indice = ctd.variable_box[2][j]
            l_var[offset+indice] = ctd.variable_box[1][j]
            u_var[offset+indice] = ctd.variable_box[3][j]
        end
    end

    return l_var, u_var
end


# IPOPT objective
function ipopt_objective(xu, ctd)

    t0 = get_initial_time(xu, ctd)
    tf = get_final_time(xu, ctd)
    N = ctd.dim_NLP_steps
    obj = 0
    
    if ctd.has_mayer_cost
        v = get_variable(xu, ctd)
        x0 = get_state_at_time_step(xu, ctd, 0)
        xf = get_state_at_time_step(xu, ctd, N)
        obj = obj + ctd.mayer(x0, xf, v)
    end
    
    if ctd.has_lagrange_cost
        obj = obj + xu[(N+1)*ctd.dim_NLP_state]
    end

    if ctd.criterion_min_max == :min
        return obj
    else
        return -obj
    end
end


# IPOPT constraints  (add bounds computation here at first call ?)
function ipopt_constraint(xu, ctd)
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
    t0 = get_initial_time(xu, ctd)
    tf = get_final_time(xu, ctd)
    N = ctd.dim_NLP_steps
    h = (tf - t0) / N
    c = zeros(eltype(xu), ctd.dim_NLP_constraints)
    v = get_variable(xu, ctd)

    # state equation
    ti = t0
    xi = get_state_at_time_step(xu, ctd, 0)
    ui = get_control_at_time_step(xu, ctd, 0)
    fi = ctd.dynamics(ti, xi, ui, v)
    if ctd.has_lagrange_cost
        xli = get_lagrange_cost_at_time_step(xu, ctd, 0)
        li = ctd.lagrange(ti, xi, ui, v)
    end

    index = 1 # counter for the constraints
    for i in 0:N-1
        tip1 = t0 + (i+1)*h
        
        # state and control at t_{i+1}
        xip1 = get_state_at_time_step(xu, ctd, i+1)
        uip1 = get_control_at_time_step(xu, ctd, i+1)
        fip1 = ctd.dynamics(tip1, xip1, uip1, v)

        # state equation
        if ctd.state_dimension == 1
            c[index] = xip1 - (xi + 0.5*h*(fi + fip1))            
        else
            c[index:index+ctd.state_dimension-1] = xip1 - (xi + 0.5*h*(fi + fip1))
        end
        if ctd.has_lagrange_cost
            xlip1 = get_lagrange_cost_at_time_step(xu, ctd, i+1)
            lip1 = ctd.lagrange(tip1, xip1, uip1, v)
            c[index+ctd.state_dimension] = xlip1 - (xli + 0.5*h*(li + lip1))
            xli = xlip1
            li = lip1
        end
        index = index + ctd.dim_NLP_state

        # path constraints
        if ctd.has_control_constraints
            c[index:index+ctd.dim_control_constraints-1] = ctd.control_constraints[2](ti, ui, v)
            index = index + ctd.dim_control_constraints
        end
        if ctd.has_state_constraints
            c[index:index+ctd.dim_state_constraints-1] = ctd.state_constraints[2](ti, xi ,v)
            index = index + ctd.dim_state_constraints
        end
        if ctd.has_mixed_constraints
            c[index:index+ctd.dim_mixed_constraints-1] = ctd.mixed_constraints[2](ti, xi, ui, v)
            index = index + ctd.dim_mixed_constraints
        end
        xi = xip1
        ui = uip1
        fi = fip1
    end

    # path constraints at final time
    xf = get_state_at_time_step(xu, ctd, N)
    uf = get_control_at_time_step(xu, ctd, N)
    if ctd.has_control_constraints
        c[index:index+ctd.dim_control_constraints-1] = ctd.control_constraints[2](tf, uf, v)      
        index = index + ctd.dim_control_constraints
    end  
    if ctd.has_state_constraints
        c[index:index+ctd.dim_state_constraints-1] = ctd.state_constraints[2](tf, xf, v)      
        index = index + ctd.dim_state_constraints
    end 
    if ctd.has_mixed_constraints
        c[index:index+ctd.dim_mixed_constraints-1] = ctd.mixed_constraints[2](tf, xf, uf, v)
        index = index + ctd.dim_mixed_constraints
    end

    # boundary conditions
    if ctd.has_boundary_conditions
        x0 = get_state_at_time_step(xu, ctd, 0)
        c[index:index+ctd.dim_boundary_conditions-1] = ctd.boundary_conditions[2](x0, xf, v)
        index = index + ctd.dim_boundary_conditions
    end

    # variable constraints
    if ctd.has_variable_constraints
        c[index:index+ctd.dim_variable_constraints-1] = ctd.variable_constraints[2](v)
        index = index + ctd.dim_variable_constraints
    end

    # null initial condition for augmented state (reformulated lagrangian cost)
    if ctd.has_lagrange_cost
        c[index] = get_lagrange_cost_at_time_step(xu, ctd, 0)
        index = index + 1
    end

    return c
end

# struct for ocop/nlp info
mutable struct CTDirect_data

    ## OCP
    # OCP variables and functions
    initial_time
    final_time
    state_dimension
    control_dimension
    dynamics
    mayer
    lagrange
    criterion_min_max
    has_free_final_time
    has_lagrange_cost
    has_mayer_cost

    # OCP constraints
    # indicators
    has_control_constraints
    has_state_constraints
    has_mixed_constraints
    has_boundary_conditions
    has_control_box
    has_state_box

    # dimensions
    dim_control_constraints
    dim_state_constraints
    dim_mixed_constraints
    dim_path_constraints
    dim_boundary_conditions
    dim_control_box
    dim_state_box

    # functions
    control_constraints
    state_constraints
    mixed_constraints
    boundary_conditions
    control_box
    state_box

    ## NLP
    # NLP problem
    dim_NLP_state
    dim_NLP_constraints
    dim_NLP_variables
    dim_NLP_steps
    dynamics_lagrange_to_mayer
    NLP_init

    # NLP solution
    NLP_solution
    NLP_objective
    NLP_sol_constraints
    NLP_constraints_violation
    NLP_iterations
    NLP_stats       # remove later ? type is https://juliasmoothoptimizers.github.io/SolverCore.jl/stable/reference/#SolverCore.GenericExecutionStats


    # put this constructor in CTDirect.jl or in utils.jl ?
    function CTDirect_data(ocp::OptimalControlModel, N::Integer, init=nothing)

        ctd = new()

        ## Optimal Control Problem OCP
        # time
        ctd.initial_time = ocp.initial_time
        ctd.final_time = ocp.final_time
        ctd.has_free_final_time = isnothing(ctd.final_time)

        # dimensions and functions
        ctd.state_dimension = ocp.state_dimension
        ctd.control_dimension = ocp.control_dimension
        ctd.dynamics = ocp.dynamics
        ctd.has_lagrange_cost = !isnothing(ocp.lagrange)
        ctd.lagrange = ocp.lagrange
        ctd.has_mayer_cost = !isnothing(ocp.mayer)
        ctd.mayer = ocp.mayer
        
        # constraints
        ctd.control_constraints, ctd.state_constraints, ctd.mixed_constraints, ctd.boundary_conditions, ctd.control_box, ctd.state_box = nlp_constraints(ocp)
        ctd.dim_control_constraints = length(ctd.control_constraints[1])
        ctd.dim_state_constraints = length(ctd.state_constraints[1])
        ctd.dim_mixed_constraints = length(ctd.mixed_constraints[1])
        ctd.dim_path_constraints = ctd.dim_control_constraints + ctd.dim_state_constraints + ctd.dim_mixed_constraints
        ctd.dim_boundary_conditions = length(ctd.boundary_conditions[1])
        ctd.dim_control_box = length(ctd.control_box[1])
        ctd.dim_state_box = length(ctd.state_box[1])
        ctd.has_control_constraints = !isempty(ctd.control_constraints[1])
        ctd.has_state_constraints = !isempty(ctd.state_constraints[1])
        ctd.has_mixed_constraints = !isempty(ctd.mixed_constraints[1])
        ctd.has_boundary_conditions = !isempty(ctd.boundary_conditions[1])
        ctd.has_control_box = !isempty(ctd.control_box[1])
        ctd.has_state_box = !isempty(ctd.state_box[1])

#=         println("ctd constructor after call to nlp_constraints")
        println("state constraints")
        println(ctd.state_constraints)
        println(ctd.state_constraints[1])
        println(ctd.dim_state_constraints)
        println("control constraints")
        println(ctd.control_constraints)
        println(ctd.control_constraints[1])
        println(ctd.dim_control_constraints)
        println("mixed constraints")
        println(ctd.mixed_constraints)
        println(ctd.mixed_constraints[1])
        println(ctd.dim_mixed_constraints)
        println("control box")
        println(ctd.control_box)
        println(ctd.control_box[1])
        println(ctd.dim_control_box)
        println("state box")
        println(ctd.state_box)
        println(ctd.state_box[1])
        println(ctd.dim_state_box) =#


        ## Non Linear Programming NLP
        ctd.dim_NLP_steps = N
        ctd.NLP_init = init

        # Mayer to Lagrange reformulation: 
        # additional state with Lagrange cost as dynamics and null initial condition
        if ctd.has_lagrange_cost
            ctd.dim_NLP_state = ctd.state_dimension + 1  
            ctd.dim_NLP_constraints = N * (ctd.dim_NLP_state + ctd.dim_path_constraints) +
            ctd.dim_path_constraints + ctd.dim_boundary_conditions + 1           
        else
            ctd.dim_NLP_state = ctd.state_dimension  
            ctd.dim_NLP_constraints = N * (ctd.dim_NLP_state + ctd.dim_path_constraints) +
            ctd.dim_path_constraints + ctd.dim_boundary_conditions
        end
        # augmented dynamics (try to evaluate the condition only once cf below)
        #ctd.dynamics_lagrange_to_mayer(t, x, u) = ctd.has_lagrange_cost ? [ctd.dynamics(t, x[1:ctd.state_dimension], u); ctd.lagrange(t, x[1:ctd.state_dimension], u)] : ctd.dynamics(t, x, u) DOES NOT COMPILE
        function f(t, x, u)
            if ctd.has_lagrange_cost
                return [ctd.dynamics(t, x[1:ctd.state_dimension], u); ctd.lagrange(t, x[1:ctd.state_dimension], u)]
            else
                return ctd.dynamics(t, x, u)
            end
        end
        ctd.dynamics_lagrange_to_mayer = f

        # min or max problem (unused ?)
        ctd.criterion_min_max = ocp.criterion

        # additional variable for free final time
        if ctd.has_free_final_time
            ctd.dim_NLP_variables = (N + 1) * (ctd.dim_NLP_state + ctd.control_dimension) + 1
        else
            ctd.dim_NLP_variables = (N + 1) * (ctd.dim_NLP_state + ctd.control_dimension)
        end

        return ctd

    end

end

function is_solvable(ocp)
    solvable = true

    # free initial time
    if isnothing(ocp.initial_time)
        solvable = false
    end
    
    return solvable
end


# bounds for the constraints
function  constraints_bounds(ctd)

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
    lb[index:index+ctd.dim_boundary_conditions-1] = ctd.boundary_conditions[1]
    ub[index:index+ctd.dim_boundary_conditions-1] = ctd.boundary_conditions[3]
    index = index + ctd.dim_boundary_conditions
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
    
    # NLP variables layout: [X0, X1 .. XN, U0, U1 .. UN]

    # state box
    if ctd.has_state_box
        index = 0
        for i in 0:N
            for j in 1:ctd.dim_state_box
                indice = ctd.state_box[2][j]
                #= l_var[index+indice] = ctd.state_box[1][indice]
                u_var[index+indice] = ctd.state_box[3][indice] =#
                l_var[index+indice] = ctd.state_box[1][j]
                u_var[index+indice] = ctd.state_box[3][j]
            end
            index = index + ctd.dim_NLP_state
        end
    end

    # control box
    if ctd.has_control_box
        index = (N+1)*ctd.dim_NLP_state 
        for i in 0:N
            for j in 1:ctd.dim_control_box
                indice = ctd.control_box[2][j]
                l_var[index+indice] = ctd.control_box[1][j]
                u_var[index+indice] = ctd.control_box[3][j]
            end
            index = index + ctd.control_dimension
        end
    end

    # free final time case
    if ctd.has_free_final_time
        l_var[end] = 1.e-3
    end

    return l_var, u_var
end


# IPOPT objective
function ipopt_objective(xu, ctd)

    t0 = ctd.initial_time
    tf = get_final_time(xu, ctd.final_time, ctd.has_free_final_time)
    N = ctd.dim_NLP_steps
    obj = 0
    
    if ctd.has_mayer_cost
        x0 = get_state_at_time_step(xu, 0, ctd.dim_NLP_state, N)
        xf = get_state_at_time_step(xu, N, ctd.dim_NLP_state, N)
        obj = obj + ctd.mayer(t0, x0[1:ctd.state_dimension], tf, xf[1:ctd.state_dimension])
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


# IPOPT constraints +++ add bounds computation here at first call
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
        u_m(t_N), ..., u_m(t_N)]
    return
    c :: 
    """
    t0 = ctd.initial_time
    tf = get_final_time(xu, ctd.final_time, ctd.has_free_final_time)
    N = ctd.dim_NLP_steps
    h = (tf - t0) / N
    c = zeros(eltype(xu), ctd.dim_NLP_constraints)

    # state equation
    ti = t0
    xi = get_state_at_time_step(xu, 0, ctd.dim_NLP_state, N)
    ui = get_control_at_time_step(xu, 0, ctd.dim_NLP_state, N, ctd.control_dimension)
    fi = ctd.dynamics_lagrange_to_mayer(ti, xi, ui)
    index = 1 # counter for the constraints
    for i in 0:N-1
        tip1 = t0 + (i+1)*h
        # state and control at t_{i+1}
        xip1 = get_state_at_time_step(xu, i+1, ctd.dim_NLP_state, N)
        uip1 = get_control_at_time_step(xu, i+1, ctd.dim_NLP_state, N, ctd.control_dimension)
        fip1 = ctd.dynamics_lagrange_to_mayer(tip1, xip1, uip1)
        # state equation
        c[index:index+ctd.dim_NLP_state-1] = xip1 - (xi + 0.5*h*(fi + fip1))
        index = index + ctd.dim_NLP_state

        # path constraints
        if ctd.has_control_constraints
            c[index:index+ctd.dim_control_constraints-1] = ctd.control_constraints[2](ti, ui)        # ui vector
            index = index + ctd.dim_control_constraints
        end
        if ctd.has_state_constraints
            c[index:index+ctd.dim_state_constraints-1] = ctd.state_constraints[2](ti, xi[1:ctd.state_dimension])
            index = index + ctd.dim_state_constraints
        end
        if ctd.has_mixed_constraints
            c[index:index+ctd.dim_mixed_constraints-1] = ctd.mixed_constraints[2](ti, xi[1:ctd.state_dimension], ui)
            index = index + ctd.dim_mixed_constraints
        end
        xi = xip1
        ui = uip1
        fi = fip1
    end

    # path constraints at final time
    if ctd.has_control_constraints
        uf = get_control_at_time_step(xu, N, ctd.dim_NLP_state, N, ctd.control_dimension)
        c[index:index+ctd.dim_control_constraints-1] = ctd.control_constraints[2](tf, uf)      
        index = index + ctd.dim_control_constraints
    end  
    if ctd.has_state_constraints
        xf = get_state_at_time_step(xu, N, ctd.dim_NLP_state, N)
        c[index:index+ctd.dim_state_constraints-1] = ctd.state_constraints[2](tf, xf[1:ctd.state_dimension])      
        index = index + ctd.dim_state_constraints
    end 
    if ctd.has_mixed_constraints
        xf = get_state_at_time_step(xu, N, ctd.dim_NLP_state, N)
        uf = get_control_at_time_step(xu, N-1, ctd.dim_NLP_state, N, ctd.control_dimension)
        c[index:index+ctd.dim_mixed_constraints-1] = ctd.mixed_constraints[2](tf, xf[1:ctd.state_dimension], uf)
        index = index + ctd.dim_mixed_constraints
    end

    # boundary conditions
    x0 = get_state_at_time_step(xu, 0, ctd.dim_NLP_state, N)
    xf = get_state_at_time_step(xu, N, ctd.dim_NLP_state, N)
    c[index:index+ctd.dim_boundary_conditions-1] = ctd.boundary_conditions[2](t0, x0[1:ctd.state_dimension], tf, xf[1:ctd.state_dimension])
    index = index + ctd.dim_boundary_conditions
    # null initial condition for augmented state (reformulated lagrangian cost)
    if ctd.has_lagrange_cost
        c[index] = xu[ctd.dim_NLP_state]
        index = index + 1
    end

    return c
end

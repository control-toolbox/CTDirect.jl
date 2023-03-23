function ADNLProblem(ocp::OptimalControlModel, N::Integer, init=nothing)

    ctd = CTDirect_data(ocp, N)

    # IPOPT objective
    function ipopt_objective(xu)
        t0 = ctd.initial_time
        tf = get_final_time(xu, ctd.final_time, ctd.has_free_final_time)
        obj = 0
        if ctd.has_mayer_cost
            x0 = get_state_at_time_step(xu, 0, ctd.dim_NLP_state, N)
            xf = get_state_at_time_step(xu, N, ctd.dim_NLP_state, N)
            obj = obj + ctd.mayer(t0, x0[1:ctd.state_dimension], tf, xf[1:ctd.state_dimension])
        end
        if ctd.has_lagrange_cost
            obj = obj + xu[(N+1)*ctd.dim_NLP_state]
        end
        return ismin(ocp) ? obj : -obj
    end

    # IPOPT constraints
    function ipopt_constraint(xu)
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

    # bounds for the constraints
    function  constraints_bounds()
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

    # box constraints
    function variables_bounds()

        l_var = -Inf*ones(ctd.dim_NLP_variables)
        u_var = Inf*ones(ctd.dim_NLP_variables)
        
        # NLP variables layout: [X0, X1 .. XN, U0, U1 .. UN]

        # state box
        if ctd.has_state_box
            index = 0
            for i in 0:N
                for j in 1:ctd.dim_state_box
                    indice = ctd.state_box[2][j]
                    l_var[index+indice] = ctd.state_box[1][indice]
                    u_var[index+indice] = ctd.state_box[3][indice]
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
                    l_var[index+indice] = ctd.control_box[1][indice]
                    u_var[index+indice] = ctd.control_box[3][indice]
                end
                index = index + ctd.control_dimension
            end
        end
        return l_var, u_var
    end

    # generate initial guess
    function set_state_at_time_step!(x, i, nx, N, xu)
        @assert i <= N "trying to set init for x(t_i) with i > N"
        xu[1+i*nx:(i+1)*nx] = x[1:nx]
    end
    
    function set_control_at_time_step!(u, i, nx, N, m, xu)
        @assert i <= N "trying to set init for u(t_i) with i > N"
        xu[1+(N+1)*nx+i*m:m+(N+1)*nx+i*m] = u[1:m]
    end

    function initial_guess()

        if init === nothing
            # default initialization (put back O.1 ?)
            xu0 = 1.1*ones(ctd.dim_NLP_variables)
        else
            if length(init) != (ctd.state_dimension + ctd.control_dimension)
                error("vector for initialization should be of size n+m",ctd.state_dimension+ctd.control_dimension)
            end
            # split state / control init values
            x_init = zeros(ctd.dim_NLP_state)
            x_init[1:ctd.state_dimension] = init[1:ctd.state_dimension]
            u_init = zeros(ctd.control_dimension)
            u_init[1:ctd.control_dimension] = init[ctd.state_dimension+1:ctd.state_dimension+ctd.control_dimension]
            
            # mayer -> lagrange additional state
            if ctd.has_lagrange_cost
                x_init[ctd.dim_NLP_state] = 0.1
            end

            # set constant initialization for state / control variables
            xu0 = zeros(ctd.dim_NLP_variables)
            for i in 0:N
                set_state_at_time_step!(x_init, i, ctd.dim_NLP_state, N, xu0)
                set_control_at_time_step!(u_init, i, ctd.dim_NLP_state, N, ctd.control_dimension, xu0)
            end
        end
        return xu0
    end

    # variables bounds   
    l_var, u_var = variables_bounds()

    # initial guess
    xu0 = initial_guess()

    # free final time case
    if ctd.has_free_final_time
      xu0[end] = 1.0
      l_var[end] = 1.e-3
    end

    lb, ub = constraints_bounds()

    nlp = ADNLPModel(ipopt_objective, xu0, l_var, u_var, ipopt_constraint, lb, ub)    

    return nlp

end

function ADNLProblem(ocp::OptimalControlModel, N::Integer, init=nothing)

    # direct_infos
    t0, tf_, n_x, m, f, control_constraints, state_constraints, mixed_constraints, boundary_conditions, control_box, state_box, dim_control_constraints, dim_state_constraints, dim_mixed_constraints, dim_boundary_conditions, dim_control_box, dim_state_box,has_control_constraints, has_state_constraints, has_mixed_constraints, has_boundary_conditions, has_control_box, has_state_box, hasLagrangeCost, hasMayerCost, dim_x, nc, dim_xu, g, f_Mayer, has_free_final_time, criterion = direct_infos(ocp, N)

    # IPOPT objective
    function ipopt_objective(xu)
        tf = get_final_time(xu, tf_, has_free_final_time)
        obj = 0
        if hasMayerCost
            x0 = get_state_at_time_step(xu, 0, dim_x, N)
            xf = get_state_at_time_step(xu, N, dim_x, N)
            obj = obj + g(t0, x0[1:n_x], tf, xf[1:n_x])
        end
        if hasLagrangeCost
            obj = obj + xu[(N+1)*dim_x]
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
        tf = get_final_time(xu, tf_, has_free_final_time)
        h = (tf-t0)/N
        c = zeros(eltype(xu),nc)
        #

        # state equation
        ti = t0
        xi = get_state_at_time_step(xu, 0, dim_x, N)
        ui = get_control_at_time_step(xu, 0, dim_x, N, m)
        fi = f_Mayer(ti, xi, ui)
        index = 1 # counter for the constraints
        for i in 0:N-1
            tip1 = t0 + (i+1)*h
            # state and control at t_{i+1}
            xip1 = get_state_at_time_step(xu, i+1, dim_x, N)
            uip1 = get_control_at_time_step(xu, i+1, dim_x, N, m)
            fip1 = f_Mayer(tip1, xip1, uip1)
            # state equation
            c[index:index+dim_x-1] = xip1 - (xi + 0.5*h*(fi + fip1))
            index = index + dim_x

            # control and mixed constraints
            if has_control_constraints
                c[index:index+dim_control_constraints-1] = control_constraints[2](ti, ui)        # ui vector
                index = index + dim_control_constraints
            end

            if has_state_constraints
                c[index:index+dim_state_constraints-1] = state_constraints[2](ti, xi[1:n_x])
                index = index + dim_state_constraints
            end

            if has_mixed_constraints
                c[index:index+dim_mixed_constraints-1] = mixed_constraints[2](ti, xi[1:n_x], ui)        # ui vector
                index = index + dim_mixed_constraints
            end
            xi = xip1
            ui = uip1
            fi = fip1
        end
        if has_control_constraints
            uf = get_control_at_time_step(xu, N, dim_x, N, m)
            c[index:index+dim_control_constraints-1] = control_constraints[2](tf, uf)      
            index = index + dim_control_constraints
        end  

        if has_state_constraints
            xf = get_state_at_time_step(xu, N, dim_x, N)
            c[index:index+dim_state_constraints-1] = state_constraints[2](tf, xf[1:n_x])      
            index = index + dim_state_constraints
        end 

        if has_mixed_constraints
            xf = get_state_at_time_step(xu, N, dim_x, N)
            uf = get_control_at_time_step(xu, N-1, dim_x, N, m)
            c[index:index+dim_mixed_constraints-1] = mixed_constraints[2](tf, xf[1:n_x], uf)        # ui is false because Euler
            index = index + dim_mixed_constraints
        end

        # boundary conditions
        # -------------------
        x0 = get_state_at_time_step(xu, 0, dim_x, N)
        xf = get_state_at_time_step(xu, N, dim_x, N)
        c[index:index+dim_boundary_conditions-1] = boundary_conditions[2](t0, x0[1:n_x], tf, xf[1:n_x])  # because Lagrange cost possible
        index = index + dim_boundary_conditions
        if hasLagrangeCost
            c[index] = xu[dim_x]
            index = index + 1
        end

        return c
    end

    # bounds for the constraints
    function  constraints_bounds()
        lb = zeros(nc)
        ub = zeros(nc)
        index = 1 # counter for the constraints
        for i in 0:N-1
            index = index + dim_x          # leave 0 for the state equation
            if has_control_constraints
                lb[index:index+dim_control_constraints-1] = control_constraints[1]
                ub[index:index+dim_control_constraints-1] = control_constraints[3]
                index = index + dim_control_constraints
            end
            if has_state_constraints
                lb[index:index+dim_state_constraints-1] = state_constraints[1]
                ub[index:index+dim_state_constraints-1] = state_constraints[3]
                index = index + dim_state_constraints
            end
            if has_mixed_constraints
                lb[index:index+dim_mixed_constraints-1] = mixed_constraints[1]
                ub[index:index+dim_mixed_constraints-1] = mixed_constraints[3]
                index = index + dim_mixed_constraints
            end
        end
        if has_control_constraints
            lb[index:index+dim_control_constraints-1] = control_constraints[1]
            ub[index:index+dim_control_constraints-1] = control_constraints[3]
            index = index + dim_control_constraints
        end
        if has_state_constraints
            lb[index:index+dim_state_constraints-1] = state_constraints[1]
            ub[index:index+dim_state_constraints-1] = state_constraints[3]
            index = index + dim_state_constraints
        end
        if has_mixed_constraints
            lb[index:index+dim_mixed_constraints-1] = mixed_constraints[1]
            ub[index:index+dim_mixed_constraints-1] = mixed_constraints[3]
            index = index + dim_mixed_constraints
        end
        # boundary conditions
        lb[index:index+dim_boundary_conditions-1] = boundary_conditions[1]
        ub[index:index+dim_boundary_conditions-1] = boundary_conditions[3]
        index = index + dim_boundary_conditions
        if hasLagrangeCost
            lb[index] = 0.
            ub[index] = 0.
            index = index + 1
        end

        return lb, ub
    end

    # todo: retrieve optional bounds from ocp parsed constraints
    function variables_bounds()
        # unbounded case
        l_var = -Inf*ones(dim_xu)
        u_var = Inf*ones(dim_xu)
        
        # state box
        if has_state_box
            index = 0
            for i in 0:N
                for j in 1:dim_state_box
                    indice = state_box[2][j]
                    l_var[index+indice] = state_box[1][indice]
                    u_var[index+indice] = state_box[3][indice]
                end
                index = index + dim_x
            end
        end
        # control box
        if has_control_box
            index = (N+1)*dim_x # the control is at the  end of xu
            for i in 0:N
                for j in 1:dim_control_box
                    indice = control_box[2][j]
                    l_var[index+indice] = control_box[1][indice]
                    u_var[index+indice] = control_box[3][control_box[2][j]]
                end
                index = index + m
            end
        end
        return l_var, u_var
    end

    # generate initial guess
    function set_state_at_time_step!(x, i, dim_x, N, xu)
        if i > N
            error("trying to set x(t_i) for i > N")
        else
            xu[1+i*dim_x:(i+1)*dim_x] = x[1:dim_x]
        end
    end
    
    function set_control_at_time_step!(u, i, dim_x, N, m, xu)
        if i > N
            error("trying to set (t_i) for i > N")
        else
            xu[1+(N+1)*dim_x+i*m:m+(N+1)*dim_x+i*m] = u[1:m]
        end
    end

    function initial_guess()
        #println("Initialization: ", init)

        if init === nothing
            # default initialization
            xu0 = 1.1*ones(dim_xu)
        else
            if length(init) != (n_x + m)
                error("vector for initialization should be of size n+m",n_x+m)
            end
            # split state / control
            x_init = zeros(dim_x)
            x_init[1:n_x] = init[1:n_x]
            u_init = zeros(m)
            u_init[1:m] = init[n_x+1:n_x+m]
            
            # mayer -> lagrange additional state
            if hasLagrangeCost
                x_init[dim_x] = 0.1
            end

            # constant initialization
            xu0 = zeros(dim_xu)
            for i in 0:N
                set_state_at_time_step!(x_init, i, dim_x, N, xu0)
                set_control_at_time_step!(u_init, i, dim_x, N, m, xu0)
            end
        end
        return xu0
    end

    # variables bounds   
    l_var, u_var = variables_bounds()
    println("lvar = ", l_var)
    println("uvar = ", u_var)
    # initial guess
    xu0 = initial_guess()

    # free final time case
    if has_free_final_time
      xu0[end] = 1.0
      l_var[end] = 1.e-3
    end

    lb, ub = constraints_bounds()

    nlp = ADNLPModel(ipopt_objective, xu0, l_var, u_var, ipopt_constraint, lb, ub)    

    return nlp

end

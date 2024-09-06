# Discretized Optimal Control Problem DOCP
# Notes:
# - for now the path constraints are checked on the time steps ie g(t_i, x_i, u_i). This requires a getter that provide a control value at each time step, which may not coincide with the actual control discretization (eg RK stages). In this case we will use the 'average' control using the same coefficients as the RK method. Later we can add an option for control discretization, step or stage. Taking a control constant per step would solve the question for the path constraints evaluation and reduces the number of variables. Compared to the more standard stage control discretization, the consistency of the costate would need to be checked, as well as the oscillations in the trajectory (which may actually be better). Further options may include CVP (control vector parametrization) on a coarser grid.
# - for the other choice of enforcing path constraints at the time stages, the symmetric question of getting state values at time stages can be solved by reusing the states used in the evaluation of the stage dynamics. This second choice (as in Bocop2) has the drawback of a larger problem size, and does not check the constraints at the points of the actual trajectory computed (including tf).

# generic discretization
abstract type Discretization end
abstract type ArgsAtStep end

"""
$(TYPEDSIGNATURES)

Struct for discretized optimal control problem DOCP

Contains:
- a copy of the original OCP
- data required to link the OCP with the discretized DOCP
"""
struct DOCP{T <: Discretization}

    ## OCP
    ocp::OptimalControlModel
    control_constraints::Any
    state_constraints::Any
    mixed_constraints::Any
    boundary_constraints::Any
    variable_constraints::Any
    control_box::Any
    state_box::Any
    variable_box::Any

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
    NLP_time_grid::Vector{Any}
    dim_NLP_variables::Int
    dim_NLP_constraints::Int

    # lower and upper bounds for variables and constraints
    var_l::Vector{Float64}
    var_u::Vector{Float64}
    con_l::Vector{Float64}
    con_u::Vector{Float64}

    # discretization scheme
    discretization::T

    # constructor
    function DOCP(ocp::OptimalControlModel, grid_size::Integer, time_grid, discretization::Discretization)

        # time grid
        if time_grid == nothing
            NLP_normalized_time_grid = collect(LinRange(0, 1, grid_size + 1))
            dim_NLP_steps = grid_size
        else
            # check strictly increasing
            if !issorted(time_grid, lt = <=)
                throw(ArgumentError("given time grid is not strictly increasing. Aborting..."))
                return nothing
            end
            # normalize input grid if needed
            if (time_grid[1] != 0) || (time_grid[end] != 1)
                #println("INFO: normalizing given time grid...")
                t0 = time_grid[1]
                tf = time_grid[end]
                NLP_normalized_time_grid = (time_grid .- t0) ./ (tf - t0)
            else
                NLP_normalized_time_grid = time_grid
            end
            dim_NLP_steps = length(time_grid) - 1
        end

        # additional indicators (faster ?)
        has_free_t0 = has_free_initial_time(ocp)
        has_free_tf = has_free_final_time(ocp)
        has_lagrange = has_lagrange_cost(ocp)
        has_mayer = has_mayer_cost(ocp)
        has_variable = is_variable_dependent(ocp)
        has_maximization = is_max(ocp)

        if has_free_t0 || has_free_tf 
            NLP_time_grid = Vector{Any}(undef, dim_NLP_steps+1)
        else 
            NLP_time_grid = ocp.initial_time .+ (NLP_normalized_time_grid .* (ocp.final_time - ocp.initial_time))
        end

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
        dim_stage = discretization.stage

        # NLP unknown (state + control + variable [+ stage])
        dim_NLP_variables = (N + 1) * dim_NLP_x + (N + discretization.additional_controls) * dim_NLP_u + dim_NLP_v + N * dim_NLP_x * dim_stage

        # NLP constraints 
        # parse NLP constraints (and initialize dimensions)
        control_constraints,
        state_constraints,
        mixed_constraints,
        boundary_constraints,
        variable_constraints,
        control_box,
        state_box,
        variable_box = CTBase.nlp_constraints!(ocp)

        dim_x_box = dim_state_range(ocp)
        dim_u_box = dim_control_range(ocp)
        dim_v_box = dim_variable_range(ocp)
        dim_path_cons = dim_path_constraints(ocp)
        dim_x_cons = dim_state_constraints(ocp)
        dim_u_cons = dim_control_constraints(ocp)
        dim_v_cons = dim_variable_constraints(ocp)
        dim_mixed_cons = dim_mixed_constraints(ocp)
        dim_boundary_cons = dim_boundary_constraints(ocp)

        # constraints (dynamics, stage, path, boundary, variable)
        dim_NLP_constraints =
            N * (dim_NLP_x + (dim_NLP_x * dim_stage) + dim_path_cons) + dim_path_cons + dim_boundary_cons + dim_v_cons
        if has_lagrange
            # add initial condition for lagrange state
            dim_NLP_constraints += 1
        end

        # call constructor with const fields
        docp = new{typeof(discretization)}(
            ocp,
            control_constraints,
            state_constraints,
            mixed_constraints,
            boundary_constraints,
            variable_constraints,
            control_box,
            state_box,
            variable_box,
            has_free_t0,
            has_free_tf,
            has_lagrange,
            has_mayer,
            has_variable,
            has_maximization,
            dim_x_box,
            dim_u_box,
            dim_v_box,
            dim_path_cons,
            dim_x_cons,
            dim_u_cons,
            dim_v_cons,
            dim_mixed_cons,
            dim_boundary_cons,
            dim_NLP_x,
            dim_NLP_u,
            dim_NLP_v,
            dim_OCP_x,
            dim_NLP_steps,
            NLP_normalized_time_grid,
            NLP_time_grid,
            dim_NLP_variables,
            dim_NLP_constraints,
            -Inf * ones(dim_NLP_variables),
            Inf * ones(dim_NLP_variables),
            zeros(dim_NLP_constraints),
            zeros(dim_NLP_constraints),
            discretization
        )

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
    for i = 0:(docp.dim_NLP_steps - 1)
        # skip (ie leave 0) for equality dynamics constraint
        index = index + docp.dim_NLP_x
        # skip (ie leave 0) for equality stage constraint (ki)
        index = index + docp.dim_NLP_x * docp.discretization.stage
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

    # first we build full ordered sets of bounds, then set them in NLP

    # state / control box
    x_lb, x_ub = build_bounds(docp.dim_OCP_x, docp.dim_x_box, docp.state_box)
    u_lb, u_ub = build_bounds(docp.dim_NLP_u, docp.dim_u_box, docp.control_box)
    for i = 0:N
        set_variables_at_t_i!(var_l, x_lb, u_lb, docp, i)
        set_variables_at_t_i!(var_u, x_ub, u_ub, docp, i)
    end

    # variable box
    if docp.has_variable
        v_lb, v_ub = build_bounds(docp.dim_NLP_v, docp.dim_v_box, docp.variable_box)
        set_optim_variable!(var_l, v_lb, docp)
        set_optim_variable!(var_u, v_ub, docp)
    end

    return var_l, var_u
end


"""
$(TYPEDSIGNATURES)

Compute the objective for the DOCP problem.
"""
function DOCP_objective(xu, docp::DOCP)

    obj = 0.0
    N = docp.dim_NLP_steps
    ocp = docp.ocp

    # optimization variables
    v = Float64[]
    docp.has_variable && (v = get_optim_variable(xu, docp))

    # final state is always needed since lagrange cost is there
    xf, uf, xlf = get_variables_at_time_step(xu, docp, N+1)

    # mayer cost
    if docp.has_mayer
        x0, u0, xl0 = get_variables_at_time_step(xu, docp, 1)
        obj = obj + ocp.mayer(x0, xf, v)
    end

    # lagrange cost
    if docp.has_lagrange
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

    N = docp.dim_NLP_steps
    if docp.has_free_t0 || docp.has_free_tf
        get_time_grid!(docp.NLP_time_grid, xu, docp)
    end

    v = Float64[]
    docp.has_variable && (v = get_optim_variable(xu, docp))

    # main loop on time steps 
    for i = 1:N
        setConstraintBlock!(docp, c, xu, v, docp.NLP_time_grid, i)
    end

    # path constraints at final time
    offset = N * (docp.dim_NLP_x*(1+docp.discretization.stage) + docp.dim_path_cons)
    tf = docp.NLP_time_grid[N+1]
    xf, uf = get_variables_at_time_step(xu, docp, N+1)
    setPathConstraints!(docp, c, tf, xf, uf, v, offset)

    # point constraints
    setPointConstraints!(docp, c, xu, v)

    return c # needed even for inplace version, AD error otherwise
end


"""
$(TYPEDSIGNATURES)

Set the path constraints at given time step: [control, state, mixed]
Convention: 1 <= i <= dim_NLP_steps+1
"""
function setPathConstraints!(docp::DOCP, c, t_i, x_i, u_i, v, offset)    

    # NB. using .= below *doubles* the allocations oO +++ later inplace
    if docp.dim_u_cons > 0
        c[offset+1:offset+docp.dim_u_cons] = docp.control_constraints[2](t_i, u_i, v)
    end
    if docp.dim_x_cons > 0 
        c[offset+docp.dim_u_cons+1:offset+docp.dim_u_cons+docp.dim_x_cons] = docp.state_constraints[2](t_i, x_i, v)
    end
    if docp.dim_mixed_cons > 0 
        c[offset+docp.dim_u_cons+docp.dim_x_cons+1:offset+docp.dim_u_cons+docp.dim_x_cons+docp.dim_mixed_cons] = docp.mixed_constraints[2](t_i, x_i, u_i, v)
    end

end


"""
$(TYPEDSIGNATURES)

Set bounds for the path constraints at given time step
"""
function setPathBounds!(docp::DOCP, index::Int, lb, ub)

    # pure control constraints
    if docp.dim_u_cons > 0
        lb[index:(index + docp.dim_u_cons - 1)] = docp.control_constraints[1]
        ub[index:(index + docp.dim_u_cons - 1)] = docp.control_constraints[3]
        index = index + docp.dim_u_cons
    end

    # pure state constraints
    if docp.dim_x_cons > 0
        lb[index:(index + docp.dim_x_cons - 1)] = docp.state_constraints[1]
        ub[index:(index + docp.dim_x_cons - 1)] = docp.state_constraints[3]
        index = index + docp.dim_x_cons
    end

    # mixed state / control constraints
    if docp.dim_mixed_cons > 0
        lb[index:(index + docp.dim_mixed_cons - 1)] = docp.mixed_constraints[1]
        ub[index:(index + docp.dim_mixed_cons - 1)] = docp.mixed_constraints[3]
        index = index + docp.dim_mixed_cons
    end

    return index
end


"""
$(TYPEDSIGNATURES)

Set the boundary and variable constraints
"""
function setPointConstraints!(docp::DOCP, c, xu, v)

    # offset
    offset = docp.dim_NLP_steps * (docp.dim_NLP_x * (1+docp.discretization.stage) + docp.dim_path_cons) + docp.dim_path_cons

    # variables
    x0, _, xl0 = get_variables_at_time_step(xu, docp, 1)
    xf, = get_variables_at_time_step(xu, docp, docp.dim_NLP_steps+1)

    # boundary constraints
    if docp.dim_boundary_cons > 0
        c[offset+1:offset+docp.dim_boundary_cons] = docp.boundary_constraints[2](x0, xf, v)
    end

    # variable constraints
    if docp.dim_v_cons > 0
        c[offset+docp.dim_boundary_cons+1:offset+docp.dim_boundary_cons+docp.dim_v_cons] = docp.variable_constraints[2](v)
    end

    # null initial condition for lagrangian cost state
    if docp.has_lagrange
        c[offset+docp.dim_boundary_cons+docp.dim_v_cons+1] = xl0
    end
end


"""
$(TYPEDSIGNATURES)

Set bounds for the boundary and variable constraints
"""
function setPointBounds!(docp::DOCP, index::Int, lb, ub)

    ocp = docp.ocp

    # boundary constraints
    if docp.dim_boundary_cons > 0
        lb[index:(index + docp.dim_boundary_cons - 1)] = docp.boundary_constraints[1]
        ub[index:(index + docp.dim_boundary_cons - 1)] = docp.boundary_constraints[3]
        index = index + docp.dim_boundary_cons
    end

    # variable constraints
    if docp.dim_v_cons > 0
        lb[index:(index + docp.dim_v_cons - 1)] = docp.variable_constraints[1]
        ub[index:(index + docp.dim_v_cons - 1)] = docp.variable_constraints[3]
        index = index + docp.dim_v_cons
    end

    # null initial condition for lagrangian cost state
    if docp.has_lagrange
        lb[index] = 0.0
        ub[index] = 0.0
        index = index + 1
    end

    return index
end


"""
$(TYPEDSIGNATURES)

Build initial guess for discretized problem
"""
function DOCP_initial_guess(docp::DOCP, init::OptimalControlInit = OptimalControlInit())

    # default initialization (internal variables such as lagrange cost, k_i for RK schemes) will keep these default values 
    NLP_X = 0.1 * ones(docp.dim_NLP_variables)

    # set variables if provided (needed first in case of free times !)
    if !isnothing(init.variable_init)
        set_optim_variable!(NLP_X, init.variable_init, docp)
    end

    # set state / control variables if provided
    time_grid = get_time_grid(NLP_X, docp)
    for i = 0:(docp.dim_NLP_steps)
        ti = time_grid[i+1]
        set_variables_at_t_i!(NLP_X, init.state_init(ti), init.control_init(ti), docp, i)
    end

    return NLP_X
end

#= OLD

    # recompute value of constraints at solution
    # NB. the constraint formulation is LB <= C <= UB
    constraints = zeros(docp.dim_NLP_constraints)
    DOCP_constraints!(constraints, solution, docp)
    # set constraint violation if needed
    # is not saved in OCP solution currently...
    if constraints_violation==nothing
        constraints_check = zeros(docp.dim_NLP_constraints)
        DOCP_constraints_check!(constraints_check, constraints, docp)
        println("Recomputed constraints violation ", norm(constraints_check, Inf))
        variables_check = zeros(docp.dim_NLP_variables)
        DOCP_variables_check!(variables_check, solution, docp)
        println("Recomputed variable bounds violation ", norm(variables_check, Inf))
        constraints_violation = norm(append!(variables_check, constraints_check), Inf)

    end

"""
$(TYPEDSIGNATURES)

Check the nonlinear constraints violation for the DOCP problem. 
"""
function DOCP_constraints_check!(cb, constraints, docp)

    # todo add a single utils function check_bounds(v,lb,ub) that returns the error vector ?

    # check constraints vs bounds
    # by construction only one of the two can be active
    for i = 1:(docp.dim_NLP_constraints)
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
    for i = 1:(docp.dim_NLP_variables)
        if variables[i] < docp.var_l[i]
            vb[i] = variables[i] - docp.var_l[i] # < 0
        end
        if variables[i] > docp.var_u[i]
            vb[i] = variables[i] - docp.var_u[i] # > 0
        end
    end
    return nothing
end
=#

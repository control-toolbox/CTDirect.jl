# todo: add more discretization schemes
# use a scheme struct to allow multiple dispatch (cf solve)
# this struct can contains useful info about the scheme
# (name, order, stage, properties etc)
# later we can add an option for control discretization: step or stage
# (meaningful only for schemes with more than 1 stage, at least I will be able to compare the two ! further options may include CVP -control vector parametrization-, and maybe even pseudo-spectral ?)


# generic discretization struct
# NB. can we mutualize common fields at the abstract level ?
abstract type DiscretizationTag end
struct TrapezeTag <: DiscretizationTag 
    stage::Int
    TrapezeTag() = new(0)
end


"""
$(TYPEDSIGNATURES)

Struct for discretized optimal control problem DOCP

Contains:
- a copy of the original OCP
- data required to link the OCP with the discretized DOCP
"""
struct DOCP

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
    dim_NLP_variables::Int
    dim_NLP_constraints::Int

    # lower and upper bounds for variables and constraints
    var_l::Vector{Float64}
    var_u::Vector{Float64}
    con_l::Vector{Float64}
    con_u::Vector{Float64}

    # discretization scheme
    discretization::DiscretizationTag

    # constructor
    function DOCP(ocp::OptimalControlModel, grid_size::Integer, time_grid, discretization::DiscretizationTag)

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
        dim_NLP_variables = (N + 1) * dim_NLP_x + N * dim_NLP_u + dim_NLP_v + N * dim_NLP_x * dim_stage

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
        docp = new(
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

    # state / control box
    x_lb, x_ub = build_bounds(docp.dim_OCP_x, docp.dim_x_box, docp.state_box)
    u_lb, u_ub = build_bounds(docp.dim_NLP_u, docp.dim_u_box, docp.control_box)
    for i = 0:N
        set_variables_at_time_step!(var_l, x_lb, u_lb, docp, i)
        set_variables_at_time_step!(var_u, x_ub, u_ub, docp, i)
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
# DOCP objective
function DOCP_objective(xu, docp::DOCP)

    obj = 0.0
    N = docp.dim_NLP_steps
    ocp = docp.ocp

    # final state is always needed since lagrange cost is there
    xf, uf, xlf = get_variables_at_time_step(xu, docp, N)

    # mayer cost
    if docp.has_mayer
        v = get_optim_variable(xu, docp)
        x0, u0, xl0 = get_variables_at_time_step(xu, docp, 0)
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

    # main loop on time steps
    index = 1 # counter for the constraints
    for i = 0:(docp.dim_NLP_steps - 1)

        # state equation
        index = setStateEquation!(docp, c, index, xu, i)
        # path constraints 
        index = setPathConstraints!(docp, c, index, xu, i)

    end

    # path constraints at final time
    index = setPathConstraints!(docp, c, index, xu, docp.dim_NLP_steps)

    # boundary conditions and variable constraints
    index = setPointConstraints!(docp, c, index, xu)

    # needed even for inplace version, AD error otherwise
    # may be because actual return would be index above ?
    return c
end



"""
$(TYPEDSIGNATURES)

Set the path constraints at given time step
"""
function setPathConstraints!(docp::DOCP, c, index::Int, xu, i::Int)

    ocp = docp.ocp

    ti = get_time_at_time_step(xu, docp, i)
    xi, ui = get_variables_at_time_step(xu, docp, i)
    v = get_optim_variable(xu, docp)

    # pure control constraints
    if docp.dim_u_cons > 0
        c[index:(index + docp.dim_u_cons - 1)] = docp.control_constraints[2](ti, ui, v)
        index = index + docp.dim_u_cons
    end

    # pure state constraints
    if docp.dim_x_cons > 0
        c[index:(index + docp.dim_x_cons - 1)] = docp.state_constraints[2](ti, xi, v)
        index = index + docp.dim_x_cons
    end

    # mixed state / control constraints
    if docp.dim_mixed_cons > 0
        c[index:(index + docp.dim_mixed_cons - 1)] = docp.mixed_constraints[2](ti, xi, ui, v)
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
function setPointConstraints!(docp::DOCP, c, index::Int, xu)

    ocp = docp.ocp

    x0, u0, xl0 = get_variables_at_time_step(xu, docp, 0)
    xf, = get_variables_at_time_step(xu, docp, docp.dim_NLP_steps)
    v = get_optim_variable(xu, docp)

    # boundary constraints
    if docp.dim_boundary_cons > 0
        c[index:(index + docp.dim_boundary_cons - 1)] = docp.boundary_constraints[2](x0, xf, v)
        index = index + docp.dim_boundary_cons
    end

    # variable constraints
    if docp.dim_v_cons > 0
        c[index:(index + docp.dim_v_cons - 1)] = docp.variable_constraints[2](v)
        index = index + docp.dim_v_cons
    end

    # null initial condition for lagrangian cost state
    if docp.has_lagrange
        c[index] = xl0
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
    for i = 0:(docp.dim_NLP_steps)
        ti = get_time_at_time_step(NLP_X, docp, i)
        set_variables_at_time_step!(NLP_X, init.state_init(ti), init.control_init(ti), docp, i)
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
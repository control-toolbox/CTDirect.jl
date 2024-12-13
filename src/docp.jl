# Discretized Optimal Control Problem DOCP

# generic discretization
abstract type Discretization end

"""
$(TYPEDSIGNATURES)

Struct for discretized optimal control problem DOCP

Contains:
- a copy of the original OCP
- data required to link the OCP with the discretized DOCP
"""
struct DOCP{T <: Discretization, O <: CTModels.Model}

    ## OCP
    ocp::O

    # constraints and their bounds
    path_constraints::Any
    boundary_constraints::Any
    variable_constraints::Any
    control_box::Any
    state_box::Any
    variable_box::Any

    # flags
    is_free_initial_time::Bool
    is_free_final_time::Bool
    is_lagrange::Bool
    is_mayer::Bool
    is_maximization::Bool

    # dimensions
    dim_x_box::Int
    dim_u_box::Int
    dim_v_box::Int
    dim_path_cons::Int
    dim_variable_cons::Int
    dim_boundary_cons::Int

    ## NLP  
    dim_NLP_x::Int  # possible lagrange cost
    dim_NLP_u::Int
    dim_NLP_v::Int
    dim_OCP_x::Int  # original OCP state
    dim_NLP_steps::Int
    NLP_normalized_time_grid::Vector{Float64}
    NLP_time_grid::Vector{Float64}
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
    function DOCP(ocp::CTModels.Model; grid_size=__grid_size(), time_grid=__time_grid(), disc_method=__disc_method())

        # time grid
        if time_grid == nothing
            NLP_normalized_time_grid = convert(Vector{Float64}, collect(LinRange(0, 1, grid_size + 1)))
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

        # additional flags
        is_free_initial_time = CTModels.has_free_initial_time(ocp)
        is_free_final_time = CTModels.has_free_final_time(ocp)
        is_lagrange = CTModels.has_lagrange_cost(ocp)
        is_mayer = CTModels.has_mayer_cost(ocp)
        is_maximization = CTModels.criterion(ocp) == :max

        # dimensions
        if is_lagrange
            dim_NLP_x = CTModels.state_dimension(ocp) + 1
        else
            dim_NLP_x = CTModels.state_dimension(ocp)
        end
        dim_NLP_u = CTModels.control_dimension(ocp)
        dim_NLP_v = CTModels.variable_dimension(ocp)
        dim_OCP_x = CTModels.state_dimension(ocp)

        # times
        if is_free_initial_time || is_free_final_time
            # time grid will be recomputed at each NLP iteration
            NLP_time_grid = Vector{Float64}(undef, dim_NLP_steps+1)
        else
            # compute time grid once for all
            t0 = CTModels.initial_time(ocp)
            tf = CTModels.final_time(ocp)
            NLP_time_grid = @. t0 + (NLP_normalized_time_grid * (tf - t0))
        end

        # NLP constraints +++ could be done in Model ?
        # parse NLP constraints (and initialize dimensions)
        path_constraints,
        variable_constraints,
        boundary_constraints,
        state_box,
        control_box,
        variable_box = CTModels.constraints(ocp)

        # get dimensions
        dim_x_box = CTModels.dim_state_cons_box(ocp)
        dim_u_box = CTModels.dim_control_cons_box(ocp)
        dim_v_box = CTModels.dim_variable_cons_box(ocp)
        dim_path_cons = CTModels.dim_path_cons_nl(ocp)
        dim_boundary_cons = CTModels.dim_boundary_cons_nl(ocp)
        dim_variable_cons = CTModels.dim_variable_cons_nl(ocp)

        # parameter: discretization method
        if disc_method == :trapeze
            discretization, dim_NLP_variables, dim_NLP_constraints = CTDirect.Trapeze(dim_NLP_steps, dim_NLP_x, dim_NLP_u, dim_NLP_v, dim_path_cons, dim_boundary_cons, dim_variable_cons)
        #=elseif disc_method == :midpoint
            discretization, dim_NLP_variables, dim_NLP_constraints = CTDirect.Midpoint(dim_NLP_steps, dim_NLP_x, dim_NLP_u, dim_NLP_v, dim_u_cons, dim_x_cons, dim_xu_cons, dim_boundary_cons, dim_variable_cons)
        elseif disc_method == :gauss_legendre_1
                discretization, dim_NLP_variables, dim_NLP_constraints = CTDirect.Gauss_Legendre_1(dim_NLP_steps, dim_NLP_x, dim_NLP_u, dim_NLP_v, dim_u_cons, dim_x_cons, dim_xu_cons, dim_boundary_cons, dim_variable_cons)
        elseif disc_method == :gauss_legendre_2
                discretization, dim_NLP_variables, dim_NLP_constraints = CTDirect.Gauss_Legendre_2(dim_NLP_steps, dim_NLP_x, dim_NLP_u, dim_NLP_v, dim_u_cons, dim_x_cons, dim_xu_cons, dim_boundary_cons, dim_variable_cons)
        elseif disc_method == :gauss_legendre_3
                discretization, dim_NLP_variables, dim_NLP_constraints = CTDirect.Gauss_Legendre_3(dim_NLP_steps, dim_NLP_x, dim_NLP_u, dim_NLP_v, dim_u_cons, dim_x_cons, dim_xu_cons, dim_boundary_cons, dim_variable_cons)=#                                 
        else           
            error("Unknown discretization method: ", disc_method, "\nValid options are disc_method={:trapeze, :midpoint, :gauss_legendre_1, :gauss_legendre_2, :gauss_legendre_3}\n", typeof(disc_method))
        end

        # add initial condition for lagrange state
        if is_lagrange
            dim_NLP_constraints += 1
        end

        # call constructor with const fields
        docp = new{typeof(discretization), typeof(ocp)}(
            ocp,
            path_constraints,
            boundary_constraints,
            variable_constraints,
            control_box,
            state_box,
            variable_box,
            is_free_initial_time,
            is_free_final_time,
            is_lagrange,
            is_mayer,
            is_maximization,
            dim_x_box,
            dim_u_box,
            dim_v_box,
            dim_path_cons,
            dim_variable_cons,
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
    for i = 1:docp.dim_NLP_steps+1
        if i <= docp.dim_NLP_steps
            # skip (ie leave 0) for state / stage equations 
            index = index + docp.discretization._state_stage_eqs_block
        end
        # path constraints
        index = setPathBounds!(docp, index, lb, ub)
    end

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
    for i = 1:N+1
        set_state_at_time_step!(var_l, x_lb, docp, i)
        set_state_at_time_step!(var_u, x_ub, docp, i)
        set_control_at_time_step!(var_l, u_lb, docp, i)
        set_control_at_time_step!(var_u, u_ub, docp, i)
    end

    # variable box
    if docp.dim_NLP_v > 0
        v_lb, v_ub = build_bounds(docp.dim_NLP_v, docp.dim_v_box, docp.variable_box)
        set_optim_variable!(var_l, v_lb, docp)
        set_optim_variable!(var_u, v_ub, docp)
    end

    return var_l, var_u
end

# Q. should we put objective and constraints *in* DOCP ?

"""
$(TYPEDSIGNATURES)

Compute the objective for the DOCP problem.
"""
function DOCP_objective(xu, docp::DOCP)

    obj = similar(xu, 1)
    N = docp.dim_NLP_steps

    # optimization variables
    v = get_OCP_variable(xu, docp)

    # final state is always needed since lagrange cost is there
    xf = get_OCP_state_at_time_step(xu, docp, N+1)

    # mayer cost
    if docp.is_mayer
        x0 = get_OCP_state_at_time_step(xu, docp, 1)
        docp.ocp.mayer(obj, x0, xf, v)
    end

    # lagrange cost
    if docp.is_lagrange
        if docp.is_mayer # NB can this actually happen in OCP (cf bolza) ?
            obj[1] = obj[1] + get_lagrange_state_at_time_step(xu, docp, docp.dim_NLP_steps+1)
        else
            obj[1] = get_lagrange_state_at_time_step(xu, docp, docp.dim_NLP_steps+1)
        end
    end

    # maximization problem
    if docp.is_maximization
        obj[1] = -obj[1]
    end
    
    return obj[1]
end


"""
$(TYPEDSIGNATURES)

Compute the constraints C for the DOCP problem (modeled as LB <= C(X) <= UB).
"""
function DOCP_constraints!(c, xu, docp::DOCP)

    # initialization
    time_grid = get_time_grid(xu, docp)
    v = get_OCP_variable(xu, docp)
    work = setWorkArray(docp, xu, time_grid, v)

    # main loop on time steps
    for i = 1:docp.dim_NLP_steps + 1
        setStepConstraints!(docp, c, xu, v, time_grid, i, work)
    end

    # point constraints (NB. view on c block could be used with offset here)
    setPointConstraints!(docp, c, xu, v)

    # NB. the function *needs* to return c for AD...
    return c
end


"""
$(TYPEDSIGNATURES)

Set path constraints at given time step
"""
function setPathConstraints!(docp, c, ti, xi, ui, v, offset)

    if docp.dim_path_cons > 0
        docp.path_constraints[2]((@view c[offset+1:offset+docp.dim_path_cons]), ti, xi, ui, v)
        offset += docp.dim_path_cons
    end

end

"""
$(TYPEDSIGNATURES)

Set bounds for the path constraints at given time step
"""
function setPathBounds!(docp::DOCP, index::Int, lb, ub)

    if docp.dim_path_cons > 0
        lb[index:(index + docp.dim_path_cons - 1)] = docp.path_constraints[1]
        ub[index:(index + docp.dim_path_cons - 1)] = docp.path_constraints[3]
        index = index + docp.dim_path_cons
    end
    return index

end

"""
$(TYPEDSIGNATURES)

Set the boundary and variable constraints
"""
function setPointConstraints!(docp::DOCP, c, xu, v)

    # offset: [state eq, stage eq, path constraints]_1..N and final path constraints
    offset = docp.dim_NLP_steps * (docp.discretization._state_stage_eqs_block + docp.discretization._step_pathcons_block) + docp.discretization._step_pathcons_block

    # variables
    x0 = get_OCP_state_at_time_step(xu, docp, 1)
    xf = get_OCP_state_at_time_step(xu, docp, docp.dim_NLP_steps+1)

    # boundary constraints
    if docp.dim_boundary_cons > 0
        docp.boundary_constraints[2]((@view c[offset+1:offset+docp.dim_boundary_cons]),x0, xf, v)
    end

    # variable constraints
    if docp.dim_variable_cons > 0
        docp.variable_constraints[2]((@view c[offset+docp.dim_boundary_cons+1:offset+docp.dim_boundary_cons+docp.dim_variable_cons]), v)
    end

    # null initial condition for lagrangian cost state
    if docp.is_lagrange
        c[offset+docp.dim_boundary_cons+docp.dim_variable_cons+1] = get_lagrange_state_at_time_step(xu, docp, 1)
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
    if docp.dim_variable_cons > 0
        lb[index:(index + docp.dim_variable_cons - 1)] = docp.variable_constraints[1]
        ub[index:(index + docp.dim_variable_cons - 1)] = docp.variable_constraints[3]
        index = index + docp.dim_variable_cons
    end

    # null initial condition for lagrangian cost state
    if docp.is_lagrange
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
function DOCP_initial_guess(docp::DOCP, init::CTBase.OptimalControlInit = CTBase.OptimalControlInit())

    # default initialization (internal variables such as lagrange cost, k_i for RK schemes) will keep these default values 
    NLP_X = 0.1 * ones(docp.dim_NLP_variables)

    # set variables if provided (needed first in case of free times !)
    if !isnothing(init.variable_init)
        set_optim_variable!(NLP_X, init.variable_init, docp)
    end

    # set state / control variables if provided (final control case handled by setter)
    time_grid = get_time_grid(NLP_X, docp)
    for i = 1:docp.dim_NLP_steps + 1
        ti = time_grid[i]
        set_state_at_time_step!(NLP_X, init.state_init(ti), docp, i)
        set_control_at_time_step!(NLP_X, init.control_init(ti), docp, i)
    end

    return NLP_X
end

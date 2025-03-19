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
    
    # discretization scheme
    discretization::T

    ## OCP
    ocp::O # parametric instead of just qualifying reduces allocations (but not time). Specialization ?

    # flags
    has_free_initial_time::Bool
    has_free_final_time::Bool
    has_lagrange::Bool
    has_mayer::Bool
    is_maximization::Bool

    # constraints dimensions
    dim_path_cons::Int
    dim_boundary_cons::Int

    ## NLP
    dim_NLP_x::Int  # possible lagrange cost
    dim_NLP_u::Int
    dim_NLP_v::Int
    dim_OCP_x::Int  # original OCP state
    dim_NLP_steps::Int
    NLP_normalized_time_grid::Vector{Float64}
    NLP_fixed_time_grid::Vector{Float64}
    dim_NLP_variables::Int
    dim_NLP_constraints::Int

    # lower and upper bounds for variables and constraints
    var_l::Vector{Float64}
    var_u::Vector{Float64}
    con_l::Vector{Float64}
    con_u::Vector{Float64}

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
        has_free_initial_time = CTModels.has_free_initial_time(ocp)
        has_free_final_time = CTModels.has_free_final_time(ocp)
        has_lagrange = CTModels.has_lagrange_cost(ocp)
        has_mayer = CTModels.has_mayer_cost(ocp)
        is_maximization = CTModels.criterion(ocp) == :max

        # dimensions
        if has_lagrange
            dim_NLP_x = CTModels.state_dimension(ocp) + 1
        else
            dim_NLP_x = CTModels.state_dimension(ocp)
        end
        dim_NLP_u = CTModels.control_dimension(ocp)
        dim_NLP_v = CTModels.variable_dimension(ocp)
        dim_OCP_x = CTModels.state_dimension(ocp)

        # times
        if !has_free_initial_time && !has_free_final_time
            # compute time grid once for all
            t0 = CTModels.initial_time(ocp)
            tf = CTModels.final_time(ocp)
            NLP_fixed_time_grid = @. t0 + (NLP_normalized_time_grid * (tf - t0))
        else
            # time grid will be recomputed at each NLP iteration
            NLP_fixed_time_grid = Vector{Float64}(undef, dim_NLP_steps+1)
        end

        # NLP constraints dimensions
        dim_path_cons = CTModels.dim_path_constraints_nl(ocp)
        dim_boundary_cons = CTModels.dim_boundary_constraints_nl(ocp)

        # parameter: discretization method
        if disc_method == :trapeze
            discretization, dim_NLP_variables, dim_NLP_constraints = CTDirect.Trapeze(dim_NLP_steps, dim_NLP_x, dim_NLP_u, dim_NLP_v, dim_path_cons, dim_boundary_cons)
        elseif disc_method == :midpoint
            discretization, dim_NLP_variables, dim_NLP_constraints = CTDirect.Midpoint(dim_NLP_steps, dim_NLP_x, dim_NLP_u, dim_NLP_v, dim_path_cons, dim_boundary_cons)
        elseif disc_method == :euler
            discretization, dim_NLP_variables, dim_NLP_constraints = CTDirect.Euler(dim_NLP_steps, dim_NLP_x, dim_NLP_u, dim_NLP_v, dim_path_cons, dim_boundary_cons)
        elseif disc_method == :euler_implicit
            discretization, dim_NLP_variables, dim_NLP_constraints = CTDirect.Euler(dim_NLP_steps, dim_NLP_x, dim_NLP_u, dim_NLP_v, dim_path_cons, dim_boundary_cons; explicit=false)
        elseif disc_method == :gauss_legendre_1
                discretization, dim_NLP_variables, dim_NLP_constraints = CTDirect.Gauss_Legendre_1(dim_NLP_steps, dim_NLP_x, dim_NLP_u, dim_NLP_v, dim_path_cons, dim_boundary_cons)
        elseif disc_method == :gauss_legendre_2
                discretization, dim_NLP_variables, dim_NLP_constraints = CTDirect.Gauss_Legendre_2(dim_NLP_steps, dim_NLP_x, dim_NLP_u, dim_NLP_v, dim_path_cons, dim_boundary_cons)
        elseif disc_method == :gauss_legendre_3
                discretization, dim_NLP_variables, dim_NLP_constraints = CTDirect.Gauss_Legendre_3(dim_NLP_steps, dim_NLP_x, dim_NLP_u, dim_NLP_v, dim_path_cons, dim_boundary_cons)
        else           
            error("Unknown discretization method: ", disc_method, "\nValid options are disc_method={:trapeze, :midpoint, :gauss_legendre_1, :gauss_legendre_2, :gauss_legendre_3}\n", typeof(disc_method))
        end

        # add initial condition for lagrange state
        if has_lagrange
            dim_NLP_constraints += 1
        end

        # call constructor with const fields
        docp = new{typeof(discretization), typeof(ocp)}(
            discretization,        
            ocp,
            has_free_initial_time,
            has_free_final_time,
            has_lagrange,
            has_mayer,
            is_maximization,
            dim_path_cons,
            dim_boundary_cons,
            dim_NLP_x,
            dim_NLP_u,
            dim_NLP_v,
            dim_OCP_x,
            dim_NLP_steps,
            NLP_normalized_time_grid,
            NLP_fixed_time_grid,
            dim_NLP_variables,
            dim_NLP_constraints,
            -Inf * ones(dim_NLP_variables),
            Inf * ones(dim_NLP_variables),
            zeros(dim_NLP_constraints),
            zeros(dim_NLP_constraints)
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

    offset = 0
    for i = 1:docp.dim_NLP_steps+1
        if i <= docp.dim_NLP_steps
            # skip (ie leave 0) for state / stage equations 
            offset = offset + docp.discretization._state_stage_eqs_block
        end
        # path constraints
        if docp.dim_path_cons > 0
            lb[offset+1:offset+docp.dim_path_cons] = CTModels.path_constraints_nl(docp.ocp)[1]
            ub[offset+1:offset+docp.dim_path_cons] = CTModels.path_constraints_nl(docp.ocp)[3]
            offset = offset + docp.dim_path_cons
        end
    end

    # boundary constraints
    if docp.dim_boundary_cons > 0
        lb[offset+1:offset+docp.dim_boundary_cons] = CTModels.boundary_constraints_nl(docp.ocp)[1]
        ub[offset+1:offset+docp.dim_boundary_cons] = CTModels.boundary_constraints_nl(docp.ocp)[3]
        offset = offset + docp.dim_boundary_cons
    end

    # null initial condition for lagrangian cost state
    if docp.has_lagrange
        lb[offset+1] = 0.0
        ub[offset+1] = 0.0
        offset = offset + 1
    end

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

    # build full ordered sets of bounds
    x_lb, x_ub = build_bounds(docp.dim_OCP_x, CTModels.dim_state_constraints_box(ocp), CTModels.state_constraints_box(docp.ocp))
    u_lb, u_ub = build_bounds(docp.dim_NLP_u, CTModels.dim_control_constraints_box(ocp), CTModels.control_constraints_box(docp.ocp))

    # set state / control box along time steps
    for i = 1:N+1
        set_state_at_time_step!(var_l, x_lb, docp, i)
        set_state_at_time_step!(var_u, x_ub, docp, i)
        set_control_at_time_step!(var_l, u_lb, docp, i)
        set_control_at_time_step!(var_u, u_ub, docp, i)
    end

    # variable box
    if docp.dim_NLP_v > 0
        v_lb, v_ub = build_bounds(docp.dim_NLP_v, CTModels.dim_variable_constraints_box(ocp), CTModels.variable_constraints_box(docp.ocp))
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

    # mayer cost
    if docp.has_mayer
        v = get_OCP_variable(xu, docp)
        x0 = get_OCP_state_at_time_step(xu, docp, 1)
        xf = get_OCP_state_at_time_step(xu, docp, docp.dim_NLP_steps+1)
        obj_mayer = CTModels.mayer(docp.ocp)(x0, xf, v)
    else
        obj_mayer = 0.0
    end

    # lagrange cost
    if docp.has_lagrange
        obj_lagrange = get_lagrange_state_at_time_step(xu, docp, docp.dim_NLP_steps+1)
    else
        obj_lagrange = 0.0
    end

    # total cost
    obj = obj_mayer + obj_lagrange

    # maximization problem
    if docp.is_maximization
        obj = -obj
    end

    return obj
end


"""
$(TYPEDSIGNATURES)

Compute the constraints C for the DOCP problem (modeled as LB <= C(X) <= UB).
"""
function DOCP_constraints!(c, xu, docp::DOCP)

    # initialization
    if docp.has_free_initial_time || docp.has_free_final_time
        time_grid = get_time_grid(xu, docp)
    else
        time_grid = docp.NLP_fixed_time_grid
    end
    v = get_OCP_variable(xu, docp)
    work = setWorkArray(docp, xu, time_grid, v)

    # main loop on time steps
    for i = 1:docp.dim_NLP_steps + 1
        setStepConstraints!(docp, c, xu, v, time_grid, i, work)
    end

    # point constraints
    setPointConstraints!(docp, c, xu, v)

    # NB. the function *needs* to return c for AD...
    return c
end


"""
$(TYPEDSIGNATURES)

Set the boundary and variable constraints
"""
function setPointConstraints!(docp::DOCP, c, xu, v)

    offset = docp.dim_NLP_constraints - docp.dim_boundary_cons
    docp.has_lagrange && (offset = offset - 1)

    # variables
    x0 = get_OCP_state_at_time_step(xu, docp, 1)
    xf = get_OCP_state_at_time_step(xu, docp, docp.dim_NLP_steps+1)

    # boundary constraints
    if docp.dim_boundary_cons > 0
        CTModels.boundary_constraints_nl(docp.ocp)[2]((@view c[offset+1:offset+docp.dim_boundary_cons]),x0, xf, v)
	offset = offset + docp.dim_boundary_cons
    end

    # null initial condition for lagrangian cost state
    if docp.has_lagrange
        c[offset+1] = get_lagrange_state_at_time_step(xu, docp, 1)
    end
end


"""
$(TYPEDSIGNATURES)

Build initial guess for discretized problem
"""
function DOCP_initial_guess(docp::DOCP, init::CTModels.Init = CTModels.Init())

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


"""
$(TYPEDSIGNATURES)

Return time grid for variable time problems (times are then dependent on NLP variables)
"""
function get_time_grid(xu, docp::DOCP)

    grid = similar(xu, docp.dim_NLP_steps+1)
    ocp = docp.ocp
    
    if docp.has_free_initial_time
        v = get_OCP_variable(xu, docp)
        t0 = CTModels.initial_time(ocp, v)
    else
        t0 = CTModels.initial_time(ocp)
    end
    if docp.has_free_final_time
        v = get_OCP_variable(xu, docp)
        tf = CTModels.final_time(ocp, v)
    else
        tf = CTModels.final_time(ocp)
    end

    @. grid = t0 + docp.NLP_normalized_time_grid * (tf - t0)
    return grid

end


"""
$(TYPEDSIGNATURES)

Build full, ordered sets of bounds for state, control or optimization variables
"""
function build_bounds(dim_var, dim_box, box_triplet)

    x_lb = -Inf * ones(dim_var)
    x_ub = Inf * ones(dim_var)
    for j = 1:(dim_box)
        indice = box_triplet[2][j]
        x_lb[indice] = box_triplet[1][j]
        x_ub[indice] = box_triplet[3][j]
    end

    return x_lb, x_ub
end

# Discretized Optimal Control Problem DOCP

# generic discretization
abstract type Discretization end
abstract type ArgsAtStep end
abstract type ScalVect end
struct ScalVariable <: ScalVect end
struct VectVariable <: ScalVect end

# use OCP model subtype later ?
abstract type TimeGrid end
struct FixedTimeGrid <: TimeGrid end
struct FreeTimeGrid <: TimeGrid end

"""
$(TYPEDSIGNATURES)

Struct for discretized optimal control problem DOCP

Contains:
- a copy of the original OCP
- data required to link the OCP with the discretized DOCP
"""
struct DOCP{T <: Discretization, X <: ScalVect, U <: ScalVect, V <: ScalVect, G <: TimeGrid}

    ## OCP
    ocp::OptimalControlModel # remove at some point ?

    # constraints and their bounds
    control_constraints::Any
    state_constraints::Any
    mixed_constraints::Any
    boundary_constraints::Any
    variable_constraints::Any
    control_box::Any
    state_box::Any
    variable_box::Any

    # flags
    time_grid_type:: G
    is_free_initial_time::Bool
    is_free_final_time::Bool
    is_lagrange::Bool
    is_mayer::Bool
    is_variable::Bool
    is_maximization::Bool

    # initial / final time
    fixed_initial_time::Float64
    fixed_final_time::Float64
    index_initial_time::Index
    index_final_time::Index

    # dimensions
    dim_x_box::Int
    dim_u_box::Int
    dim_v_box::Int
    dim_x_cons::Int
    dim_u_cons::Int
    dim_v_cons::Int
    dim_xu_cons::Int
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

    # discretization scheme
    discretization::T

    # scalar / vector aux function
    _type_x::X
    _type_u::U  
    _type_v::V

    # constructor
    function DOCP(ocp::OptimalControlModel; grid_size=__grid_size(), time_grid=__time_grid(), disc_method=__disc_method(), control_type=__control_type())

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
        is_free_initial_time = has_free_initial_time(ocp)
        is_free_final_time = has_free_final_time(ocp)
        is_lagrange = has_lagrange_cost(ocp)
        is_mayer = has_mayer_cost(ocp)
        is_variable = is_variable_dependent(ocp)
        is_maximization = is_max(ocp)

        # dimensions
        if is_lagrange
            dim_NLP_x = ocp.state_dimension + 1
        else
            dim_NLP_x = ocp.state_dimension
        end
        dim_NLP_u = ocp.control_dimension
        if is_variable
            dim_NLP_v = ocp.variable_dimension
        else
            dim_NLP_v = 0 # dim in ocp would be Nothing
        end
        dim_OCP_x = ocp.state_dimension

        # times
        # use 2 different variables for value / index (type stability)
        if is_free_initial_time
            index_initial_time = ocp.initial_time
            fixed_initial_time = 0. # unused
        else
            fixed_initial_time = ocp.initial_time
            index_initial_time = Index(1) # unused
        end
        if is_free_final_time
            index_final_time = ocp.final_time
            fixed_final_time = 0. # unused
        else
            fixed_final_time = ocp.final_time
            index_final_time = Index(1) # unused
        end

        if !is_free_initial_time && !is_free_final_time
            # compute time grid once for all
            NLP_fixed_time_grid = @. fixed_initial_time + (NLP_normalized_time_grid * (fixed_final_time - fixed_initial_time))
            time_grid_type = FixedTimeGrid()
        else
            # unused
            NLP_fixed_time_grid = Vector{Float64}(undef, dim_NLP_steps+1)
            time_grid_type = FreeTimeGrid()
        end

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

        # get dimensions
        dim_x_box = dim_state_range(ocp)
        dim_u_box = dim_control_range(ocp)
        dim_v_box = dim_variable_range(ocp)
        dim_x_cons = dim_state_constraints(ocp)
        dim_u_cons = dim_control_constraints(ocp)
        dim_v_cons = dim_variable_constraints(ocp)
        dim_xu_cons = dim_mixed_constraints(ocp)
        dim_boundary_cons = dim_boundary_constraints(ocp)

        # parameter: discretization method
        if disc_method == :trapeze
            discretization, dim_NLP_variables, dim_NLP_constraints = CTDirect.Trapeze(dim_NLP_steps, dim_NLP_x, dim_NLP_u, dim_NLP_v, dim_u_cons, dim_x_cons, dim_xu_cons, dim_boundary_cons, dim_v_cons)
        elseif disc_method == :trapeze_stage
                discretization, dim_NLP_variables, dim_NLP_constraints = CTDirect.Trapeze_stage(dim_NLP_steps, dim_NLP_x, dim_NLP_u, dim_NLP_v, dim_u_cons, dim_x_cons, dim_xu_cons, dim_boundary_cons, dim_v_cons)            
        elseif disc_method == :midpoint
            discretization, dim_NLP_variables, dim_NLP_constraints = CTDirect.Midpoint(dim_NLP_steps, dim_NLP_x, dim_NLP_u, dim_NLP_v, dim_u_cons, dim_x_cons, dim_xu_cons, dim_boundary_cons, dim_v_cons)
        elseif disc_method == :midpoint_stage
            discretization, dim_NLP_variables, dim_NLP_constraints = CTDirect.Midpoint_stage(dim_NLP_steps, dim_NLP_x, dim_NLP_u, dim_NLP_v, dim_u_cons, dim_x_cons, dim_xu_cons, dim_boundary_cons, dim_v_cons)
        elseif disc_method == :midpoint_nowork
            discretization, dim_NLP_variables, dim_NLP_constraints = CTDirect.Midpoint_nowork(dim_NLP_steps, dim_NLP_x, dim_NLP_u, dim_NLP_v, dim_u_cons, dim_x_cons, dim_xu_cons, dim_boundary_cons, dim_v_cons)            
        elseif disc_method == :euler
            discretization, dim_NLP_variables, dim_NLP_constraints = CTDirect.Euler(dim_NLP_steps, dim_NLP_x, dim_NLP_u, dim_NLP_v, dim_u_cons, dim_x_cons, dim_xu_cons, dim_boundary_cons, dim_v_cons)
        elseif disc_method == :euler_implicit
            discretization, dim_NLP_variables, dim_NLP_constraints = CTDirect.Euler(dim_NLP_steps, dim_NLP_x, dim_NLP_u, dim_NLP_v, dim_u_cons, dim_x_cons, dim_xu_cons, dim_boundary_cons, dim_v_cons; explicit=false)
        elseif disc_method == :gauss_legendre_1
                discretization, dim_NLP_variables, dim_NLP_constraints = CTDirect.Gauss_Legendre_1(dim_NLP_steps, dim_NLP_x, dim_NLP_u, dim_NLP_v, dim_u_cons, dim_x_cons, dim_xu_cons, dim_boundary_cons, dim_v_cons)
        elseif disc_method == :gauss_legendre_2
                discretization, dim_NLP_variables, dim_NLP_constraints = CTDirect.Gauss_Legendre_2(dim_NLP_steps, dim_NLP_x, dim_NLP_u, dim_NLP_v, dim_u_cons, dim_x_cons, dim_xu_cons, dim_boundary_cons, dim_v_cons, control_type)
        elseif disc_method == :gauss_legendre_3
                discretization, dim_NLP_variables, dim_NLP_constraints = CTDirect.Gauss_Legendre_3(dim_NLP_steps, dim_NLP_x, dim_NLP_u, dim_NLP_v, dim_u_cons, dim_x_cons, dim_xu_cons, dim_boundary_cons, dim_v_cons, control_type)
        elseif disc_method == :euler_irk
                discretization, dim_NLP_variables, dim_NLP_constraints = CTDirect.Euler_explicit(dim_NLP_steps, dim_NLP_x, dim_NLP_u, dim_NLP_v, dim_u_cons, dim_x_cons, dim_xu_cons, dim_boundary_cons, dim_v_cons)
        elseif disc_method == :euler_implicit_irk
                discretization, dim_NLP_variables, dim_NLP_constraints = CTDirect.Euler_implicit(dim_NLP_steps, dim_NLP_x, dim_NLP_u, dim_NLP_v, dim_u_cons, dim_x_cons, dim_xu_cons, dim_boundary_cons, dim_v_cons)
        elseif disc_method == :trapeze_irk
                discretization, dim_NLP_variables, dim_NLP_constraints = CTDirect.Trapeze_implicit(dim_NLP_steps, dim_NLP_x, dim_NLP_u, dim_NLP_v, dim_u_cons, dim_x_cons, dim_xu_cons, dim_boundary_cons, dim_v_cons, control_type)           
        else           
            error("Unknown discretization method: ", disc_method, "\nValid options are disc_method={:euler, :euler_implicit, :trapeze, :trapeze_irk, :midpoint, :gauss_legendre_2, :gauss_legendre_3, (:euler_irk, :euler_implicit_irk, :gauss_legendre_1, :trapeze_irk)}\n", typeof(disc_method))
        end

        # add initial condition for lagrange state
        if is_lagrange
            dim_NLP_constraints += 1
        end

        # parameters: scalar / vector for state, control, variable 
        dim_OCP_x == 1 ? _type_x = ScalVariable() : _type_x = VectVariable()
        dim_NLP_u == 1 ? _type_u = ScalVariable() : _type_u = VectVariable()
        dim_NLP_v == 1 ? _type_v = ScalVariable() : _type_v = VectVariable()

        # call constructor with const fields
        docp = new{typeof(discretization), typeof(_type_x), typeof(_type_u), typeof(_type_v), typeof(time_grid_type)}(
            ocp,
            control_constraints,
            state_constraints,
            mixed_constraints,
            boundary_constraints,
            variable_constraints,
            control_box,
            state_box,
            variable_box,
            time_grid_type,
            is_free_initial_time,
            is_free_final_time,
            is_lagrange,
            is_mayer,
            is_variable,
            is_maximization,
            fixed_initial_time,
            fixed_final_time,
            index_initial_time,
            index_final_time,
            dim_x_box,
            dim_u_box,
            dim_v_box,
            dim_x_cons,
            dim_u_cons,
            dim_v_cons,
            dim_xu_cons,
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
            zeros(dim_NLP_constraints),
            discretization,
            _type_x,
            _type_u,
            _type_v
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
        #+++ setStepConstraintsBounds in disc/  different for trapeze_stage !
        if i <= docp.dim_NLP_steps
            # skip (ie leave 0) for state / stage equations 
            offset = offset + docp.discretization._state_stage_eqs_block
        end
        # path constraints
        offset = setPathConstraintsBounds!(docp, lb, ub, offset)
    end

    # boundary and variable constraints
    offset = setPointConstraintsBounds!(docp, lb, ub)
    if offset != docp.dim_NLP_constraints
        error("Mismatch for last index in constraints: ", offset, " instead of ", docp.dim_NLP_constraints)
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
    x_lb, x_ub = build_bounds(docp.dim_OCP_x, docp.dim_x_box, docp.state_box)
    u_lb, u_ub = build_bounds(docp.dim_NLP_u, docp.dim_u_box, docp.control_box)
    
    # set state / control box along time steps
    for i = 1:N+1
        set_state_at_time_step!(var_l, x_lb, docp, i)
        set_state_at_time_step!(var_u, x_ub, docp, i)
        set_control_at_time_step!(var_l, u_lb, docp, i)
        set_control_at_time_step!(var_u, u_ub, docp, i)
    end

    # variable box
    if docp.is_variable
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
        if docp.is_mayer # bolza case
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

    # main loop on time steps (using c block view seems similar)
    for i = 1:docp.dim_NLP_steps + 1
        setStepConstraints!(docp, c, xu, v, time_grid, i, work)
        #offset = (i-1)*(docp.discretization._state_stage_eqs_block + docp.discretization._step_pathcons_block)
        #setStepConstraints!(docp, (@view c[offset+1:offset+docp.dim_NLP_x+docp.discretization._step_pathcons_block]), xu, v, time_grid, i, work)
    end

    # point constraints
    setPointConstraints!(docp, c, xu, v)

    # NB. the function *needs* to return c for AD...
    return c
end


"""
$(TYPEDSIGNATURES)

Set path constraints at given time step
"""
function setPathConstraints!(docp, c, ti, xi, ui, v, offset)

    # control constraints
    if docp.dim_u_cons > 0
        docp.control_constraints[2]((@view c[offset+1:offset+docp.dim_u_cons]),ti, ui, v)
        offset += docp.dim_u_cons
    end

    # state constraints
    if docp.dim_x_cons > 0 
        docp.state_constraints[2]((@view c[offset+1:offset+docp.dim_x_cons]),ti, xi, v)
        offset += docp.dim_x_cons
    end

    # mixed constraints
    if docp.dim_xu_cons > 0
        docp.mixed_constraints[2]((@view c[offset+1:offset+docp.dim_xu_cons]), ti, xi, ui, v)
        offset += docp.dim_xu_cons
    end
end

"""
$(TYPEDSIGNATURES)

Set bounds for the path constraints at given time step
"""
function setPathConstraintsBounds!(docp::DOCP, lb, ub, offset)

    # pure control constraints
    if docp.dim_u_cons > 0
        lb[offset+1:offset+docp.dim_u_cons] = docp.control_constraints[1]
        ub[offset+1:offset+docp.dim_u_cons] = docp.control_constraints[3]
        offset = offset + docp.dim_u_cons
    end

    # pure state constraints
    if docp.dim_x_cons > 0
        lb[offset+1:offset+docp.dim_x_cons] = docp.state_constraints[1]
        ub[offset+1:offset+docp.dim_x_cons] = docp.state_constraints[3]
        offset = offset + docp.dim_x_cons
    end

    # mixed state / control constraints
    if docp.dim_xu_cons > 0
        lb[offset+1:offset+docp.dim_xu_cons] = docp.mixed_constraints[1]
        ub[offset+1:offset+docp.dim_xu_cons] = docp.mixed_constraints[3]
        offset = offset + docp.dim_xu_cons
    end

    return offset
end

"""
$(TYPEDSIGNATURES)

Set the boundary and variable constraints
"""
function setPointConstraints!(docp::DOCP, c, xu, v)

    offset = docp.dim_NLP_constraints - (docp.dim_boundary_cons + docp.dim_v_cons)
    docp.is_lagrange && (offset = offset - 1)

    # variables
    x0 = get_OCP_state_at_time_step(xu, docp, 1)
    xf = get_OCP_state_at_time_step(xu, docp, docp.dim_NLP_steps+1)

    # boundary constraints
    if docp.dim_boundary_cons > 0
        docp.boundary_constraints[2]((@view c[offset+1:offset+docp.dim_boundary_cons]),x0, xf, v)
        offset = offset + docp.dim_boundary_cons
    end

    # variable constraints
    if docp.dim_v_cons > 0
        docp.variable_constraints[2]((@view c[offset+1:offset+docp.dim_v_cons]), v)
        offset = offset + docp.dim_v_cons
    end

    # null initial condition for lagrangian cost state
    if docp.is_lagrange
        c[offset+1] = get_lagrange_state_at_time_step(xu, docp, 1)
    end
end


"""
$(TYPEDSIGNATURES)

Set bounds for the boundary and variable constraints
"""
function setPointConstraintsBounds!(docp::DOCP, lb, ub)
    
    offset = docp.dim_NLP_constraints - (docp.dim_boundary_cons + docp.dim_v_cons)
    docp.is_lagrange && (offset = offset - 1)

    # boundary constraints
    if docp.dim_boundary_cons > 0
        lb[offset+1:offset+docp.dim_boundary_cons] = docp.boundary_constraints[1]
        ub[offset+1:offset+docp.dim_boundary_cons] = docp.boundary_constraints[3]
        offset = offset + docp.dim_boundary_cons
    end

    # variable constraints
    if docp.dim_v_cons > 0
        lb[offset+1:offset+docp.dim_v_cons] = docp.variable_constraints[1]
        ub[offset+1:offset+docp.dim_v_cons] = docp.variable_constraints[3]
        offset = offset + docp.dim_v_cons
    end

    # null initial condition for lagrangian cost state
    if docp.is_lagrange
        lb[offset+1] = 0.0
        ub[offset+1] = 0.0
        offset = offset + 1
    end

    return offset
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

    # set state / control variables if provided (final control case handled by setter)
    time_grid = get_time_grid(NLP_X, docp)
    for i = 1:docp.dim_NLP_steps + 1
        ti = time_grid[i]
        set_state_at_time_step!(NLP_X, init.state_init(ti), docp, i)
        set_control_at_time_step!(NLP_X, init.control_init(ti), docp, i)
    end

    return NLP_X
end


function get_time_grid(xu, docp::DOCP{<:Discretization, <: ScalVect, <: ScalVect, <: ScalVect, FixedTimeGrid})
    return docp.NLP_fixed_time_grid
end
function get_time_grid(xu, docp::DOCP{<:Discretization, <: ScalVect, <: ScalVect, <: ScalVect, FreeTimeGrid})

    if docp.is_free_initial_time
        v = get_OCP_variable(xu, docp)
        t0 = v[docp.index_initial_time]
    else
        t0 = docp.fixed_initial_time
    end

    if docp.is_free_final_time
        v = get_OCP_variable(xu, docp)
        tf = v[docp.index_final_time]
    else
        tf = docp.fixed_final_time
    end

    grid = similar(xu, docp.dim_NLP_steps+1)
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

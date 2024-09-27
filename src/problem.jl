# Discretized Optimal Control Problem DOCP
# Notes:
# - for now the path constraints are checked on the time steps ie g(t_i, x_i, u_i). This requires a getter that provide a control value at each time step, which may not coincide with the actual control discretization (eg RK stages). In this case we will use the 'average' control using the same coefficients as the RK method. Later we can add an option for control discretization, step or stage. Taking a control constant per step would solve the question for the path constraints evaluation and reduces the number of variables. Compared to the more standard stage control discretization, the consistency of the costate would need to be checked, as well as the oscillations in the trajectory (which may actually be better). Further options may include CVP (control vector parametrization) on a coarser grid.
# - for the other choice of enforcing path constraints at the time stages, the symmetric question of getting state values at time stages can be solved by reusing the states used in the evaluation of the stage dynamics. This second choice (as in Bocop2) has the drawback of a larger problem size, and does not check the constraints at the points of the actual trajectory computed (including tf).

# generic discretization
abstract type Discretization end
abstract type ArgsAtStep end
abstract type ScalVect end
struct ScalVariable <: ScalVect end
struct VectVariable <: ScalVect end

"""
$(TYPEDSIGNATURES)

Struct for discretized optimal control problem DOCP

Contains:
- a copy of the original OCP
- data required to link the OCP with the discretized DOCP
"""
struct DOCP{T <: Discretization, X <: ScalVect, V <: ScalVect}

    ## OCP
    ocp::OptimalControlModel # remove at some point ?

    # functions
    dynamics_ext!::Function

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
    is_free_initial_time::Bool
    is_free_final_time::Bool
    is_lagrange::Bool
    is_mayer::Bool
    is_variable::Bool
    is_maximization::Bool
    is_inplace::Bool

    # initial / final time
    fixed_initial_time::Float64
    fixed_final_time::Float64
    index_initial_time::Index
    index_final_time::Index

    # dimensions
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

    # scalar / vector aux function
    _vec::Function
    _x::Function
    _u::Function
    _v::Function
    _type_x::X    
    _type_v::V

    # constructor
    function DOCP(ocp::OptimalControlModel; grid_size=__grid_size(), time_grid=__time_grid(), disc_method="trapeze")

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
        is_inplace = is_in_place(ocp)

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
        if is_free_initial_time || is_free_final_time
            # time grid will be recomputed at each NLP iteration
            NLP_time_grid = Vector{Any}(undef, dim_NLP_steps+1)
        else
            # compute time grid once for all 
            NLP_time_grid = @. ocp.initial_time + (NLP_normalized_time_grid * (ocp.final_time - ocp.initial_time))
        end

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
        dim_path_cons = dim_path_constraints(ocp)
        dim_x_cons = dim_state_constraints(ocp)
        dim_u_cons = dim_control_constraints(ocp)
        dim_v_cons = dim_variable_constraints(ocp)
        dim_mixed_cons = dim_mixed_constraints(ocp)
        dim_boundary_cons = dim_boundary_constraints(ocp)

        # encapsulated OCP functions
        _vec(f::AbstractVector) = f
        _vec(f::Number) = [f]


        # x[1:dim_OCP_x] adds allocations even with a view
        # returning just x also allocates, although a bit less
        _x(x) = (dim_OCP_x == 1) ? x[1] : x[1:dim_OCP_x]
        _u(u) = (dim_NLP_u == 1) ? u[1] : u
        _v(v) = (dim_NLP_v == 1) ? v[1] : v
        # NB. defining 2 functions inside the if is not better
        if dim_NLP_v == 1
            _type_v = ScalVariable()
        else
            _type_v = VectVariable()
        end
        if dim_OCP_x == 1
            _type_x = ScalVariable()
        else
            _type_x = VectVariable()
        end
        #= mayer cost
        function mayer(obj, x0, xf, v) 
        
            if dim_OCP_x == 1
                _x0 = x0[1]
                _xf = xf[1]
            else
                _x0 = x0
                _xf = xf
            end
            if dim_NLP_v == 1
                _v = v[1]
            else
                _v = v
            end
            if is_inplace
                ocp.mayer(obj, _x0, _xf, _v)
            else
                obj[1] = ocp.mayer(_x0, _xf, _v)
            end
            return
        end=#

        # extended dynamics with lagrange cost
        # use only one function ?
        if is_inplace
            if is_lagrange
                dynamics_ext! = function (f, t, x, u, v)
                    ocp.dynamics((@view f[1:dim_OCP_x]), t, _x(x), _u(u), _v(v))
                    ocp.lagrange((@view f[dim_NLP_x:dim_NLP_x]), t, _x(x), _u(u), _v(v))
                    return
                end

            else
                dynamics_ext! = (f, t, x, u, v) -> ocp.dynamics(f, t, _x(x), _u(u), _v(v))
            end
        else
            if is_lagrange
                dynamics_ext! = function (f, t, x, u, v)
                    f[1:docp.dim_OCP_x] .= ocp.dynamics(t, _x(x), _u(u), _v(v)) # .= required for scalar case
                    f[docp.dim_NLP_x] = ocp.lagrange(t, _x(x), _u(u), _v(v))
                    return
                end                
            else
                dynamics_ext! = function (f, t, x, u, v)
                    f[:] = ocp.dynamics(t, _x(x), _u(u), _v(v))
                    return
                end
            end
        end


        # discretization
        if disc_method == "midpoint"
            discretization = CTDirect.Midpoint(dim_NLP_x, dim_NLP_u, dim_NLP_steps)
        elseif disc_method == "trapeze"
            discretization = CTDirect.Trapeze(dim_NLP_x, dim_NLP_u)
        else
            error("Unknown discretization method:", disc_method)
        end

        # NLP variables size (state, control, variable, stage)
        dim_stage = discretization.stage
        dim_NLP_variables = (dim_NLP_steps + 1) * dim_NLP_x + (dim_NLP_steps + discretization.additional_controls) * dim_NLP_u + dim_NLP_v + dim_NLP_steps * dim_NLP_x * dim_stage

        # NLP constraints size (dynamics, stage, path, boundary, variable)
        dim_NLP_constraints =
        dim_NLP_steps * (dim_NLP_x + (dim_NLP_x * dim_stage) + dim_path_cons) + dim_path_cons + dim_boundary_cons + dim_v_cons
        if is_lagrange
            # add initial condition for lagrange state
            dim_NLP_constraints += 1
        end

        # call constructor with const fields
        docp = new{typeof(discretization), typeof(_type_x), typeof(_type_v)}(
            ocp,
            dynamics_ext!,
            control_constraints,
            state_constraints,
            mixed_constraints,
            boundary_constraints,
            variable_constraints,
            control_box,
            state_box,
            variable_box,
            is_free_initial_time,
            is_free_final_time,
            is_lagrange,
            is_mayer,
            is_variable,
            is_maximization,
            is_inplace,
            fixed_initial_time,
            fixed_final_time,
            index_initial_time,
            index_final_time,
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
            discretization,
            _vec,
            _x,
            _u,
            _v,
            _type_x,
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
function DOCP_objective_param(xu, docp::DOCP)

    obj = similar(xu, 1)
    N = docp.dim_NLP_steps

    # optimization variables
    v = get_OCP_variable_param(xu, docp)

    # final state is always needed since lagrange cost is there
    xf = get_OCP_state_at_time_step_param(xu, docp, N+1)

    # mayer cost
    if docp.is_mayer
        x0 = get_OCP_state_at_time_step_param(xu, docp, 1)
        if docp.is_inplace
            docp.ocp.mayer(obj, x0, xf, v)
        else
            obj[1] = docp.ocp.mayer(x0, xf, v)
        end
    end

    # lagrange cost
    if docp.is_lagrange
        error("get lagrange state")
        if docp.is_mayer # NB can this actually happen in OCP (cf bolza) ?
            obj[1] = obj[1] + xf[end]
        else
            obj[1] = xf[end]
        end
    end

    # maximization problem
    if docp.is_maximization
        obj[1] = -obj[1]
    end
    
    return obj[1]
end

function DOCP_objective_OCP(xu, docp::DOCP)

    obj = similar(xu, 1)
    N = docp.dim_NLP_steps

    # optimization variables
    v = get_OCP_variable(xu, docp)

    # final state is always needed since lagrange cost is there
    xf = get_OCP_state_at_time_step(xu, docp, N+1)

    # mayer cost
    if docp.is_mayer
        x0 = get_OCP_state_at_time_step(xu, docp, 1)
        if docp.is_inplace
            docp.ocp.mayer(obj, x0, xf, v)
        else
            obj[1] = docp.ocp.mayer(x0, xf, v)
        end
    end

    # lagrange cost
    if docp.is_lagrange
        error("get lagrange state")
        if docp.is_mayer # NB can this actually happen in OCP (cf bolza) ?
            obj[1] = obj[1] + xf[end]
        else
            obj[1] = xf[end]
        end
    end

    # maximization problem
    if docp.is_maximization
        obj[1] = -obj[1]
    end
    
    return obj[1]
end

function DOCP_objective(xu, docp::DOCP)

    obj = similar(xu, 1)
    N = docp.dim_NLP_steps

    # optimization variables
    v = get_optim_variable(xu, docp)

    # final state is always needed since lagrange cost is there
    xf = get_state_at_time_step(xu, docp, N+1)

    # mayer cost
    if docp.is_mayer
        x0 = get_state_at_time_step(xu, docp, 1)
        if docp.is_inplace
            docp.ocp.mayer(obj, docp._x(x0), docp._x(xf), docp._v(v))
        else
            obj[1] = docp.ocp.mayer(docp._x(x0), docp._x(xf), docp._v(v))
        end
    end

    # lagrange cost
    if docp.is_lagrange
        if docp.is_mayer # NB can this actually happen in OCP (cf bolza) ?
            obj[1] = obj[1] + xf[end]
        else
            obj[1] = xf[end]
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
    if docp.is_free_initial_time || docp.is_free_final_time
        get_time_grid!(xu, docp)
    end
    v = get_optim_variable(xu, docp)
    # NB using inplace here did not seem to give better results
    work = setWorkArray(docp, xu, docp.NLP_time_grid, v)

    # main loop on time steps 
    for i = 1:docp.dim_NLP_steps+1
        setConstraintBlock!(docp, c, xu, v, docp.NLP_time_grid, i, work)
    end

    # point constraints
    setPointConstraints!(docp, c, xu, v)

    # NB. the function *needs* to return c for AD...
    return c
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
    x0 = get_state_at_time_step(xu, docp, 1)
    xf = get_state_at_time_step(xu, docp, docp.dim_NLP_steps+1)

    # boundary constraints
    if docp.dim_boundary_cons > 0
        if docp.is_inplace
            docp.boundary_constraints[2]((@view c[offset+1:offset+docp.dim_boundary_cons]), docp._x(x0), docp._x(xf), v)
        else
            c[offset+1:offset+docp.dim_boundary_cons] = docp.boundary_constraints[2](docp._x(x0), docp._x(xf), v)
        end
    end

    # variable constraints
    if docp.dim_v_cons > 0
        if docp.is_inplace
            docp.variable_constraints[2]((@view c[offset+docp.dim_boundary_cons+1:offset+docp.dim_boundary_cons+docp.dim_v_cons]), v)
        else
            c[offset+docp.dim_boundary_cons+1:offset+docp.dim_boundary_cons+docp.dim_v_cons] = docp.variable_constraints[2](v)
        end
    end

    # null initial condition for lagrangian cost state
    if docp.is_lagrange
        c[offset+docp.dim_boundary_cons+docp.dim_v_cons+1] = x0[end]
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
function DOCP_initial_guess(docp::DOCP, init::OptimalControlInit = OptimalControlInit())

    # default initialization (internal variables such as lagrange cost, k_i for RK schemes) will keep these default values 
    NLP_X = 0.1 * ones(docp.dim_NLP_variables)

    # set variables if provided (needed first in case of free times !)
    if !isnothing(init.variable_init)
        set_optim_variable!(NLP_X, init.variable_init, docp)
    end

    # set state / control variables if provided
    get_time_grid!(NLP_X, docp)
    for i = 1:docp.dim_NLP_steps+1
        ti = docp.NLP_time_grid[i]
        set_state_at_time_step!(NLP_X, init.state_init(ti), docp, i)
        set_control_at_time_step!(NLP_X, init.control_init(ti), docp, i)
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

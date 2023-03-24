module CTDirect

# using
#
using CTBase

# nlp modeling and resolution
using NLPModelsIpopt, ADNLPModels

# Other declarations
const nlp_constraints = CTBase.nlp_constraints
const __grid_size_direct = CTBase.__grid_size_direct
const __print_level_ipopt = CTBase.__print_level_ipopt
const __mu_strategy_ipopt = CTBase.__mu_strategy_ipopt
const __display = CTBase.__display
const matrix2vec = CTBase.matrix2vec

# includes
include("utils.jl")
include("problem.jl")
include("solution.jl")
include("solve.jl")

# struct for ocop/nlp info
mutable struct CTDirect_data

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

    # NLP
    dim_NLP_state
    dim_NLP_constraints
    dim_NLP_variables
    dim_NLP_steps
    dynamics_lagrange_to_mayer
    NLP_init

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

# export functions only for user
export solve

end
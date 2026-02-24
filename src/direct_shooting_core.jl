# Direct Shooting discretizer: core functions

# ==========================================================================================
# Build core DOCP structure with discretization information (ADNLP)
# ==========================================================================================
function get_docp(discretizer::DirectShooting, ocp::AbstractModel)
    
    # recover discretization scheme and options
    scheme = Strategies.options(discretizer)[:scheme]
    grid_size = Strategies.options(discretizer)[:grid_size]
    control_steps = Strategies.options(discretizer)[:control_steps]

    # initialize DOCP
    docp = DOCP_DirectShooting(ocp, grid_size, control_steps, scheme)

    # set bounds in DOCP
    variables_bounds!(docp)
    constraints_bounds!(docp)

    return docp
end

# could we reuse the same DOCP as collocation
# pass a keyword for the discretizer, as for scheme ?
mutable struct DOCP_DirectShooting{
    D<:CTDirect.Discretization, O<:CTModels.Model
    }

    # discretization scheme
    discretization::D

    # OCP (for OCP functions calls and some getters eg free times)
    # parametric instead of just qualifying reduces allocations (but not time). Specialization ?
    ocp::O

    # boolean flags
    flags::DOCPFlags

    # dimensions
    dims::DOCPdims

    # time grid
    time::DOCPtime

    # lower and upper bounds for variables and constraints
    bounds::DOCPbounds

    # NLP variables and constraints
    dim_NLP_variables::Int
    dim_NLP_constraints::Int

    # constructor
    function DOCP_DirectShooting(ocp::CTModels.Model, grid_size, control_steps, scheme)

        # boolean flags
        flags = DOCPFlags(ocp)

        # dimensions
        dims = DOCPdims(ocp)

        # time grid
        time = DOCPtime(ocp, grid_size, nothing)

        # discretization method 
        disc_args = [:direct_shooting, dims, time.steps, control_steps]
        
        if scheme == :trapeze
            discretization, dim_NLP_variables, dim_NLP_constraints = 
            CTDirect.Trapeze(disc_args...)

        elseif scheme == :midpoint
            discretization, dim_NLP_variables, dim_NLP_constraints = 
            CTDirect.Midpoint(disc_args...)

        elseif scheme == :euler || scheme == :euler_explicit || scheme == :euler_forward
            discretization, dim_NLP_variables, dim_NLP_constraints = 
            CTDirect.Euler(disc_args...)
        elseif scheme == :euler_implicit || scheme == :euler_backward
            discretization, dim_NLP_variables, dim_NLP_constraints = 
            CTDirect.Euler(disc_args...; explicit=false)

        elseif scheme == :gauss_legendre_2
            discretization, dim_NLP_variables, dim_NLP_constraints = 
            CTDirect.Gauss_Legendre_2(disc_args...)

        elseif scheme == :gauss_legendre_3
            discretization, dim_NLP_variables, dim_NLP_constraints = 
            CTDirect.Gauss_Legendre_3(disc_args...)

        else
            error(
                "Unknown discretization method: ",
                scheme,
                "\nValid options are scheme={:trapeze, :midpoint, :euler | :euler_explicit | :euler_forward, :euler_implicit | :euler_backward, :gauss_legendre_2, :gauss_legendre_3}\n",
                typeof(scheme),
            )
        end

        # lower and upper bounds for variables and constraints
        bounds = DOCPbounds(
            -Inf * ones(dim_NLP_variables),
            Inf * ones(dim_NLP_variables),
            zeros(dim_NLP_constraints),
            zeros(dim_NLP_constraints),
        )

        # call constructor with const fields
        docp = new{typeof(discretization),typeof(ocp)}(
            discretization, ocp, flags, dims, time, bounds, dim_NLP_variables, dim_NLP_constraints,
        )

        return docp
    end
end

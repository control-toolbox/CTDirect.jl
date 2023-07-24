# NB. to be moved up to CTBase later

# init struct (to be passed to solve)
# functions of time for state and control variables
# values for optimization variables (NB. free t0/tf in particular !)
# later for indirect methods add costate and structure information

# constructors
# - default: keep empty for method-specific handling
# - from constant vector (dim x + dim u + dim v)
# - from existing solution

mutable struct OptimalControlInit

    state_init
    control_init
    variable_init
    info

    function OptimalControlInit(; x_init=nothing, u_init=nothing, v_init=nothing)
        init = new()
        init.info          = :constant
        init.state_init    = t -> x_init
        init.control_init  = t -> u_init
        init.variable_init = v_init
        return init
    end

    function OptimalControlInit(sol::OptimalControlSolution)
        init = new()
        init.info          = :solution
        init.state_init    = t -> sol.state(t)  # scalar init if dim=1
                             # t -> sol.state(t)[1:sol.state_dimension] # remove possible additional state for Lagrange cost
        init.control_init  = t -> sol.control(t) # scalar init if dim=1
        init.variable_init = sol.variable # scalar init if dim=1
        return init
    end
end

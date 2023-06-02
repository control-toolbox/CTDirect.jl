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

    state_dimension
    control_dimension
    variable_dimension
    state_init
    control_init
    variable_init
    info

    function OptimalControlInit()
        init = new()
        init.info = :undefined
        return init
    end

    function OptimalControlInit(x_init, u_init, v_init=Real[])
        init = new()
        init.info = :constant
        init.state_dimension = length(x_init)
        init.control_dimension = length(u_init)
        init.variable_dimension = length(v_init)
        #println("Init dims x u v: ",init.state_dimension," ",init.control_dimension," ",init.variable_dimension)
        init.state_init = t -> x_init
        init.control_init = t -> u_init
        init.variable_init = v_init
        return init
    end

    function OptimalControlInit(sol::OptimalControlSolution)
        init = new()
        init.info = :solution
        init.state_dimension = sol.state_dimension
        init.control_dimension = sol.control_dimension
        init.variable_dimension = sol.variable_dimension
        init.state_init = t-> sol.state(t)[1:sol.state_dimension]
        init.control_init = sol.control
        init.variable_init = sol.infos[:variable]  #+++ need proper variable field in ocp solution !
        return init
    end
end
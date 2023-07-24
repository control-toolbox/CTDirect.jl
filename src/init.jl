# NB. to be moved up to CTBase later

# init struct (to be passed to solve)
# functions of time for state and control variables
# values for optimization variables (NB. free t0/tf in particular !)
# later for indirect methods add costate and structure information

# constructors
# - default: keep empty for method-specific handling
# - from constant vector or function handles
# - from existing solution
"""
$(TYPEDSIGNATURES)

Initialization of the optimal control problem solution for the NLP problem.

# Constructors:

- `OptimalControlInit()`: default initialization
- `OptimalControlInit(x_init, u_init, v_init)`: constant vector or function handles
- `OptimalControlInit(sol)`: from existing solution

# Examples

```julia-repl
julia> init = OptimalControlInit()
julia> init = OptimalControlInit(x_init=[0.1, 0.2], u_init=0.3)
julia> init = OptimalControlInit(x_init=[0.1, 0.2], u_init=0.3, v_init=0.5)
julia> init = OptimalControlInit(x_init=[0.1, 0.2], u_init=t->sin(t), v_init=0.5)
julia> init = OptimalControlInit(sol)
```

"""
mutable struct OptimalControlInit

    state_init::Function
    control_init::Function
    variable_init::Union{Nothing, ctVector}
    info::Symbol

    # constructor from constant vector or function handles
    function OptimalControlInit(; x_init::Union{Nothing, ctVector, Function}=nothing, 
                                  u_init::Union{Nothing, ctVector, Function}=nothing, 
                                  v_init::Union{Nothing, ctVector}=nothing)
        init = new()
        init.info          = :constant_or_function
        init.state_init    = (x_init isa Function) ? t -> x_init(t) : t -> x_init
        init.control_init  = (u_init isa Function) ? t -> u_init(t) : t -> u_init
        init.variable_init = v_init
        return init
    end

    # constructor from existing solution
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

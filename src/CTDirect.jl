module CTDirect

using DocStringExtensions
import ADNLPModels
import ExaModels
import CTModels
import CTSolvers, CTSolvers.Strategies, CTSolvers.Options
import SolverCore
import SparseArrays
# using CTBase
# using NLPModels

# ----------------------------------------------------------------------
# TYPES
# ----------------------------------------------------------------------
const AbstractModel = CTModels.AbstractModel

# ---------------------------------------------------------------------------
# Abstract discretizer type
# ---------------------------------------------------------------------------
"""
$(TYPEDEF)

Abstract type for optimal control problem discretization strategies.

This type defines the interface for converting continuous optimal control problems
into discrete optimization problems suitable for numerical solvers. Subtypes must
implement the callable interface to perform the actual discretization.

# Interface Requirements

Subtypes must implement:
- `(discretizer::AbstractDiscretizer)(ocp::AbstractModel)`: Perform discretization

# Example
```julia-repl
julia> using CTDirect

julia> MyDiscretizer <: AbstractDiscretizer end

julia> # Implement the required interface
julia> (d::MyDiscretizer)(ocp) = discretize_collocation(ocp)
```

See also: `Collocation`, `discretize`
"""
abstract type AbstractDiscretizer <: Strategies.AbstractStrategy end

"""
$(TYPEDSIGNATURES)

Discretize an optimal control problem using the specified discretizer.

This function applies a discretization strategy to convert a continuous optimal
control problem into a form suitable for numerical optimization.

# Arguments
- `ocp::AbstractModel`: The optimal control problem to discretize
- `discretizer::AbstractDiscretizer`: The discretization strategy to apply

# Returns
- The discretized problem (type depends on the specific discretizer)

# Example
```julia-repl
julia> using CTDirect

julia> ocp = create_ocp()  # Create your OCP
julia> discretized = discretize(ocp, Collocation())
```

See also: `AbstractDiscretizer`, `Collocation`
"""
function discretize(
    ocp::AbstractModel, 
    discretizer::AbstractDiscretizer
)
    return discretizer(ocp)
end

"""
$(TYPEDSIGNATURES)

Returns the default discretizer instance used by the convenience method.

This function provides the default collocation discretizer that is used when
no specific discretizer is provided to the `discretize` function.

# Returns
- `AbstractDiscretizer`: The default collocation discretizer instance

# Notes
- This function is internal and primarily used for providing default values
- The default is currently `Collocation()` but may change in future versions

See also: `discretize`, `Collocation`
"""
__discretizer()::AbstractDiscretizer = Collocation()

"""
$(TYPEDSIGNATURES)

Discretize an optimal control problem using the default discretizer.

This is a convenience method that uses the default collocation discretizer
when no specific discretizer is provided.

# Arguments
- `ocp::AbstractModel`: The optimal control problem to discretize
- `discretizer::AbstractDiscretizer`: The discretization strategy (default: Collocation())

# Returns
- The discretized problem using the specified or default discretizer

# Example
```julia-repl
julia> using CTDirect

julia> ocp = create_ocp()  # Create your OCP
julia> discretized = discretize(ocp)  # Uses default Collocation
julia> discretized_custom = discretize(ocp, MyCustomDiscretizer())
```

# Notes
- The default discretizer is `Collocation()`
- For custom discretization, provide a specific discretizer instance

See also: `AbstractDiscretizer`, `Collocation`
"""
function discretize(
    ocp::AbstractModel;
    discretizer::AbstractDiscretizer=__discretizer(),
)
    return discretize(ocp, discretizer)
end

# ---------------------------------------------------------------------------
# Discretization schemes: see disc/
# ---------------------------------------------------------------------------
"""
$(TYPEDEF)

Abstract type representing a discretization strategy for an optimal
control problem.  

Concrete subtypes of `Discretization` define specific schemes for
transforming a continuous-time problem into a discrete-time
representation suitable for numerical solution.

# Example

```julia-repl
julia> struct MyDiscretization <: Discretization end
MyDiscretization
```
"""
abstract type Discretization end


# includes
include("collocation.jl")
include("collocation_core.jl")
include("collocation_variables.jl")
include("collocation_functions.jl")
include("disc/common.jl")
include("disc/euler.jl")
include("disc/irk.jl")
include("disc/midpoint.jl")
include("disc/trapeze.jl")

end
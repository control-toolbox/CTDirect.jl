module CTDirect

using DocStringExtensions
import ADNLPModels
import ExaModels
import CTModels
import CTSolvers, CTSolvers.Strategies, CTSolvers.Options
import SolverCore
import SparseArrays

# ---------------------------------------------------------------------------
# Abstract discretizer type
# ---------------------------------------------------------------------------
const AbstractModel = CTModels.AbstractModel
abstract type AbstractDiscretizer <: Strategies.AbstractStrategy end

__discretizer()::AbstractDiscretizer = Collocation()

"""
$(TYPEDSIGNATURES)

Discretize an optimal control problem using the specified discretizer.

# Arguments
- `ocp::AbstractModel`: The optimal control problem to discretize
- `discretizer::AbstractDiscretizer`: The discretization strategy to apply

# Returns
- The discretized problem representation
"""
function discretize(ocp::AbstractModel, discretizer::AbstractDiscretizer)
    return discretizer(ocp)
end
"""
$(TYPEDSIGNATURES)

Discretize an optimal control problem using the default discretizer.

This is a convenience method that uses the default discretizer (Collocation).

# Arguments
- `ocp::AbstractModel`: The optimal control problem to discretize
- `discretizer::AbstractDiscretizer`: Optional discretization strategy (default: Collocation)

# Returns
- The discretized problem representation
"""
function discretize(ocp::AbstractModel; discretizer::AbstractDiscretizer=__discretizer())
    return discretize(ocp, discretizer)
end

# ---------------------------------------------------------------------------
# Discretization schemes: see ode/
# ---------------------------------------------------------------------------
"""
$(TYPEDEF)

Abstract type representing a discretization scheme strategy for an optimal
control problem.  

Concrete subtypes of `Scheme` define specific schemes for
transforming a continuous-time problem into a discrete-time
representation suitable for numerical solution.

# Example

```julia-repl
julia> struct MyScheme <: Scheme end
MyScheme
```
"""
abstract type Scheme end


# includes
include("DOCP_data.jl")
include("DOCP_variables.jl")
include("DOCP_functions.jl")

include("ode/common.jl")
include("ode/euler.jl")
include("ode/irk.jl")
include("ode/irk_stagewise.jl")
include("ode/midpoint.jl")
include("ode/trapeze.jl")

#include("ode/variable.jl")

include("collocation.jl")
include("direct_shooting.jl")

end

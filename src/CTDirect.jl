module CTDirect

using ADNLPModels
using ExaModels
using CTBase
using CTModels
import CTSolvers, CTSolvers.Strategies, CTSolvers.Options
using DocStringExtensions
using SparseArrays
using SolverCore
using NLPModels

# ----------------------------------------------------------------------
# TYPES
const AbstractModel = CTModels.AbstractModel

# ---------------------------------------------------------------------------
# Abstract discretizer type
# ---------------------------------------------------------------------------
abstract type AbstractDiscretizer <: Strategies.AbstractStrategy end

function discretize(
    ocp::AbstractModel, 
    discretizer::AbstractDiscretizer
)
    return discretizer(ocp)
end

__discretizer()::AbstractDiscretizer = Collocation()

function discretize(
    ocp::AbstractModel;
    discretizer::AbstractDiscretizer=__discretizer(),
)
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
include("ode/common.jl")
include("ode/euler.jl")
include("ode/irk.jl")
include("ode/midpoint.jl")
include("ode/trapeze.jl")

include("collocation.jl")
include("collocation_variables.jl")
include("collocation_functions.jl")

# ode/variable_step
#include("direct_shooting.jl")
#include("direct_shooting_core.jl")

end

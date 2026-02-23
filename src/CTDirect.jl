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
include("ode/common.jl")
include("ode/euler.jl")
include("ode/irk.jl")
include("ode/midpoint.jl")
include("ode/trapeze.jl")
# ode/variable_step

end

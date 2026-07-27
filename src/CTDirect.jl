module CTDirect

using DocStringExtensions
import ADNLPModels
import ExaModels
import CTModels
import CTSolvers
import CTBase
import CTBase.Strategies
import CTBase.Options
import CTBase.Core
import CTBase.Exceptions
import SolverCore
import SparseArrays

# ---------------------------------------------------------------------------
# Discretizers
#
# `AbstractDiscretizer` is owned by CTSolvers (CTSolvers.DOCP). Concrete
# discretizers (Collocation, DirectShooting) implement the CTSolvers contract:
# `CTSolvers.discretize`, `CTSolvers.build_model`, `CTSolvers.build_solution`.
# ---------------------------------------------------------------------------
const AbstractModel = CTModels.AbstractModel

__discretizer()::CTSolvers.AbstractDiscretizer = Collocation()

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
include("DOCP_cache.jl")
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

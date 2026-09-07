"""
Direct transcription methods for optimal control problems.

CTDirect turns a continuous OCP into a nonlinear program by discretizing the state and
control on a time grid — `CTDirect.Collocation` (collocation schemes, `ADNLP` and
`Exa` backends) or `CTDirect.DirectShooting` (sequential shooting, `ADNLP` only) —
and implements the CTSolvers `discretize` / `build_model` / `build_solution` contract so
the resulting NLP can be handed to a CTSolvers modeler and solver.
"""
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
"""
Alias for `CTModels.AbstractModel`, the continuous-time OCP model type consumed by the
discretizers.
"""
const AbstractModel = CTModels.AbstractModel

"""
$(TYPEDSIGNATURES)

Default discretizer used when none is given: `CTDirect.Collocation()`.
"""
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

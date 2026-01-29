module CTDirect

using ADNLPModels
using ExaModels
using CTBase
using CTModels: CTModels
using DocStringExtensions
using SparseArrays
using SolverCore: SolverCore
using NLPModels: NLPModels

# ----------------------------------------------------------------------
# TYPES
const AbstractOptimalControlProblem = CTModels.AbstractModel

# includes
include("core_types.jl")
include("discretization_api.jl")
include("collocation.jl")
include("collocation_core.jl")
include("disc/common.jl")
include("disc/euler.jl")
include("disc/irk.jl")
include("disc/midpoint.jl")
include("disc/trapeze.jl")

end

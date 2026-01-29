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

# # ----------------------------------------------------------------------
# # EXTENSIONS

# # NLP model backend extensions
#abstract type AbstractNLPModelBackend end
#struct ADNLPBackend <: AbstractNLPModelBackend end
#struct ExaBackend <: AbstractNLPModelBackend end

# includes
include("utils.jl")
include("core_types.jl")
include("discretization_api.jl")
include("collocation.jl")
include("collocation_core.jl")

include("disc/common.jl")
include("disc/euler.jl")
include("disc/irk.jl")
include("disc/midpoint.jl")
include("disc/trapeze.jl")

include("solution.jl")


end
